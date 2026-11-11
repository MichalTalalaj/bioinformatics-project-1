"""Write a Markdown gene report using NCBI, UniProtKB and Ensembl."""

import argparse
import re
import time
from pathlib import Path
from urllib.parse import quote

import requests


EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
ENSEMBL_URL = "https://rest.ensembl.org/lookup/symbol/"


# endpoint specifies the NCBI E-utilities function (esearch, elink)
# parameters is a dictionary of query parameters for the request (db, term, rettype, etc.)
# returns a dictionary parsed from the JSON response
def get_ncbi_json(endpoint, parameters):
    """Request JSON and reject HTTP failures or NCBI error responses."""
    # Space consecutive NCBI calls to stay below three requests per second.
    time.sleep(0.34)
    response = requests.get(EUTILS_URL + endpoint, params=parameters, timeout=30)
    response.raise_for_status()
    data = response.json()
    if not isinstance(data, dict) or data.get("error"):
        raise ValueError("NCBI returned an error or an unexpected response.")
    return data



# takes symbol (e.g. APOE) and organism (e.g. Homo sapiens) as input
# returns the NCBI Gene ID
def find_gene_id(symbol, organism):
    """Search NCBI Gene and require exactly one result."""
    search = get_ncbi_json("esearch.fcgi", {
        "db": "gene",
        "term": f'"{symbol}"[Gene Name] AND "{organism}"[Organism]',
        "retmode": "json",
        "retmax": 2,
    })
    result = search.get("esearchresult", {})
    if result.get("errorlist"):
        raise ValueError("NCBI could not interpret the gene or organism query.")
    if "count" not in result or "idlist" not in result:
        raise ValueError("NCBI search response is missing its count or Gene IDs.")
    count = int(result["count"])
    if count == 0:
        raise ValueError("No matching NCBI Gene record was found.")
    if count != 1:
        raise ValueError("Multiple NCBI Gene records matched; use a specific gene symbol and organism.")
    if len(result["idlist"]) != 1:
        raise ValueError("NCBI search returned an inconsistent list of Gene IDs.")
    return result["idlist"][0]

# simillar to find_gene_id, but fetches the gene summary and checks that it matches the requested symbol and organism
# takes gene_id (from find_gene_id), symbol, and organism as input
def fetch_gene_summary(gene_id, symbol, organism):
    """Retrieve a summary and check that it identifies the requested gene."""
    summary = get_ncbi_json("esummary.fcgi", {
        "db": "gene", "id": gene_id, "retmode": "json",
    })
    gene = summary.get("result", {}).get(gene_id)
    if not isinstance(gene, dict) or gene.get("error"):
        raise ValueError("NCBI did not return a usable gene summary.")

    returned_symbol = gene.get("nomenclaturesymbol") or gene.get("name", "")
    returned_organism = gene.get("organism", {}).get("scientificname", "")
    # Gene searches may match aliases or related taxa; do not silently accept those.
    if returned_symbol.casefold() != symbol.casefold():
        raise ValueError("The returned gene symbol does not match the requested symbol.")
    if returned_organism.casefold() != organism.casefold():
        raise ValueError("The returned organism does not match the requested scientific name.")
    if gene.get("uid") != gene_id:
        raise ValueError("The returned Gene ID does not match the search result.")

    return gene


# based on fetch_gene_summary returns dictionary with gene information
# required by the task pdf
# takes paresed user input (symbol and organism)
def get_ncbi_gene(symbol, organism):
    """Find the gene, retrieve its summary and extract readable report fields."""
    gene_id = find_gene_id(symbol, organism)
    gene = fetch_gene_summary(gene_id, symbol, organism)

    mim_references = []
    for mim_id in gene.get("mim") or []:
        mim_references.append(str(mim_id))

    return {
        "Gene ID": gene_id,
        "Symbol": gene.get("nomenclaturesymbol") or gene.get("name", ""),
        "Organism": gene.get("organism", {}).get("scientificname", ""),
        "Full name": gene.get("nomenclaturename") or gene.get("description") or "Unavailable",
        "Chromosome": gene.get("chromosome") or "Unavailable",
        "Map location": gene.get("maplocation") or "Unavailable",
        "Aliases": gene.get("otheraliases") or "Unavailable",
        "Summary": gene.get("summary") or "Unavailable",
        "MIM references": ", ".join(mim_references) or "Unavailable",
    }

# based on gene_id finds all RefSeq RNA records linked to this Gene ID
def find_refseq_rna_ids(gene_id):
    """Find all RefSeq RNA records linked to this Gene ID."""
    data = get_ncbi_json("elink.fcgi", {
        "dbfrom": "gene", "db": "nuccore", "id": gene_id,
        "linkname": "gene_nuccore_refseqrna", "retmode": "json",
    })
    if "linksets" not in data:
        raise ValueError("NCBI link response is missing its link sets.")

    rna_ids = []
    for linkset in data["linksets"]:
        if linkset.get("error"):
            raise ValueError("NCBI could not retrieve linked RefSeq RNA records.")
        for database in linkset.get("linksetdbs", []):
            if database.get("linkname") == "gene_nuccore_refseqrna":
                for rna_id in database.get("links", []):
                    if str(rna_id) not in rna_ids:
                        rna_ids.append(str(rna_id))
    return rna_ids


# extracts CDS coordinates, protein_id, and transl_table from
# a single RefSeq mRNA record in GenBank format (which is the input)
def parse_cds_features(record):
    """Read CDS locations and the two required qualifiers from GenBank text."""
    coding_sequences = []
    current_cds = None
    reading_location = False

    for line in record.splitlines():
        if line.startswith("ORIGIN"):
            break
        if not line.startswith("     "):
            continue

        feature_name = line[5:21].strip()
        value = line[21:].strip()
        if feature_name:
            current_cds = None
            reading_location = False
            if feature_name == "CDS":
                current_cds = {
                    "Location": value,
                    "protein_id": "Unavailable",
                    "transl_table": "Unavailable",
                }
                coding_sequences.append(current_cds)
                reading_location = True
        elif current_cds is not None:
            if value.startswith("/"):
                reading_location = False
                qualifier, separator, content = value[1:].partition("=")
                if separator and qualifier in ("protein_id", "transl_table"):
                    current_cds[qualifier] = content.strip('"') or "Unavailable"
            elif reading_location:
                # Keep wrapped join/complement locations in GenBank notation.
                current_cds["Location"] += value

    return coding_sequences

# simillar to parse_cds_features
# extracts RefSeq mRNA accession and version, length, last modification date
# asserts that the record is correctly linked to the requested Gene ID and organism

# returns both the transcript information and the CDS information

def parse_refseq_mrna(record, gene_id, organism):
    """Extract required fields from one mRNA record; skip noncoding RNAs."""
    record = record.replace("\r\n", "\n")
    locus = re.search(r"^LOCUS\s+(.+)$", record, re.MULTILINE)
    version = re.search(r"^VERSION\s+(\S+)", record, re.MULTILINE)
    if locus is None or version is None:
        raise ValueError("NCBI returned a record without a LOCUS or VERSION header.")

    locus_fields = locus.group(1).split()
    accession = version.group(1)
    if "mRNA" not in locus_fields or not accession.startswith(("NM_", "XM_")):
        return None
    if f'/db_xref="GeneID:{gene_id}"' not in record:
        raise ValueError("A linked mRNA record does not identify the requested Gene ID.")
    source_organism = re.search(r"^  ORGANISM\s+(.+)$", record, re.MULTILINE)
    if source_organism is None or source_organism.group(1).strip().casefold() != organism.casefold():
        raise ValueError("A linked mRNA record does not match the requested organism.")

    modified = re.search(r"\b\d{2}-[A-Z]{3}-\d{4}$", locus.group(1))
    return {
        "Accession/version": accession,
        "Length (nt)": int(locus_fields[1]),
        "Last modified": modified.group(0) if modified else "Unavailable",
        "CDS": parse_cds_features(record),
    }

# combines 3 functions: find_refseq_rna_ids, parse_refseq_mrna, and get_ncbi_json
# gets all linked RefSeq mRNAs for a given gene_id and organism
# returns a list of dictionaries with transcript information and CDS information

def get_ncbi_nucleotide(gene_id, organism):
    """Retrieve all linked RefSeq mRNAs, without choosing one transcript."""
    rna_ids = find_refseq_rna_ids(gene_id)
    if not rna_ids:
        return []

    time.sleep(0.34)
    response = requests.get(EUTILS_URL + "efetch.fcgi", params={
        "db": "nuccore", "id": ",".join(rna_ids),
        "rettype": "gb", "retmode": "text",
    }, timeout=30)
    response.raise_for_status()
    if not response.text.strip().endswith("//"):
        raise ValueError("NCBI returned an incomplete GenBank response.")
    records = re.split(r"(?m)^//\s*$", response.text)
    transcripts = []
    record_count = 0
    for record in records:
        if not record.strip():
            continue
        record_count += 1
        transcript = parse_refseq_mrna(record, gene_id, organism)
        if transcript is not None:
            for previous in transcripts:
                if previous["Accession/version"] == transcript["Accession/version"]:
                    raise ValueError("NCBI returned a duplicate mRNA record.")
            transcripts.append(transcript)
    if record_count != len(rna_ids):
        raise ValueError("NCBI did not return all linked RefSeq RNA records.")
    transcripts.sort(key=lambda transcript: transcript["Accession/version"])
    return transcripts



# collect all protein entries from uniprot matchin the symbol and organism
# responses are limited to either reviewed or unreviewed entries
# returns a list of dictionaries with protein information and annotations

def search_uniprot_entries(symbol, organism, reviewed):
    """Retrieve all exact gene/organism matches for one review status."""
    status = "true" if reviewed else "false"
    parameters = {
        "query": f'gene_exact:"{symbol}" AND organism_name:"{organism}" AND reviewed:{status}',
        "format": "json",
        "fields": (
            "accession,reviewed,protein_name,gene_names,organism_name,organism_id,"
            "cc_function,length,ft_binding,cc_cofactor,cc_disease,xref_pdb,keyword"
        ),
        "size": 500,
    }
    url = UNIPROT_URL
    entries = []
    while url:
        response = requests.get(url, params=parameters, timeout=30)
        response.raise_for_status()
        data = response.json()
        if not isinstance(data, dict) or not isinstance(data.get("results"), list):
            raise ValueError("UniProt returned an unexpected search response.")

        for entry in data["results"]:
            returned_organism = entry.get("organism", {}).get("scientificName", "")
            gene_matches = False
            for gene in entry.get("genes", []):
                returned_symbol = gene.get("geneName", {}).get("value", "")
                if returned_symbol.casefold() == symbol.casefold():
                    gene_matches = True
            if gene_matches and returned_organism.casefold() == organism.casefold():
                if not entry.get("primaryAccession"):
                    raise ValueError("UniProt returned a matching entry without an accession.")
                entries.append(entry)

        url = response.links.get("next", {}).get("url")
        # The next-page URL already contains the original query and selected fields.
        parameters = None
    return entries

# helper function to extract the text values from dictionaries
# returns a list of strings

def uniprot_texts(items):
    """Read the text values used in UniProt comments and notes."""
    texts = []
    for item in items:
        if item.get("value"):
            texts.append(item["value"])
    return texts


# devides the comments into three categories: Function, Cofactors, and Diseases
# returns a dictionary with three lists of strings

def parse_uniprot_comments(entry):
    """Keep function descriptions, cofactor annotations and disease context."""
    annotations = {"Function": [], "Cofactors": [], "Diseases": []}
    for comment in entry.get("comments", []):
        comment_type = comment.get("commentType")
        if comment_type == "FUNCTION":
            annotations["Function"].extend(uniprot_texts(comment.get("texts", [])))
        elif comment_type == "COFACTOR":
            for cofactor in comment.get("cofactors", []):
                if cofactor.get("name"):
                    annotations["Cofactors"].append(cofactor["name"])
            notes = comment.get("note", {}).get("texts", [])
            annotations["Cofactors"].extend(uniprot_texts(notes))
        elif comment_type == "DISEASE":
            disease = comment.get("disease", {})
            description = []
            for field in ("diseaseId", "diseaseAccession", "description"):
                if disease.get(field):
                    description.append(disease[field])
            # Preserve notes distinguishing susceptibility from disease causation.
            notes = comment.get("note", {}).get("texts", [])
            description.extend(uniprot_texts(notes))
            if description:
                annotations["Diseases"].append(" | ".join(description))
    return annotations

# based on previous functions, extracts the required protein fields f0m one entry
# returns a dictionary with protein information and annotations

def parse_uniprot_entry(entry):
    """Extract the required protein fields from one verified UniProt entry."""
    description = entry.get("proteinDescription", {})
    protein_name = description.get("recommendedName", {}).get("fullName", {}).get("value")
    if not protein_name:
        for submitted_name in description.get("submissionNames", []):
            protein_name = submitted_name.get("fullName", {}).get("value")
            if protein_name:
                break

    ligands = []
    for feature in entry.get("features", []):
        if feature.get("type") == "Binding site":
            name = feature.get("ligand", {}).get("name")
            if name and name not in ligands:
                ligands.append(name)

    structures = []
    for reference in entry.get("uniProtKBCrossReferences", []):
        if reference.get("database") == "PDB" and reference.get("id"):
            structures.append(reference["id"])

    keywords = []
    for keyword in entry.get("keywords", []):
        if keyword.get("name"):
            keywords.append(keyword["name"])

    comments = parse_uniprot_comments(entry)
    return {
        "Accession": entry["primaryAccession"],
        "Review status": entry.get("entryType") or "Unavailable",
        "Protein name": protein_name or "Unavailable",
        "Function": "\n".join(comments["Function"]) or "Unavailable",
        "Length (aa)": entry.get("sequence", {}).get("length") or "Unavailable",
        "Ligands": ", ".join(ligands) or "Unavailable",
        "Cofactors": "; ".join(comments["Cofactors"]) or "Unavailable",
        "Diseases": "\n".join(comments["Diseases"]) or "Unavailable",
        "PDB structures": ", ".join(structures) or "Unavailable",
        "Keywords": ", ".join(keywords) or "Unavailable",
    }


# searches for all UniProt entries matching the gene symbol and organism
# returns a list of dictionaries with protein information and annotations
# in case there are no reviewed entries, it will return unreviewed entries instead

def get_uniprot_proteins(symbol, organism):
    """Prefer all reviewed matches; use unreviewed matches only if none exist."""
    entries = search_uniprot_entries(symbol, organism, reviewed=True)
    if not entries:
        entries = search_uniprot_entries(symbol, organism, reviewed=False)
    entries.sort(key=lambda entry: entry["primaryAccession"])

    proteins = []
    for entry in entries:
        proteins.append(parse_uniprot_entry(entry))
    return proteins


# makes a call to the Ensembl REST API 
# fetches the gene record for the given symbol and organism
# checks that the returned record matches the requested symbol and organism
# returns the gene record as a dictionary

def fetch_ensembl_gene(symbol, organism):
    """Retrieve the expanded gene record and check its symbol and species."""
    species = "_".join(organism.lower().split())
    url = ENSEMBL_URL + quote(species, safe="") + "/" + quote(symbol, safe="")
    response = requests.get(url, params={
        "expand": 1, "content-type": "application/json",
    }, timeout=(10, 90))
    response.raise_for_status()
    gene = response.json()
    if not isinstance(gene, dict) or gene.get("error"):
        raise ValueError("Ensembl returned an error or an unexpected response.")
    if gene.get("object_type") != "Gene":
        raise ValueError("Ensembl did not return a gene record.")
    if gene.get("display_name", "").casefold() != symbol.casefold():
        raise ValueError("The Ensembl gene symbol does not match the requested symbol.")
    if gene.get("species", "").casefold() != species.casefold():
        raise ValueError("The Ensembl species does not match the requested organism.")
    return gene

# parses the Ensembl gene record and extracts the required fields
# counts the number of transcripts returned in the expanded record
# checks that the canonical transcript is present in the transcript list (if both canonical and transcript list are present)

def parse_ensembl_gene(gene):
    """Extract gene fields and count the returned expanded transcripts."""
    transcripts = gene.get("Transcript")
    transcript_count = "Unavailable"
    if transcripts is not None:
        if not isinstance(transcripts, list):
            raise ValueError("Ensembl returned an unexpected transcript list.")
        transcript_count = len(transcripts)

    canonical = gene.get("canonical_transcript")
    if canonical and transcripts is not None:
        # Transcript IDs omit the version suffix used by canonical_transcript.
        canonical_id = canonical.split(".")[0]
        canonical_found = False
        for transcript in transcripts:
            if transcript.get("id") == canonical_id:
                canonical_found = True
        if not canonical_found:
            raise ValueError("The canonical transcript is absent from Ensembl's transcript list.")

    return {
        "Gene ID": gene.get("id") or "Unavailable",
        "Assembly": gene.get("assembly_name") or "Unavailable",
        "Sequence region": gene.get("seq_region_name") or "Unavailable",
        "Start (1-based, inclusive)": gene.get("start") or "Unavailable",
        "End (1-based, inclusive)": gene.get("end") or "Unavailable",
        "Strand (+1 forward, -1 reverse)": gene.get("strand") or "Unavailable",
        "Transcript count": transcript_count,
        "Canonical transcript": canonical or "Unavailable",
        "Biotype": gene.get("biotype") or "Unavailable",
    }


##############################

# this section is responsible for formatting the report and writing it to a file
# no bioinformatics is being done here, just formatting and writing to a file

def markdown_text(value):
    """Keep annotation punctuation literal when displayed as Markdown."""
    text = str(value)
    for character in ("\\", "`", "*", "_", "[", "]", "<", ">"):
        text = text.replace(character, "\\" + character)
    return text


def format_markdown_fields(fields):
    """Display fields as labelled paragraphs, keeping separate source paragraphs."""
    paragraphs = []
    for field, value in fields.items():
        text = markdown_text(value).replace("\n", "\n\n")
        paragraphs.append(f"**{field}:** {text}")
    return "\n\n".join(paragraphs)


def build_report(symbol, organism, gene, transcripts, proteins, ensembl_fields):
    """Assemble the four database sections without dropping any extracted fields."""
    sections = [
        f"# Gene report: {markdown_text(symbol)}",
        f"**Organism:** {markdown_text(organism)}",
        "Unavailable means the requested annotation was not returned; it does not prove biological absence.",
        "## NCBI Gene",
        format_markdown_fields(gene),
        "## NCBI Nucleotide",
        "All linked RefSeq mRNA variants are included. CDS locations use 1-based inclusive transcript "
        "coordinates; GenBank location notation is preserved.",
    ]
    if not transcripts:
        sections.append("**RefSeq mRNA:** Unavailable")
    for transcript in transcripts:
        sections.append(f"### {markdown_text(transcript['Accession/version'])}")
        fields = {
            "Accession/version": transcript["Accession/version"],
            "Length (nt)": transcript["Length (nt)"],
            "Last modified": transcript["Last modified"],
        }
        sections.append(format_markdown_fields(fields))
        if not transcript["CDS"]:
            sections.append("**CDS:** Unavailable")
        for number, cds in enumerate(transcript["CDS"], start=1):
            sections.append(f"#### CDS {number}")
            sections.append(format_markdown_fields(cds))

    sections.append("## UniProtKB")
    sections.append("All exact reviewed matches are preferred; unreviewed matches are included only "
                    "when no exact reviewed match is available.")
    if not proteins:
        sections.append("**UniProtKB entry:** Unavailable")
    for protein in proteins:
        sections.append(f"### {markdown_text(protein['Accession'])}")
        sections.append(format_markdown_fields(protein))

    sections.append("## Ensembl")
    sections.append(format_markdown_fields(ensembl_fields))
    return "\n\n".join(sections) + "\n"


def write_report(output_directory, symbol, organism, report):
    """Create the output directory and write or replace one UTF-8 Markdown file."""
    filename = re.sub(r"[^A-Za-z0-9_-]+", "_", f"{symbol}_{organism}") + ".md"
    output_directory.mkdir(parents=True, exist_ok=True)
    report_path = output_directory / filename
    report_path.write_text(report, encoding="utf-8")
    return report_path


##############################

# parser is added to interpret the command line arguments
# then the main function calls all the other functions in order to generate the report

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("symbol", help="Official gene symbol, for example APOE.")
    parser.add_argument("--organism", required=True, help="Scientific name, for example Homo sapiens.")
    parser.add_argument(
        "--out", required=True, type=Path,
        help="Output directory; the gene/organism Markdown report is replaced on reruns.",
    )
    args = parser.parse_args()
    symbol = args.symbol.strip()
    organism = args.organism.strip()
    if not symbol or not organism:
        parser.error("Gene symbol and organism must not be empty.")

    try:
        gene = get_ncbi_gene(symbol, organism)
        transcripts = get_ncbi_nucleotide(gene["Gene ID"], organism)
        proteins = get_uniprot_proteins(symbol, organism)
        ensembl_gene = fetch_ensembl_gene(symbol, organism)
        ensembl_fields = parse_ensembl_gene(ensembl_gene)
        report = build_report(symbol, organism, gene, transcripts, proteins, ensembl_fields)
        report_path = write_report(args.out, symbol, organism, report)
    except requests.RequestException as error:
        parser.exit(2, f"Error: Database request failed ({type(error).__name__}).\n")
    except (ValueError, OSError) as error:
        parser.exit(2, f"Error: {error}\n")

    print(f"Wrote report to {report_path}.")


if __name__ == "__main__":
    main()
