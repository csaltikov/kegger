import shutil
import tempfile
import requests
import requests_cache
import pandas as pd
import io
import os
from collections import defaultdict
from datetime import timedelta


def initialize_kegger(cache_path: str | None = None, expire_days: int = 30):
    """
        Sets up a persistent local cache for KEGG API requests.

        This function initializes a SQLite database to store API responses. Subsequent
        calls to the same KEGG URL will pull data from the local cache instead of
        the internet, significantly speeding up data processing and reducing
        server load.

        Args:
            cache_path (str | None): The filename or path for the SQLite cache.
                Defaults to "kegg_cache" (which creates 'kegg_cache.sqlite').
            expire_days (int): How many days a cached response remains valid
                before a fresh request is forced. Defaults to 30.

        Note:
            If a cache file already exists at the specified path, this function
            will automatically load and reuse it.
    """
    if cache_path is None:
        cache_path = "kegg_cache"
    requests_cache.install_cache(cache_path,
                                 backend="sqlite",
                                 expire_days=timedelta(days=expire_days))


def get_url(url: str, request_type="text") -> str | dict:
    response = requests.get(url)
    response.raise_for_status()
    if request_type == "json":
        result = response.json()
    else:
        result = response.text
    return result


def get_kegg(results, org: str, save=False) -> pd.DataFrame:
    cols = ["pathid", "description"]
    pathways_df = pd.read_csv(results, sep="\t", header=None, names=cols)
    if save:
        pathways_df.to_csv(f"{org}_pathways.tsv", sep="\t", index=False)
    return pathways_df


def list_all_pathways(org: str) -> pd.DataFrame:
    """
        Retrieves a list of all KEGG pathways for a specific organism.

        This function queries the KEGG REST API to find every metabolic and signaling
        pathway associated with the provided organism code. It automatically cleans
        the 'path:' prefix from the results to simplify downstream data merging.

        Args:
            org (str): The 3-4 letter KEGG organism code (e.g., 'eco' for E. coli,
                'hsa' for humans, or 'mmu' for mouse).

        Returns:
            pd.DataFrame: A DataFrame with two columns:
                - 'pathid': The unique KEGG pathway identifier (e.g., 'eco00010').
                - 'description': The human-readable name of the pathway.

        Example:
            >>> import kegger
            >>> pathways = kegger.list_all_pathways('eco')
            >>> print(pathways.head())
    """
    url = f'https://rest.kegg.jp/list/pathway/{org}'
    response = get_url(url)
    record = io.StringIO(response)
    col_names = ["pathid", "description"]
    df = pd.read_csv(record, sep="\t", header=None, names=col_names)
    df["pathid"] = df["pathid"].str.replace("path:", "", regex=False)
    return df


def genes_to_pathways(org: str) -> pd.DataFrame:
    """
    Retrieves the mapping between genes and their associated pathways for an organism.

    Queries the KEGG 'link' endpoint to produce a many-to-many map of pathways
    and gene identifiers. This is useful for enrichment analysis or finding
    all genes within a specific biological process.

    Args:
        org: The KEGG organism code (e.g., 'shn' or 'eco').

    Returns:
        pd.DataFrame: A two-column DataFrame:
            - 'pathid': The KEGG pathway identifier (e.g., 'shn00010').
            - 'gene': The specific gene identifier (e.g., 'Shewana3_0001').

    Example:
        >>> df_map = genes_to_pathways('shn')
        >>> # To find all genes in a specific pathway:
        >>> glycolysis = df_map[df_map['pathid'] == 'shn00010']
    """
    url = f'https://rest.kegg.jp/link/{org}/pathway'
    response = get_url(url)
    record = io.StringIO(response)
    col_names = ["pathid", "gene"]
    df = pd.read_csv(record, sep="\t", header=None, names=col_names)
    # Removing the 'org' and 'path:' prefixes
    df["pathid"] = df["pathid"].str.replace("path:", "", regex=False)
    df["gene"] = df["gene"].str.replace(f"{org}:", "", regex=False)
    return df


def get_module(mdid: str):
    url = f"https://rest.kegg.jp/get/md:{mdid}"
    request = get_url(url, request_type="text")
    return request


def get_path(pathid: str):
    url = f"https://rest.kegg.jp/get/{pathid}"
    request = get_url(url, request_type="text")
    return request


def get_entry(entry_id: str):
    url = f"https://rest.kegg.jp/get/{entry_id}"
    request = get_url(url, request_type="text")
    return request


def get_org(org: str) -> pd.DataFrame:
    """
    Retrieve and parse a KEGG organism genome list.

    Connects to the KEGG REST API 'list' endpoint to retrieve gene-level
    metadata and converts the tab-delimited response into a structured DataFrame.

    Parameters
    ----------
    org : str
        The three- or four-letter KEGG organism identifier.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with the following columns:
        - gene: KEGG gene identifier (e.g., 'Shewana3_0001')
        - feature: Biological category (e.g., 'CDS', 'RNA')
        - position: Chromosomal coordinates
        - annotation: Functional description/gene name
    """
    url = f"https://rest.kegg.jp/list/{org}"
    response = requests.get(url)
    record = io.StringIO(response.text)
    cols = ["gene", "feature", "position", "annotation"]
    df = pd.read_csv(record, sep="\t", header=None, names=cols)
    # Removing the 'org' prefix
    df["gene"] = df["gene"].str.replace(f"{org}:", "", regex=False)
    return df


def clean_entry(entry: dict) -> dict:
    """
        Standardizes and cleans raw KEGG record fields.

        This internal helper takes a dictionary of raw tags and values and applies
        specific parsing logic based on the KEGG field type (e.g., splitting GENE
        identifiers from their descriptions or parsing PATHWAY_MAP strings).

        Args:
            entry (dict): A dictionary where keys are KEGG tags (e.g., 'GENE')
                and values are lists of raw string lines.

        Returns:
            dict: A cleaned dictionary where values are strings or lists
                of strings, structured for easier data analysis.

        Note:
            Special handling is applied to 'GENE' fields to separate gene IDs
            from their associated 'ORTHOLOG' identifiers.
    """
    cleaned_entry = defaultdict(list)
    for tag, value in entry.items():
        if tag == "ENTRY":
            cleaned_entry[tag] = value[0].split()
        elif tag in ("NAME", "ORGANISM"):
            cleaned_entry[tag] = value[0].strip()
        elif tag == "GENE":
            cleaned_entry[tag], cleaned_entry["ORTHOLOG"] = map(list, zip(*[v.split(None, 1) for v in value]))
        elif tag in ["REL_PATHWAY"]:
            cleaned_entry[tag] = [v.strip() for v in value if v.strip()]
        elif tag == "PATHWAY_MAP":
            cleaned_entry[tag] = value[0].split("  ")
        elif tag in ("PATHWAY", "GENES", "REACTION", "MODULE"):
            for v in value:
                cleaned_entry[tag].append(v.strip())
        else:
            cleaned_entry[tag] = value[0].strip()

    return cleaned_entry


def kegg_parser(request_text: str) -> dict:
    """
        Parses a raw KEGG REST API response into a structured dictionary.

        This function reads the flat-file format used by KEGG, identifying tags
        (like ENTRY, NAME, PATHWAY) and capturing their associated data. It
        utilizes a temporary file for memory-efficient processing of large records.

        Args:
            request_text (str): The raw text response from a KEGG REST API call.

        Returns:
            dict: A processed dictionary containing the parsed and cleaned
                fields of the KEGG record.
    """
    res = io.StringIO(request_text)
    with tempfile.NamedTemporaryFile(delete=False, mode="w+", encoding="utf-8") as file_path:
        shutil.copyfileobj(res, file_path)
        file_path_name = file_path.name

    current_key = None
    saved_rec = dict()

    try:
        with open(file_path_name, "r") as entry_file:
            for line in entry_file:
                if line.startswith("///"):
                    break
                tag = line[:12].strip()
                value = line[12:].strip()
                if tag:
                    current_key = tag
                    saved_rec[current_key] = [value]
                else:
                    saved_rec[current_key].append(value)
        cleaned_recs = clean_entry(saved_rec)
    finally:
        # Putting this in 'finally' ensures the temp file is deleted
        # even if the parsing logic above crashes.
        if os.path.exists(file_path_name):
            os.remove(file_path_name)

    return cleaned_recs
