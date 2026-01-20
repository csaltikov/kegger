import shutil
import tempfile
import requests
import requests_cache
import pandas as pd
import io
import os
from collections import defaultdict
from datetime import timedelta

requests_cache.install_cache("kegg_cache", expire_after=timedelta(days=30))


def get_url(url, type="text"):
    response = requests.get(url)
    response.raise_for_status()
    if type == "json":
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


def list_all_pathways(org) -> pd.DataFrame:
    url = f'https://rest.kegg.jp/list/pathway/{org}'
    response = get_url(url)
    record = io.StringIO(response)
    col_names = ["pathid", "description"]
    df = pd.read_csv(record, sep="\t", header=None, names=col_names)
    return df


def genes_to_pathways(org) -> pd.DataFrame:
    url = f'https://rest.kegg.jp/link/{org}/pathway'
    response = get_url(url)
    record = io.StringIO(response)
    col_names = ["pathid", "gene"]
    df = pd.read_csv(record, sep="\t", header=None, names=col_names)
    return df


def get_module(mdid):
    url = f"https://rest.kegg.jp/get/md:{mdid}"
    request = get_url(url, type="text")
    return request


def get_path(pathid):
    url = f"https://rest.kegg.jp/get/{pathid}"
    request = get_url(url, type="text")
    return request


def get_entry(entry_id):
    url = f"https://rest.kegg.jp/get/{entry_id}"
    request = get_url(url, type="text")
    return request


def get_org(org: str) -> pd.DataFrame:
    url = f"https://rest.kegg.jp/list/{org}"
    response = requests.get(url)
    record = io.StringIO(response.text)
    cols = ["gene", "feature", "position", "annotation"]
    df = pd.read_csv(record, sep="\t", header=None, names=cols)
    return df


def clean_entry(entry: dict) -> dict:
    cleaned_entry = defaultdict(list)
    for tag, value in entry.items():
        if tag == "ENTRY":
            cleaned_entry[tag] = value[0].split() #re.sub(r"(\s+)Pathway", "", value[0])
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


def kegg_parser(request_text: str, force_refresh=False) -> dict:
    # Hashing ensures the key is always the same length regardless of input size
    res = io.StringIO(request_text)
    with tempfile.NamedTemporaryFile(delete=False, mode="w+", encoding="utf-8") as file_path:
        shutil.copyfileobj(res, file_path)
        file_path_name = file_path.name

    current_key = None
    saved_rec = dict()
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
    os.remove(file_path_name)
    return cleaned_recs


if __name__ == "__main__":
    ##
    org = "eco"
    kd = "00190"
    path_record = get_path(f"{org}{kd}")
    path_dict = kegg_parser(path_record)
    print(path_dict.keys())
    print(path_dict.get("GENE"))

    ##
    org_pathways = list_all_pathways(org)
    print(org_pathways.head())
    pathway_mapper = org_pathways.set_index("pathid")["description"]

    ##
    genes_paths = genes_to_pathways(org)
    # genes_paths["gene"] = genes_paths["gene"].str.replace(f"{org}:", "")
    genes_paths["pathid"] = genes_paths["pathid"].str.replace("path:", "")
    genes_paths["description"] = genes_paths["pathid"].map(pathway_mapper)
    print(genes_paths.head())

    ##
    org_annotations = get_org(org)
    print(org_annotations.head())

    ##
    res = genes_paths.merge(org_annotations, on="gene", how="left")
    print(res.shape)
