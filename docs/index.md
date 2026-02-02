# Welcome to Kegger
[![PyPI version](https://img.shields.io/pypi/v/kegger.svg)](https://pypi.org/project/kegger/)

A lightweight tool for interfacing with the KEGG (Kyoto Encyclopedia of Genes and Genomes) database.

`kegger` allows you to retrieve biological pathways, genome annotations, and gene-to-pathway mappings for use in Pandas DataFrames or downstream databases (like Django or SQLAlchemy).

## Key Features
* **Automated Parsing:** Converts tricky KEGG flat-files into clean Python dictionaries.
* **Pandas Integration:** Native support for converting KEGG lists into DataFrames.
* **Smart Caching:** Optional SQLite caching to speed up repeated calls and reduce server load.

## Quick Start
To get started, you can initialize the optional cache and fetch a specific pathway record:

```python
import kegger as kg

# Optional: Enable caching for 30 days
kg.initialize_kegger(cache_path="my_kegg_cache", expire_days=30)

# Fetch Glycolysis for E. coli
org = "eco"
path_id = "eco00010"  # or 'path:eco00010'
path_record = kg.get_path(path_id)

# Parse the raw text into a structured dictionary
path_dict = kg.kegg_parser(path_record)

print(path_dict.keys())
print(path_dict.get("NAME"))
print(path_dict.get("GENE")[:5]) # Show first 5 genes
```
```text
dict_keys(['ENTRY', 'NAME', 'DESCRIPTION', 'CLASS', 'PATHWAY_MAP', 'MODULE',
            'DBLINKS', 'ORGANISM', 'GENE', 'ORTHOLOG', 'COMPOUND', 'REFERENCE',
            'AUTHORS', 'TITLE', 'JOURNAL', 'REL_PATHWAY', 'KO_PATHWAY'])

Glycolysis / Gluconeogenesis - Escherichia coli K-12 MG1655

['b0114', 'b0115', 'b0116', 'b0325', 'b0356']
```

## Project layout

    kegger/
        kegg_tools.py    # Core library functions
    examples/
        basic_usage.py   # Starter script
        genome_mapping.py # Advanced DataFrame examples

## Data Attribution and Citations

This package retrieves data from the **KEGG (Kyoto Encyclopedia of Genes and Genomes)** REST API. 

* **Data Source:** [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html)
* **Terms of Use:** Please note that KEGG is for academic use by academic users at academic institutions. For-profit users may require a license. See the [KEGG Terms of Use](https://www.kegg.jp/kegg/legal.html) for details.

If you use data retrieved via `kegger` in a publication, please cite the primary KEGG reference:

> Kanehisa, M. and Goto, S.; **KEGG: kyoto encyclopedia of genes and genomes.** *Nucleic Acids Res.* 28, 27-30 (2000). [doi:10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)

## Site Downloads
[![Downloads](https://static.pepy.tech/badge/kegger/month)](https://pepy.tech/project/kegger)
[![Total Downloads](https://static.pepy.tech/badge/kegger)](https://pepy.tech/project/kegger)
