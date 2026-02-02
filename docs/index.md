# Welcome to Kegger

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
kg.initialize_kegger(cache_name="my_kegg_cache", expire_days=30)

# Fetch Glycolysis for E. coli
org = "eco"
path_id = "eco00010"  # or 'path:eco00010'
path_record = kg.get_path(path_id)

# Parse the raw text into a structured dictionary
path_dict = kg.kegg_parser(path_record)

print(path_dict.get("NAME"))
print(path_dict.get("GENE")[:5]) # Show first 5 genes
```

## Project layout

    kegger/
        kegg_tools.py    # Core library functions
    examples/
        basic_usage.py   # Starter script
        genome_mapping.py # Advanced DataFrame examples
