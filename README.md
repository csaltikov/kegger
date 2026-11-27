# kegger

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://csaltikov.github.io/kegger/)
[![PyPI](https://img.shields.io/pypi/v/kegger)](https://pypi.org/project/kegger/)
[![Python Version](https://img.shields.io/pypi/pyversions/kegger)](https://pypi.org/project/kegger/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Interact with the KEGG (Kyoto Encyclopedia of Genes and Genomes) API using Python and Pandas.

Retrieve genes, pathways, and associated data for organisms of interest in the KEGG database. Originally developed for RNA-seq analysis to organize gene expression patterns by KEGG pathways.

## Features

- 🧬 Retrieve pathway information for any organism in KEGG
- 📊 Returns data as Pandas DataFrames for easy analysis
- 🔍 Parse KEGG records into structured dictionaries
- 💾 Built-in request caching for faster repeated queries
- 🐍 Simple, intuitive Python API

## Installation

Install from PyPI using pip:

```bash
pip install kegger
```

Or using uv:

```bash
uv add kegger
```

## Quick Start

### Get all pathways for an organism

```python
import kegger as kg

org = "eco"  # E. coli K-12 MG1655
org_pathways = kg.list_all_pathways(org)

print(org_pathways.head())
```

**Output:**
```
     pathid                                        description
0  eco01100  Metabolic pathways - Escherichia coli K-12 MG1655
1  eco01110  Biosynthesis of secondary metabolites - Escher...
2  eco01120  Microbial metabolism in diverse environments -...
3  eco01200   Carbon metabolism - Escherichia coli K-12 MG1655
4  eco01210  2-Oxocarboxylic acid metabolism - Escherichia ...
```
### Get all genes and their pathway associations

```python
genes_pathids = kg.genes_to_pathways(org)
genes_pathids = genes_pathids.merge(org_pathways, how="outer")
print(genes_pathids.head())
```
**Output**
```
     pathid   gene                                        description
0  eco00010  b0114  Glycolysis / Gluconeogenesis - Escherichia col...
1  eco00010  b0115  Glycolysis / Gluconeogenesis - Escherichia col...
2  eco00010  b0116  Glycolysis / Gluconeogenesis - Escherichia col...
3  eco00010  b0325  Glycolysis / Gluconeogenesis - Escherichia col...
4  eco00010  b0356  Glycolysis / Gluconeogenesis - Escherichia col...
```

### Get details for a specific pathway from a KEGG Entry
If you need more details about a pathway, use the get_path() function to retrieve the raw text from the KEGG request and kegg_parser() to convert the string into a dictionary.

```python
pathid = "eco00190"  # Oxidative phosphorylation
path_record = kg.get_path(pathid) # the kegg request as a string
path_dict = kg.kegg_parser(path_record)

# Available data fields
print(path_dict.keys())
# dict_keys(['ENTRY', 'NAME', 'CLASS', 'PATHWAY_MAP', 'MODULE', 
#            'ORGANISM', 'GENE', 'ORTHOLOG', 'COMPOUND', 
#            'REFERENCE', 'AUTHORS', 'TITLE', 'JOURNAL', 'KO_PATHWAY'])

# Show the pathway name
print(path_dict.get("NAME"))
# Oxidative phosphorylation - Escherichia coli K-12 MG1655

# Show associated modules
print(path_dict.get("MODULE"))
# [
# 'eco_M00144  NADH:quinone oxidoreductase, prokaryotes [PATH:eco00190]',
# 'eco_M00149  Succinate dehydrogenase, prokaryotes [PATH:eco00190]',
# 'eco_M00150  Fumarate reductase, prokaryotes [PATH:eco00190]',
# 'eco_M00153  Cytochrome bd ubiquinol oxidase [PATH:eco00190]',
# 'eco_M00157  F-type ATPase, prokaryotes and chloroplasts [PATH:eco00190]', 
# 'eco_M00417  Cytochrome o ubiquinol oxidase [PATH:eco00190]'
# ]

# Get genes in this pathway
genes = path_dict.get("GENE")
print(f"Found {len(genes)} genes in pathway")
print(genes[:5])  # First 5 genes
# ['b0428', 'b0429', 'b0430', 'b0431', 'b0432']
```

### Get annotation and KEGG pathways for all genes

```python
import kegger as kg

org = "eco"
annotated_org = kg.get_org(org)
genes_pathways_org = annotated_org.merge(org_pathways, how="outer")
print()

```
## Documentation

Full documentation is available at [https://csaltikov.github.io/kegger/](https://csaltikov.github.io/kegger/)

## Requirements

- Python ≥ 3.13
- pandas ≥ 2.2.3
- requests ≥ 2.32.5
- requests-cache ≥ 1.2.1

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use kegger in your research, please cite:

```
Saltikov, C. (2025). kegger: A Python package for interacting with the KEGG API. 
GitHub repository: https://github.com/csaltikov/kegger
```

## Acknowledgments

This package interfaces with the [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html). Please cite KEGG appropriately in your work:

> Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000).
