# kegger
Interact with the KEGG API using Python and Pandas

Retrieve genes and KEGG pathway for your favorite organisms in the KEGG database

I developed this primarily for use in RNAseq analysis and organizing gene expression patterns by KEGG pathways

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://csaltikov.github.io/kegger/)

### Here are some quick examples:
Get all pathways from E. coli
```python
>>> import kegger as kg
>>> org = "eco"

>>> org_pathways = kg.list_all_pathways(org)
>>> org_pathways.head()
     pathid                                        description
0  eco01100  Metabolic pathways - Escherichia coli K-12 MG1655
1  eco01110  Biosynthesis of secondary metabolites - Escher...
2  eco01120  Microbial metabolism in diverse environments -...
3  eco01200   Carbon metabolism - Escherichia coli K-12 MG1655
4  eco01210  2-Oxocarboxylic acid metabolism - Escherichia ...
```
Get a specific pathway:

```python
>>> import kegger as kg
>>> pathid = "eco00190" # Oxidative phosphorylation - Escherichia coli K-12 MG1655
>>> org = "eco"
>>> path_record = kg.get_path(pathid)
>>> path_dict = kg.kegg_parser(path_record)

>>> print(path_dict.keys())
dict_keys(['ENTRY', 'NAME', 'CLASS', 'PATHWAY_MAP', 'MODULE', 'ORGANISM', 'GENE', 'ORTHOLOG', 'COMPOUND', 
           'REFERENCE', 'AUTHORS', 'TITLE', 'JOURNAL', 'KO_PATHWAY'])

>>> path_dict.get("GENE")
['b0428', 'b0429', 'b0430', 'b0431', 'b0432', 'b0721', 'b0722', 'b0723', 'b0724', 'b0733', 'b0734', 'b0978', 
 'b0979', 'b1109', 'b2276', 'b2277', 'b2278', 'b2279', 'b2280', 'b2281', 'b2282', 'b2283', 'b2284', 'b2285', 
 'b2286', 'b2287', 'b2288', 'b2501', 'b3731', 'b3732', 'b3733', 'b3734', 'b3735', 'b3736', 'b3737', 'b3738', 
 'b4151', 'b4152', 'b4153', 'b4154', 'b4226', 'b4515', 'b4592']
```
