# kegger
Light weight python tool for interacting with KEGG API.

Retrieve genes and KEGG pathway for you favorite organisms in the KEGG database

```
org = "eco"
kd = "00190"  # Oxidative phosphorylation - Escherichia coli K-12 MG1655
path_record = get_path(f"{org}{kd}")
path_dict = kegg_parser(path_record)
print(path_dict.get("GENE"))

['b0428', 'b0429', 'b0430', 'b0431', 'b0432', ...]
```
