# Welcome to Kegger

A ligthweight tool for interfacing with the KEGG database

## Commands
Easy example to get the pathway for E. coli's eco00010, glycolysis

    import kegger as kg
    org = "eco"
    kd = "00010"
    path_record = kg.get_path(f"{org}{kd}")
    path_dict = kg.kegg_parser(path_record)
    print(path_dict.keys())
    print(path_dict.get("GENE"))


## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.
