# Welcome to Kegger

A ligthweight tool for interfacing with the KEGG database.

Retrieve KEGG pathways by organism id code (e.g. eco, hsa) for use in dataframes or databases.

Uses caching for speeding up repeated calls for the same KEGG id.

## Commands
Easy example to get the pathway for E. coli's eco00010, glycolysis pathway. Below get_path() calls the KEGG API and returns the raw entry.

The next function, kegg_parser() loops through the request and converts the entry into a dictionary.

You can use the dictionary as an input into a dataframe or whatever data system you are interested (e.g Django model)

    import kegger as kg
    org = "eco"
    kd = "eco00010"
    path_record = kg.get_path(kd)
    path_dict = kg.kegg_parser(path_record)
    print(path_dict.keys())
    print(path_dict.get("GENE"))

Once you run the code, you may notice a new file being created, kegg_cache.sqlite.
This will hold all the previous requests to the KEGG API for up to 30 days.

## Project layout

    kegger/
        kegg_tools.py # the main library
