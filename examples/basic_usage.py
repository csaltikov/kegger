from kegger import kegg_tools as kg

# 1. Initialize the cache so we don't spam KEGG
kg.initialize_kegger(cache_path="kegg_cache", expire_days=7)


def main():
    org = "eco"
    print(f"--- Fetching pathway data for {org} ---")

    # Get pathways
    org_pathways = kg.list_all_pathways(org)
    print(org_pathways.head())

    # Map genes to pathways
    genes_paths = kg.genes_to_pathways(org)
    print("\n--- Gene to Pathway Mapping ---")
    print(genes_paths.head())


if __name__ == "__main__":
    main()
