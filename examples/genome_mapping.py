from kegger import kegg_tools as kg
import pandas as pd


def main():
    org = "eco"  # E. coli
    print(f"--- Annotating Genome for {org} ---")

    # 1. Get the list of all genes and their annotations
    # Returns: gene, feature, position, annotation
    genome_df = kg.get_org(org)

    # 2. Get the mapping of genes to pathways
    # Returns: pathid, gene
    gene_to_path = kg.genes_to_pathways(org)

    # 3. Get the human-readable names for those pathways
    # Returns: pathid, description
    path_names = kg.list_all_pathways(org)

    # 4. Merge them all together
    # First, link pathway IDs to their names
    path_map = gene_to_path.merge(path_names, on="pathid")

    # Then, link those pathways to the genome annotations
    final_df = genome_df.merge(path_map, on="gene", how="left")

    print("\nTop 5 annotated genes with their pathways:")
    print(final_df[['gene', 'annotation', 'description']].head())

    # Optional: Save to a file for use in Excel/R
    # final_df.to_csv(f"{org}_full_annotation.tsv", sep="\t", index=False)

    return final_df


if __name__ == "__main__":
    res = main()
