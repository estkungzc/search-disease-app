from pandas import DataFrame
import pandas as pd


def get_pathway_by_genes(
    pathway_kegg_data: DataFrame,  # lookup path:hsa <-> hsa
    disease_pathway_kegg_data: DataFrame,  # look up path:hsa <-> disease
    gene_ncbi_info_selected: DataFrame,  # to searching path:hsa by gene_id
):
    pathway_by_genes_result = pathway_kegg_data[
        pathway_kegg_data["path:hsa"].isin(disease_pathway_kegg_data["path:hsa"])
        & pathway_kegg_data["hsa"].isin(gene_ncbi_info_selected["Gene Id"])
    ]
    pathway_with_info = pd.merge(
        pathway_by_genes_result, disease_pathway_kegg_data, on=["path:hsa"]
    )

    pathway_genes_with_info = pd.merge(
        pathway_with_info,
        gene_ncbi_info_selected,
        how="left",
        left_on=["hsa"],
        right_on=["Gene Id"],
    )

    pathway_genes_with_info = pathway_genes_with_info.rename(
        columns={
            "name": "pathway_name",
            "Gene Id": "gene_id",
            "Official Symbol": "gene_symbol",
        }
    )
    pathway_genes_with_info = pathway_genes_with_info[
        ["disease", "gene_id", "gene_symbol", "path:hsa", "pathway_name", "path_type",]
    ]

    return pathway_genes_with_info
