from multiprocessing.pool import ThreadPool
import streamlit as st
import numpy as np
import pandas as pd

from collections import OrderedDict

from src.data.load_data import (
    load_disease_from_huge,
    load_disease_from_kegg,
    load_disease_pathway_from_kegg,
    load_mapping_dataset,
    load_pathway_from_kegg,
)
from src.data.ncbi_scraping import get_all_genes_info

from st_aggrid import AgGrid
from src.processing.ncbi_info import merge_gene_symbol
from src.processing.pathway import get_pathway_by_genes

from src.utils.data_table import get_base_grid_options

st.set_page_config(
    page_title="Search Disease", page_icon="🍫",
)


def _max_width_():
    max_width_str = f"max-width: 1400px;"
    st.markdown(
        f"""
    <style>
    .reportview-container .main .block-container{{
        {max_width_str}
    }}
    </style>    
    """,
        unsafe_allow_html=True,
    )


_max_width_()

c30, c31, c32 = st.columns([2.5, 1, 3])

with c30:
    # st.image("logo.png", width=400)
    st.title("🍫 Search Disease")
    st.header("")


with st.expander("ℹ️ - About this app", expanded=True):

    st.write(
        """     
-   To determine whether SNPs contain genes that appear in research on the disease of interest.
-   How?
    -   a
    -   b
	    """
    )

    st.markdown("")

st.markdown("")
st.markdown("## **📌 Paste data **")
with st.form(key="my_form"):

    ce, c1, ce, c2, c3 = st.columns([0.07, 2, 0.07, 5, 0.07])
    with c1:
        st.header("Parameters")

        st.subheader("Region")
        upstream_threshold = st.number_input(
            "Upstream",
            value=2_000,
            min_value=0,
            max_value=5_000_000,
            step=1000,
            help="Upstream threshold (default: 0)",
        )
        downstream_threshold = st.number_input(
            "Downstream",
            value=500,
            min_value=0,
            max_value=5_000_000,
            step=1000,
            help="Downstream threshold (default: 0)",
        )

        st.subheader("Disease")
        if not "diseases_list" in st.session_state:
            st.session_state.diseases_list = OrderedDict(
                [
                    ("Bipolar disorder", True),
                    ("Coronary artery disease", True),
                    ("Crohn disease", True),
                    ("High blood pressure", True),
                    ("Rheumatoid arthritis", True),
                    ("Type 1 diabetes mellitus", True),
                    ("Type 2 diabetes mellitus", True),
                ]
            )
        for key, value in st.session_state.diseases_list.items():
            st.checkbox(key, value)

    with c2:
        uploaded_file = st.file_uploader(
            "",
            type=["txt", "csv", "tsv", "xls", "xlsx", "json"],
            key="uploaded_file",
            help="Upload a text file or a CSV/TSV/XLS/XLSX/JSON file.",
        )

        st.markdown("")
        submit_button = st.form_submit_button(label="✨ Get me the data!")

        if uploaded_file is not None:
            file_container = st.expander("Check your uploaded .csv")
            shows = pd.read_csv(uploaded_file)
            uploaded_file.seek(0)
            file_container.write(shows)
            st.info(f"Found : {len(shows)} records")


if not (submit_button and uploaded_file):
    st.stop()

st.header("")
st.markdown("## **🎈 Results **")

data = load_mapping_dataset()
disease_huge_data = load_disease_from_huge()
disease_kegg_data = load_disease_from_kegg()

pathway_kegg_data = load_pathway_from_kegg()
disease_pathway_kegg_data = load_disease_pathway_from_kegg()


if all(item in shows.columns for item in ["Probe Set ID", "dbSNP RS ID"]):
    selected_data = data[
        data["Probe Set ID"].isin(shows["Probe Set ID"])
        | data["dbSNP RS ID"].isin(shows["dbSNP RS ID"])
    ]
elif {"Probe Set ID"}.issubset(shows.columns):
    selected_data = data[data["Probe Set ID"].isin(shows["Probe Set ID"])]
elif {"dbSNP RS ID"}.issubset(shows.columns):
    selected_data = data[data["dbSNP RS ID"].isin(shows["dbSNP RS ID"])]
else:
    st.error("No Probe Set ID or dbSNP RS ID found in the data")
    st.stop()

# filter bound of region (upstream and downstream)
is_inbound_upstream = (selected_data["region"] == "upstream") & (
    selected_data["distance"] <= upstream_threshold
)
is_inbound_downstream = (selected_data["region"] == "downstream") & (
    selected_data["distance"] <= downstream_threshold
)
selected_data = selected_data[
    is_inbound_upstream
    | is_inbound_downstream
    | (~selected_data["region"].isin(["upstream", "downstream"]))
]
selected_data["gene_id"] = selected_data["gene_id"].astype(str)


# selected_data = selected_data.append(
#     [{"gene_id": g} for g in (["3123", "3119", "3117", "3630"] + ["1585"])],
#     ignore_index=True,
# )

st.subheader("Mapping from selected to dataset 🔎")
AgGrid(
    selected_data,
    gridOptions=get_base_grid_options(selected_data, "Probe Set ID"),
    enable_enterprise_modules=True,
)


def update_discontinued_genes_info():
    # search discontinue gene
    replaced_genes_df = st.session_state.gene_nbci_info[
        "Replaced with Gene ID"
    ].dropna()
    search_replaced_genes_df = replaced_genes_df[
        ~replaced_genes_df.isin(st.session_state.gene_nbci_info["Gene Id"])
    ]
    if not replaced_genes_df.empty and not search_replaced_genes_df.empty:
        updated_genes_info = get_all_genes_info(search_replaced_genes_df.to_list())

        st.session_state.gene_nbci_info = st.session_state.gene_nbci_info.append(
            updated_genes_info
        )


# caching gene info from NCBI searching
if "gene_interest" not in st.session_state:
    st.session_state.gene_interest = list(selected_data["gene_id"].unique())
    st.session_state.gene_nbci_info = get_all_genes_info(st.session_state.gene_interest)
    update_discontinued_genes_info()
else:
    # check new selected_data is in gene_interest
    new_gene_interest = list(
        np.setdiff1d(selected_data["gene_id"].unique(), st.session_state.gene_interest)
    )

    if len(new_gene_interest) > 0:
        new_gene_ncbi_info = get_all_genes_info(new_gene_interest)

        st.session_state.gene_interest = (
            st.session_state.gene_interest + new_gene_interest
        )
        st.session_state.gene_nbci_info = st.session_state.gene_nbci_info.append(
            new_gene_ncbi_info
        )
        update_discontinued_genes_info()


st.subheader("List of genes selected that search on NCBI 🔍")
gene_ncbi_info_selected = st.session_state.gene_nbci_info[
    st.session_state.gene_nbci_info["Gene Id"].isin(
        selected_data["gene_id"].astype(str)
    )
]
AgGrid(
    gene_ncbi_info_selected,
    gridOptions=get_base_grid_options(gene_ncbi_info_selected),
    enable_enterprise_modules=True,
)


# genes list for search on HuGE
genes_list = gene_ncbi_info_selected.apply(merge_gene_symbol, axis=1).explode().unique()

st.subheader("Gene that found on Disease (HuGE) 😎")
disease_huge_include_gene = disease_huge_data[
    disease_huge_data["Gene"].isin(genes_list)
]
AgGrid(
    disease_huge_include_gene,
    gridOptions=get_base_grid_options(disease_huge_include_gene, "Disease"),
    enable_enterprise_modules=True,
)


st.subheader("Gene that found on Disease (KEGG) 🤧")
disease_kegg_include_gene = disease_kegg_data[
    disease_kegg_data["gene_symbol"].isin(genes_list)
]
AgGrid(
    disease_kegg_include_gene,
    gridOptions=get_base_grid_options(disease_kegg_include_gene, "disease"),
    enable_enterprise_modules=True,
)


st.subheader("Gene that found on Pathway (KEGG) 🪶")
pathway_genes_with_info = get_pathway_by_genes(
    pathway_kegg_data, disease_pathway_kegg_data, gene_ncbi_info_selected
)
AgGrid(
    pathway_genes_with_info,
    gridOptions=get_base_grid_options(pathway_genes_with_info, "disease"),
    enable_enterprise_modules=True,
)

with st.expander("View raw data"):
    st.subheader("HuGE navigator (Disease <-> Gene)")
    AgGrid(
        disease_huge_data,
        gridOptions=get_base_grid_options(disease_huge_data, "Disease"),
        enable_enterprise_modules=True,
    )

    st.subheader("KEGG (Disease <-> Gene)")
    AgGrid(
        disease_kegg_data,
        gridOptions=get_base_grid_options(disease_kegg_data, "disease"),
        enable_enterprise_modules=True,
    )

    st.subheader("KEGG (Disease <-> Pathway)")
    AgGrid(
        disease_pathway_kegg_data,
        gridOptions=get_base_grid_options(disease_pathway_kegg_data, "disease"),
        enable_enterprise_modules=True,
        key="disease_pathway_kegg_data",
    )

    st.subheader("KEGG (Pathway <-> Gene)")
    AgGrid(
        pathway_kegg_data,
        gridOptions=get_base_grid_options(pathway_kegg_data),
        enable_enterprise_modules=True,
        key="pathway_kegg_data",
    )
