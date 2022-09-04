import streamlit as st
import numpy as np
import pandas as pd
import re


from src.data.load_data import (
    load_disease_from_huge,
    load_disease_from_kegg,
    load_disease_pathway_from_kegg,
    load_mapping_dataset,
    load_pathway_from_kegg,
)
from src.data.ncbi_scraping import get_all_genes_info

from st_aggrid import AgGrid

from src.utils.data_table import get_base_grid_options

import streamlit_authenticator as stauth

st.set_page_config(
    page_title="Search Disease", page_icon="🍫",
)


names = [st.secrets["app_credentials"]["name"]]
usernames = [st.secrets["app_credentials"]["username"]]
passwords = [st.secrets["app_credentials"]["password"]]


hashed_passwords = stauth.Hasher(passwords).generate()
authenticator = stauth.Authenticate(
    names,
    usernames,
    hashed_passwords,
    "some_cookie_name",
    "some_signature_key",
    cookie_expiry_days=30,
)

with st.sidebar:
    name, authentication_status, username = authenticator.login("Login", "main")
    if authentication_status:
        authenticator.logout("Logout", "main")
        st.write("Welcome *%s*" % (name))
    elif authentication_status == False:
        st.error("Username/password is incorrect")
    elif authentication_status == None:
        st.warning("Please enter your username and password")


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
    -   Enter SNP id, e.g. rs1234 then enter threshold value for upstream and downstream region.
    -   Then, the SNP will categorize the SNPs into the following categories: SNP near or within a gene, SNP outside genes.
    -   Finally, SNP that category is "SNP near or within a gene"  will use to search disease on HUGE, KEGG or KEGG pathways.
	    """
    )

    st.markdown("")

if not authentication_status:
    st.stop()

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
        st.write(
            """     
            ✅ Bipolar disorder \n
            ✅ Coronary artery disease \n
            ✅ Crohn disease \n
            ✅ High blood pressure \n
            ✅ Rheumatoid arthritis \n
            ✅ Type 1 diabetes mellitus \n
            ✅ Type 2 diabetes mellitus
            """
        )

    with c2:
        uploaded_file = st.file_uploader(
            "",
            type=["txt", "csv", "tsv", "xls", "xlsx", "json"],
            key="uploaded_file",
            help="Upload a text file or a CSV/TSV/XLS/XLSX/JSON file.",
        )

        st.markdown("")
        snp_input = st.text_area("Input rs ID", height=500)
        st.markdown("")
        submit_button = st.form_submit_button(label="✨ Get me the data!")

        shows = None
        if uploaded_file is not None:
            file_container = st.expander("Check your uploaded .csv")
            shows = pd.read_csv(uploaded_file)
            uploaded_file.seek(0)
            file_container.write(shows)
            st.info(f"Found : {len(shows)} records")


if not (submit_button and (uploaded_file or snp_input)):
    st.stop()


# ░█▀▀█ █▀▀ █▀▀ █──█ █── ▀▀█▀▀ █▀▀
# ░█▄▄▀ █▀▀ ▀▀█ █──█ █── ──█── ▀▀█
# ░█─░█ ▀▀▀ ▀▀▀ ─▀▀▀ ▀▀▀ ──▀── ▀▀▀

st.header("")
st.markdown("## **🎈 Results **")
st.markdown("")

# load all local data
data = load_mapping_dataset()
disease_huge_data = load_disease_from_huge()
disease_kegg_data = load_disease_from_kegg()

pathway_kegg_data = load_pathway_from_kegg()
disease_pathway_kegg_data = load_disease_pathway_from_kegg()


filter_data = None

if shows is not None:
    if all(item in shows.columns for item in ["Probe Set ID", "dbSNP RS ID"]):
        # selected_data = data[
        #     data["Probe Set ID"].isin(shows["Probe Set ID"])
        #     | data["dbSNP RS ID"].isin(shows["dbSNP RS ID"])
        # ]
        filter_data = data["Probe Set ID"].isin(shows["Probe Set ID"]) | data[
            "dbSNP RS ID"
        ].isin(shows["dbSNP RS ID"])
    elif {"Probe Set ID"}.issubset(shows.columns):
        # selected_data = data[data["Probe Set ID"].isin(shows["Probe Set ID"])]
        filter_data = data["Probe Set ID"].isin(shows["Probe Set ID"])
    elif {"dbSNP RS ID"}.issubset(shows.columns):
        # selected_data = data[data["dbSNP RS ID"].isin(shows["dbSNP RS ID"])]
        filter_data = data["dbSNP RS ID"].isin(shows["dbSNP RS ID"])
    else:
        st.error("No Probe Set ID or dbSNP RS ID found in the data")
        st.stop()

if snp_input:
    snp_list = re.split("; |, |\n| ", snp_input)

    if snp_list:
        snp_list = set(snp_list)
        selected_from_snp_data = data["dbSNP RS ID"].isin(snp_list) | data[
            "Probe Set ID"
        ].isin(snp_list)
        if filter_data is not None:
            filter_data = filter_data | selected_from_snp_data
        else:
            filter_data = selected_from_snp_data

selected_data = data[filter_data]

# print(selected_data)


selected_data["gene_id"] = selected_data["gene_id"].astype(str)


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

# update discontinued genes info in raw data
v = st.session_state.gene_nbci_info[
    ~st.session_state.gene_nbci_info["Replaced with Gene ID"].isna()
]
for index, row in v.iterrows():
    new_data = st.session_state.gene_nbci_info[
        st.session_state.gene_nbci_info["Gene Id"] == row["Replaced with Gene ID"]
    ].iloc[0]
    selected_data["gene_symbol"] = np.where(
        selected_data["gene_id"] == row["Gene Id"],
        new_data["Official Symbol"],
        selected_data["gene_symbol"],
    )
    selected_data["gene_id"] = np.where(
        selected_data["gene_id"] == row["Gene Id"],
        new_data["Gene Id"],
        selected_data["gene_id"],
    )


# ░█▀▄▀█ █▀▀█ █▀▀█ █▀▀█ ─▀─ █▀▀▄ █▀▀▀ 　 █▀▀ █▀▀█ █▀▀█ █▀▄▀█ 　 █▀▀ █▀▀ █── █▀▀ █▀▀ ▀▀█▀▀ █▀▀ █▀▀▄ 　 ▀▀█▀▀ █▀▀█
# ░█░█░█ █▄▄█ █──█ █──█ ▀█▀ █──█ █─▀█ 　 █▀▀ █▄▄▀ █──█ █─▀─█ 　 ▀▀█ █▀▀ █── █▀▀ █── ──█── █▀▀ █──█ 　 ──█── █──█
# ░█──░█ ▀──▀ █▀▀▀ █▀▀▀ ▀▀▀ ▀──▀ ▀▀▀▀ 　 ▀── ▀─▀▀ ▀▀▀▀ ▀───▀ 　 ▀▀▀ ▀▀▀ ▀▀▀ ▀▀▀ ▀▀▀ ──▀── ▀▀▀ ▀▀▀─ 　 ──▀── ▀▀▀▀

# █▀▀▄ █▀▀█ ▀▀█▀▀ █▀▀█ █▀▀ █▀▀ ▀▀█▀▀
# █──█ █▄▄█ ──█── █▄▄█ ▀▀█ █▀▀ ──█──
# ▀▀▀─ ▀──▀ ──▀── ▀──▀ ▀▀▀ ▀▀▀ ──▀──

gene_ncbi_info_repr = pd.merge(
    st.session_state.gene_nbci_info,
    selected_data,
    left_on=["Gene Id"],
    right_on=["gene_id"],
)


def is_snp_within_gene(snp):

    is_within_upstream_threshold = (
        snp["region"] == "upstream" and snp["distance"] <= upstream_threshold
    )
    is_within_downstream_threshold = (
        snp["region"] == "downstream" and snp["distance"] <= downstream_threshold
    )

    return (
        "Y"
        if (
            (is_within_upstream_threshold or is_within_downstream_threshold)
            or (snp["region"] not in ["upstream", "downstream"])
        )
        else "N"
    )


gene_ncbi_info_repr["Is SNP's near or within a gene"] = gene_ncbi_info_repr.apply(
    is_snp_within_gene, axis=1
)

gene_ncbi_info_repr = gene_ncbi_info_repr.drop("gene_id", axis=1)
gene_ncbi_info_repr = gene_ncbi_info_repr[
    [
        "Probe Set ID",
        "dbSNP RS ID",
        "Gene Id",
        "Official Symbol",
        "region",
        "distance",
        "Official Full Name",
        "Also known as",
        "Is SNP's near or within a gene",
    ]
]

snp_within_genes = gene_ncbi_info_repr[
    gene_ncbi_info_repr["Is SNP's near or within a gene"] == "Y"
]

# Special case: SNPs can be found inside and outside the genes
snp_groups = gene_ncbi_info_repr.groupby(["Probe Set ID"])

snp_in_out_genes = []  # list of probe set id
for group_name, df_group in snp_groups:
    gene_with_info_in_out = {}
    for row_index, row in df_group.iterrows():
        # print(row_index)
        if row["Official Symbol"] not in gene_with_info_in_out:
            gene_with_info_in_out[row["Official Symbol"]] = {
                row["Is SNP's near or within a gene"]
            }
        else:
            gene_with_info_in_out[row["Official Symbol"]].add(
                row["Is SNP's near or within a gene"]
            )
    for gene, info in gene_with_info_in_out.items():
        if len(info) > 1:
            snp_in_out_genes.append(group_name)

gene_ncbi_info_snp_in_out = gene_ncbi_info_repr[
    gene_ncbi_info_repr["Probe Set ID"].isin(snp_in_out_genes)
]

gene_ncbi_info_repr = gene_ncbi_info_repr[
    ~gene_ncbi_info_repr["Probe Set ID"].isin(snp_in_out_genes)
]

# columns for represent report HuGE, KEGG
COLUMNS_LIST_VIEW = [
    "Probe Set ID",
    "dbSNP RS ID",
    "Disease",
    "Gene Id",
    "Official Symbol",
    "region",
    "distance",
    "Also known as",
    "Gene (report on disease)",
    # "Is SNP's near or within a gene",
]

# ░█▀▀█ █▀▀ █▀▀▄ █▀▀ 　 ▀▀█▀▀ █──█ █▀▀█ ▀▀█▀▀ 　 █▀▀ █▀▀ █▀▀█ █▀▀█ █▀▀ █──█ 　 █▀▀█ █▀▀▄
# ░█─▄▄ █▀▀ █──█ █▀▀ 　 ──█── █▀▀█ █▄▄█ ──█── 　 ▀▀█ █▀▀ █▄▄█ █▄▄▀ █── █▀▀█ 　 █──█ █──█
# ░█▄▄█ ▀▀▀ ▀──▀ ▀▀▀ 　 ──▀── ▀──▀ ▀──▀ ──▀── 　 ▀▀▀ ▀▀▀ ▀──▀ ▀─▀▀ ▀▀▀ ▀──▀ 　 ▀▀▀▀ ▀──▀

# ░█▀▀▄ ─▀─ █▀▀ █▀▀ █▀▀█ █▀▀ █▀▀ 　 ▄▀ ░█─░█ █──█ ░█▀▀█ ░█▀▀▀ ▀▄
# ░█─░█ ▀█▀ ▀▀█ █▀▀ █▄▄█ ▀▀█ █▀▀ 　 █─ ░█▀▀█ █──█ ░█─▄▄ ░█▀▀▀ ─█
# ░█▄▄▀ ▀▀▀ ▀▀▀ ▀▀▀ ▀──▀ ▀▀▀ ▀▀▀ 　 ▀▄ ░█─░█ ─▀▀▀ ░█▄▄█ ░█▄▄▄ ▄▀

genes_info_lookup = pd.DataFrame(
    {c: pd.Series(dtype=t) for c, t in [("gene_id", str), ("gene_symbol", str)]}
)
for index, row in snp_within_genes.iterrows():
    result = {"gene_id": row["Gene Id"], "gene_symbol": row["Official Symbol"]}
    genes_info_lookup = genes_info_lookup.append(result, ignore_index=True)

    if row["Also known as"]:
        for gene_name in row["Also known as"]:
            result = {"gene_id": row["Gene Id"], "gene_symbol": gene_name.upper()}
            genes_info_lookup = genes_info_lookup.append(result, ignore_index=True)

genes_info_lookup = genes_info_lookup.drop_duplicates(keep=False)

genes_disease_info = pd.merge(
    genes_info_lookup,
    disease_huge_data,
    # how="left",
    left_on=["gene_symbol"],
    right_on=["Gene"],
)

genes_disease_info = genes_disease_info[genes_disease_info["Gene"].notna()]
genes_disease_info = genes_disease_info[["gene_id", "Gene", "Disease"]]

disease_huge_repr = pd.merge(
    snp_within_genes,
    genes_disease_info,
    how="left",
    left_on=["Gene Id"],
    right_on=["gene_id"],
)

disease_huge_repr.rename(columns={"Gene": "Gene (report on disease)"}, inplace=True)
disease_huge_repr = disease_huge_repr[COLUMNS_LIST_VIEW]

# ░█▀▀█ █▀▀ █▀▀▄ █▀▀ 　 ▀▀█▀▀ █──█ █▀▀█ ▀▀█▀▀ 　 █▀▀ █▀▀ █▀▀█ █▀▀█ █▀▀ █──█ 　 █▀▀█ █▀▀▄
# ░█─▄▄ █▀▀ █──█ █▀▀ 　 ──█── █▀▀█ █▄▄█ ──█── 　 ▀▀█ █▀▀ █▄▄█ █▄▄▀ █── █▀▀█ 　 █──█ █──█
# ░█▄▄█ ▀▀▀ ▀──▀ ▀▀▀ 　 ──▀── ▀──▀ ▀──▀ ──▀── 　 ▀▀▀ ▀▀▀ ▀──▀ ▀─▀▀ ▀▀▀ ▀──▀ 　 ▀▀▀▀ ▀──▀

# ░█▀▀▄ ─▀─ █▀▀ █▀▀ █▀▀█ █▀▀ █▀▀ 　 ▄▀ ░█─▄▀ ░█▀▀▀ ░█▀▀█ ░█▀▀█ ▀▄
# ░█─░█ ▀█▀ ▀▀█ █▀▀ █▄▄█ ▀▀█ █▀▀ 　 █─ ░█▀▄─ ░█▀▀▀ ░█─▄▄ ░█─▄▄ ─█
# ░█▄▄▀ ▀▀▀ ▀▀▀ ▀▀▀ ▀──▀ ▀▀▀ ▀▀▀ 　 ▀▄ ░█─░█ ░█▄▄▄ ░█▄▄█ ░█▄▄█ ▄▀

genes_disease_kegg_info = pd.merge(
    snp_within_genes,
    disease_kegg_data,
    how="left",
    left_on=["Gene Id"],
    right_on=["gene_id"],
)
genes_disease_kegg_info.drop(["gene_id", "gene_name"], axis=1, inplace=True)
genes_disease_kegg_info.rename(
    columns={"disease": "Disease", "gene_symbol": "Gene (report on disease)"},
    inplace=True,
)

genes_disease_kegg_info = genes_disease_kegg_info[COLUMNS_LIST_VIEW]


# ░█▀▀█ █▀▀ █▀▀▄ █▀▀ 　 ▀▀█▀▀ █──█ █▀▀█ ▀▀█▀▀ 　 █▀▀ █▀▀ █▀▀█ █▀▀█ █▀▀ █──█ 　 █▀▀█ █▀▀▄
# ░█─▄▄ █▀▀ █──█ █▀▀ 　 ──█── █▀▀█ █▄▄█ ──█── 　 ▀▀█ █▀▀ █▄▄█ █▄▄▀ █── █▀▀█ 　 █──█ █──█
# ░█▄▄█ ▀▀▀ ▀──▀ ▀▀▀ 　 ──▀── ▀──▀ ▀──▀ ──▀── 　 ▀▀▀ ▀▀▀ ▀──▀ ▀─▀▀ ▀▀▀ ▀──▀ 　 ▀▀▀▀ ▀──▀

# ░█▀▀█ █▀▀█ ▀▀█▀▀ █──█ █───█ █▀▀█ █──█ 　 ▄▀ ░█─▄▀ ░█▀▀▀ ░█▀▀█ ░█▀▀█ ▀▄
# ░█▄▄█ █▄▄█ ──█── █▀▀█ █▄█▄█ █▄▄█ █▄▄█ 　 █─ ░█▀▄─ ░█▀▀▀ ░█─▄▄ ░█─▄▄ ─█
# ░█─── ▀──▀ ──▀── ▀──▀ ─▀─▀─ ▀──▀ ▄▄▄█ 　 ▀▄ ░█─░█ ░█▄▄▄ ░█▄▄█ ░█▄▄█ ▄▀

pathway_kegg_data_used = pathway_kegg_data[
    pathway_kegg_data["path:hsa"].isin(disease_pathway_kegg_data["path:hsa"])
]
pathway_with_info = pd.merge(
    pathway_kegg_data_used, disease_pathway_kegg_data, on=["path:hsa"]
)

pathway_kegg_genes_with_info = pd.merge(
    snp_within_genes,
    pathway_with_info,
    how="left",
    left_on=["Gene Id"],
    right_on=["hsa"],
)
pathway_kegg_genes_with_info.drop("hsa", axis=1, inplace=True)
pathway_kegg_genes_with_info.rename(
    columns={"name": "Pathway Name", "disease": "Disease", "path_type": "Path Type"},
    inplace=True,
)

# ░█▀▀▀█ █──█ █▀▄▀█ █▀▄▀█ █▀▀█ █▀▀█ █──█ 　 █▀▀█ █▀▀ █▀▀█ █▀▀█ █▀▀█ ▀▀█▀▀ 　 █▀▀█ █▀▀
# ─▀▀▀▄▄ █──█ █─▀─█ █─▀─█ █▄▄█ █▄▄▀ █▄▄█ 　 █▄▄▀ █▀▀ █──█ █──█ █▄▄▀ ──█── 　 █──█ █▀▀
# ░█▄▄▄█ ─▀▀▀ ▀───▀ ▀───▀ ▀──▀ ▀─▀▀ ▄▄▄█ 　 ▀─▀▀ ▀▀▀ █▀▀▀ ▀▀▀▀ ▀─▀▀ ──▀── 　 ▀▀▀▀ ▀──

# █▀▀ █▀▀ █── █▀▀ █▀▀ ▀▀█▀▀ █▀▀ █▀▀▄ 　 █▀▀▀ █▀▀ █▀▀▄ █▀▀ █▀▀ 　 █▀▀█ █▀▀▄ █▀▀▄
# ▀▀█ █▀▀ █── █▀▀ █── ──█── █▀▀ █──█ 　 █─▀█ █▀▀ █──█ █▀▀ ▀▀█ 　 █▄▄█ █──█ █──█
# ▀▀▀ ▀▀▀ ▀▀▀ ▀▀▀ ▀▀▀ ──▀── ▀▀▀ ▀▀▀─ 　 ▀▀▀▀ ▀▀▀ ▀──▀ ▀▀▀ ▀▀▀ 　 ▀──▀ ▀──▀ ▀▀▀─

# █▀▀▄ ─▀─ █▀▀ █▀▀ █▀▀█ █▀▀ █▀▀ █▀▀
# █──█ ▀█▀ ▀▀█ █▀▀ █▄▄█ ▀▀█ █▀▀ ▀▀█
# ▀▀▀─ ▀▀▀ ▀▀▀ ▀▀▀ ▀──▀ ▀▀▀ ▀▀▀ ▀▀▀

summary_df = snp_within_genes.copy()
summary_df.drop("Is SNP's near or within a gene", axis=1, inplace=True)

summary_df["Is reported on HuGE"] = np.where(
    snp_within_genes["Gene Id"].isin(genes_disease_info["gene_id"]), "Y", "N",
)

genes_disease_kegg_info_found_lookup = genes_disease_kegg_info[
    ~genes_disease_kegg_info["Disease"].isnull()
]
summary_df["Is reported on KEGG"] = np.where(
    snp_within_genes["Gene Id"].isin(genes_disease_kegg_info_found_lookup["Gene Id"]),
    "Y",
    "N",
)
pathway_kegg_genes_with_info_found_lookup = pathway_kegg_genes_with_info[
    ~pathway_kegg_genes_with_info["Disease"].isnull()
]
summary_df["Is reported on KEGG pathways"] = np.where(
    snp_within_genes["Gene Id"].isin(
        pathway_kegg_genes_with_info_found_lookup["Gene Id"]
    ),
    "Y",
    "N",
)

summary_df["Is reported"] = np.where(
    snp_within_genes["Gene Id"].isin(genes_disease_info["gene_id"])
    | snp_within_genes["Gene Id"].isin(genes_disease_kegg_info_found_lookup["Gene Id"])
    | snp_within_genes["Gene Id"].isin(
        pathway_kegg_genes_with_info_found_lookup["Gene Id"]
    ),
    "Y",
    "N",
)


# ░█─── ─▀─ █▀▀ ▀▀█▀▀ 　 █▀▀█ █▀▀ 　 █▀▀▄ █▀▀█ ▀▀█▀▀ █▀▀█ 　 ▀▀█▀▀ █▀▀█ █▀▀▄ █── █▀▀
# ░█─── ▀█▀ ▀▀█ ──█── 　 █──█ █▀▀ 　 █──█ █▄▄█ ──█── █▄▄█ 　 ──█── █▄▄█ █▀▀▄ █── █▀▀
# ░█▄▄█ ▀▀▀ ▀▀▀ ──▀── 　 ▀▀▀▀ ▀── 　 ▀▀▀─ ▀──▀ ──▀── ▀──▀ 　 ──▀── ▀──▀ ▀▀▀─ ▀▀▀ ▀▀▀

with st.expander("List of genes selected that search on NCBI 🔍"):
    AgGrid(
        st.session_state.gene_nbci_info,
        reload_data=True,
        gridOptions=get_base_grid_options(st.session_state.gene_nbci_info),
        enable_enterprise_modules=True,
    )

st.subheader("Mapping from selected to dataset 🔎")
AgGrid(
    gene_ncbi_info_repr,
    reload_data=True,
    gridOptions=get_base_grid_options(gene_ncbi_info_repr, "Probe Set ID"),
    enable_enterprise_modules=True,
    key="gene_ncbi_info_repr",
)

st.subheader("Special case: SNPs can be found inside and outside the genes")
AgGrid(
    gene_ncbi_info_snp_in_out,
    reload_data=True,
    gridOptions=get_base_grid_options(gene_ncbi_info_snp_in_out, "Probe Set ID"),
    enable_enterprise_modules=True,
    key="gene_ncbi_info_snp_in_out",
)

st.subheader("Summary report of selected genes and diseases 📊")
AgGrid(
    summary_df,
    reload_data=True,
    gridOptions=get_base_grid_options(summary_df, ["Probe Set ID"]),
    enable_enterprise_modules=True,
)

st.subheader("Gene that search on Disease (HuGE) 😎")
AgGrid(
    disease_huge_repr,
    reload_data=True,
    gridOptions=get_base_grid_options(disease_huge_repr, ["Disease"]),
    enable_enterprise_modules=True,
)

st.subheader("Gene that search on Disease (KEGG) 🤧")
AgGrid(
    genes_disease_kegg_info,
    reload_data=True,
    gridOptions=get_base_grid_options(genes_disease_kegg_info, "Disease"),
    enable_enterprise_modules=True,
    key="search_disease_kegg_info",
)

st.subheader("Gene that search on Pathway (KEGG) 🪶")
AgGrid(
    pathway_kegg_genes_with_info,
    reload_data=True,
    gridOptions=get_base_grid_options(pathway_kegg_genes_with_info, "Disease"),
    enable_enterprise_modules=True,
)


# ░█──░█ ─▀─ █▀▀ █───█ 　 █▀▀█ █▀▀█ █───█ 　 █▀▀▄ █▀▀█ ▀▀█▀▀ █▀▀█
# ─░█░█─ ▀█▀ █▀▀ █▄█▄█ 　 █▄▄▀ █▄▄█ █▄█▄█ 　 █──█ █▄▄█ ──█── █▄▄█
# ──▀▄▀─ ▀▀▀ ▀▀▀ ─▀─▀─ 　 ▀─▀▀ ▀──▀ ─▀─▀─ 　 ▀▀▀─ ▀──▀ ──▀── ▀──▀

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
