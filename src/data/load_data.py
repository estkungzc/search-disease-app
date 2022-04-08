import streamlit as st
import pandas as pd

from pathlib import Path

data_folder = Path(__file__).parent.parent.parent / "data"


@st.cache(persist=True)
def load_mapping_dataset():
    return pd.read_csv(data_folder / "mapping-dataset.csv")


@st.cache(persist=True)
def load_disease_from_huge():
    return pd.read_csv(data_folder / "disease-mapping-huge.csv")


@st.cache(persist=True)
def load_disease_from_kegg():
    df = pd.read_csv(data_folder / "disease-mapping-kegg.csv")
    df["gene_id"] = df["gene_id"].astype(str)
    return df


@st.cache(persist=True)
def load_pathway_from_kegg():
    df = pd.read_csv(data_folder / "kegg-pathway.csv")
    df["hsa"] = df["hsa"].astype(str)
    return df


@st.cache(persist=True)
def load_disease_pathway_from_kegg():
    return pd.read_csv(data_folder / "kegg-pathway-disease.csv")
