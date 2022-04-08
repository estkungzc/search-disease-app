from datetime import datetime
from multiprocessing.pool import ThreadPool
from numpy import datetime64, ndarray
import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
import re
import streamlit as st


def extract_gene_info(soup):
    # get updated date
    # pattern::: Gene ID: 440107, updated on 25-Jan-2022
    # pattern::: Gene ID: 338809, discontinued on 8-Jan-2020
    geneField = soup.find("span", class_="geneid")
    gene_id, status_date = geneField.text.split(", ")
    status, date = [s.strip() for s in status_date.split(" on")]

    # results
    usages = {"Status": status, "Date": date, "Replaced with Gene ID": None}

    # New pattern::: This record was replaced with Gene ID: 440107
    replaced = soup.find(
        lambda tag: tag.name == "span" and "This record was replaced with" in tag.text
    )
    if replaced:
        r = re.findall("\d+", replaced.text)
        new_gene_id = r[-1] if len(r) > 0 else None
        usages["Replaced with Gene ID"] = new_gene_id

    # get summary details
    results = soup.find("dl", id="summaryDl")

    # pre clean
    keys_list = []
    values_list = []
    for i in results.find_all("dt"):
        t = i.get_text(strip=True)
        t = re.sub(r"[\ \n]{2,}", " ", t)
        keys_list.append(t)
    for i in results.find_all("dd"):
        t = i.get_text(strip=True)
        t = re.sub(r"provided.\w*", "", t)
        values_list.append(t)

    # post clean -> structure, filter -> usages
    interest_keys = ["Official Symbol", "Also known as"]

    #     for k, v in zip(keys_list, values_list):
    #         if k in interest_keys:
    #             usages[k] = v

    for k, v in zip(keys_list, values_list):
        usages[k] = v

    usages["Also known as"] = (
        usages["Also known as"].split("; ") if "Also known as" in usages else None
    )
    return usages


def extract(gene_id: str):
    print(f"{gene_id} : start...")
    # request
    url = f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"
    page = requests.get(url)

    # get content
    soup = BeautifulSoup(page.content, "html.parser")

    # extract gene_info
    try:
        result = extract_gene_info(soup)
        result["Gene Id"] = str(gene_id)
        return result
    except:
        print("==========Occur an error on: ", gene_id)
        raise


NCBI_COLUMNS_INFO = [
    ("Status", str),
    ("Date", "datetime64[ns]"),
    ("Gene Id", int),
    ("Replaced with Gene ID", int),
    ("Official Symbol", str),
    ("Also known as", str),
    ("Primary source", str),
    ("See related", str),
    ("Gene type", str),
    ("RefSeq status", str),
    ("Organism", str),
    ("Lineage", str),
    ("Summary", str),
    ("Expression", str),
    ("Orthologs", str),
    ("NEW", str),
]


@st.cache(
    persist=True, allow_output_mutation=True, ttl=60 * 60 * 24, suppress_st_warning=True
)
def get_all_genes_info(selected_genes):
    # st.write("Cache miss: on call get_all_genes_info()")
    df = pd.DataFrame({c: pd.Series(dtype=t) for c, t in NCBI_COLUMNS_INFO})
    max_thread = min(len(selected_genes), 20)
    with ThreadPool(max_thread) as pool:
        for result in pool.map(extract, selected_genes):
            result["Date"] = datetime.strptime(result["Date"], "%d-%b-%Y").date()
            df = df.append(result, ignore_index=True)
    return df
