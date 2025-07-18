import pandas as pd

from src.processing.ncbi_info import merge_gene_symbol


def test_official_symbol_present_alias_none():
    row = pd.Series({
        "Official Symbol": "ABC",
        "Gene symbol": "GENE",
        "Also known as": None,
    })
    assert merge_gene_symbol(row) == ["ABC"]


def test_official_symbol_absent_alias_none():
    row = pd.Series({
        "Official Symbol": pd.NA,
        "Gene symbol": "GENE",
        "Also known as": None,
    })
    assert merge_gene_symbol(row) == ["GENE"]


def test_official_symbol_present_alias_list():
    row = pd.Series({
        "Official Symbol": "ABC",
        "Gene symbol": "GENE",
        "Also known as": ["alias1", "alias2"],
    })
    assert merge_gene_symbol(row) == ["ABC", "alias1", "alias2"]


def test_official_symbol_absent_alias_list():
    row = pd.Series({
        "Official Symbol": pd.NA,
        "Gene symbol": "GENE",
        "Also known as": ["alias1"],
    })
    assert merge_gene_symbol(row) == ["GENE", "alias1"]
