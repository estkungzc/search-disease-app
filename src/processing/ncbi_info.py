from pandas import DataFrame
import pandas as pd


def merge_gene_symbol(row: DataFrame) -> list:
    official_symbol = (
        row["Official Symbol"]
        if not pd.isna(row["Official Symbol"])
        else row["Gene symbol"]
    )
    usage_symbol = [official_symbol] if official_symbol else []

    also_known_as = row["Also known as"] if row["Also known as"] is not None else []
    return usage_symbol + also_known_as
