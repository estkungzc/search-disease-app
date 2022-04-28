from typing import List, Union
from pandas import DataFrame
from st_aggrid.grid_options_builder import GridOptionsBuilder


def get_base_grid_options(df: DataFrame, columnGroup: Union[List, str] = None):
    gb = GridOptionsBuilder.from_dataframe(df)
    if columnGroup:
        if type(columnGroup) == str:
            columnGroup = [columnGroup]
        for col in columnGroup:
            gb.configure_column(col, rowGroup=True, aggFunc="sum")
    gb.configure_pagination()
    gb.configure_side_bar()
    gb.configure_default_column(
        groupable=True, value=True, enableRowGroup=True, aggFunc="sum"
    )
    return gb.build()
