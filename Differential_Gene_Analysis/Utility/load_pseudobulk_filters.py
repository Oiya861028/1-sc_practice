import pandas as pd
from pathlib import Path
from typing import List, Tuple, Union

def load_filters(merge = True) -> Union[List[pd.DataFrame], Tuple[List[pd.DataFrame], List[pd.DataFrame]]]:
    '''
    Load the pre and post psedubolk filter present in the DataFrames folder. If merge is True, will return one list with both pre and post filter dataframe.
    Else will return both separately

    Parameters
    ----------
    merge: Whether to merge the pre and post filter into one List. Default to True

    '''
    pre_filter_path = Path("/project/imoskowitz/yubin/1-sc_practice/Differential_Gene_Analysis/DataFrames/Pre-pseudobulk-filter")
    post_filter_path = Path("/project/imoskowitz/yubin/1-sc_practice/Differential_Gene_Analysis/DataFrames/Post-pseudobulk-filter")
    pre_filter_path = sorted(pre_filter_path.iterdir())
    post_filter_path = sorted(post_filter_path.iterdir())
    pre_filters = [] # list of pre filter DEG ranging from genes present in 0-5 cells filter
    for file_path in pre_filter_path:
        if file_path.is_file():  
            pre_filters.append(pd.read_csv(file_path))

    post_filters = [] # list of post filter DEG ranging from genes present in 0-5 cells filter
    for file_path in post_filter_path:
        if file_path.is_file():  
            post_filters.append(pd.read_csv(file_path))
    
    if merge:
        return pre_filters + post_filters
    else:
        return pre_filters, post_filters

    