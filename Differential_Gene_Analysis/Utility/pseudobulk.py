import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def run_pseudobulk(subset: sc.AnnData, by: str, metadatas: str | list[str], min_cells: None | int) -> sc.AnnData:
    """
    Perform pseudobulk aggregation on a subset of AnnData, grouping by a specified condition.

    This function aggregates gene expression counts across cells within each unique group
    defined by the 'by' parameter, creating a new AnnData object with pseudobulk data.
    Metadata from the original subset is transferred to the resulting pseudobulk AnnData.
    Also, filter out genes present in less than specified cell if requested

    Parameters:
    -----------
    subset : sc.AnnData
        The AnnData object containing the subset of cells to perform pseudobulk on.
        It should include raw count data in the .X attribute.
    by : str
        The column name in subset.obs to group cells by (e.g., 'Sample').
        Each unique value in this column will result in one pseudobulk sample.
    metadatas : list[str]
        A list of column names from subset.obs to transfer as metadata to the pseudobulk AnnData.
        For each metadata column, the value from the first cell in the group is used.
    min_cells: None | int
        if int, then will filter genes with less than int cells expressed prior to pseudobulking
    Returns:
    --------
    sc.AnnData
        A new AnnData object containing the pseudobulk data, with aggregated gene counts
        and transferred metadata.
    """

    # Subset the subset to filter out genes
    if min_cells is not None:
        sc.pp.filter_genes(subset, min_cells = min_cells)
    # Initialize a list to hold pseudobulk AnnData objects for each group
    pbs = []

    # Iterate over each unique value in the 'by' column
    for sample in subset.obs[by].unique():
        # Subset the data to only include cells from the current group
        samp_cell_subset = subset[subset.obs[by] == sample]

        # Sum gene counts across all cells in the group (axis=0 for genes)
        count = samp_cell_subset.X.sum(axis=0)

        # Ensure count is in a compatible format for AnnData
        if isinstance(count, np.matrix):
            count = np.asarray(count)
        if sparse.issparse(count):
            count = count.A  # Convert sparse matrix to dense array if needed

        # Create a new AnnData object for this pseudobulk sample
        rep_adata = sc.AnnData(X=count, var=samp_cell_subset.var.copy())
        
        # Set observation names based on condition and replicate
        rep_adata.obs_names = [samp_cell_subset.obs["Condition"].iloc[0] + "_" + samp_cell_subset.obs["Replicate"].iloc[0]]
        
        # Add core metadata
        rep_adata.obs["Condition"] = samp_cell_subset.obs["Condition"].iloc[0]
        rep_adata.obs["Replicate"] = samp_cell_subset.obs["Replicate"].iloc[0]

        # Transfer additional metadata specified in metadatas
        if type(metadatas) == str:
            rep_adata.obs[metadatas] = samp_cell_subset.obs[metadatas].iloc[0]
        else:
            for metadata in metadatas:
                rep_adata.obs[metadata] = samp_cell_subset.obs[metadata].iloc[0]
        
        # Append the pseudobulk AnnData to the list
        pbs.append(rep_adata)
    
    # Concatenate all pseudobulk AnnData objects into a single object
    pb = sc.concat(pbs)
    return pb


# TODO: Run with independent_filter off 
def run_deseq(adata: sc.Adata, design: str | list[str], min_cells: int, independent_filter= True):
    """
    Run DESeq2 differential expression analysis on pseudobulk data.

    This function prepares the pseudobulk AnnData for DESeq2 by converting it to a DataFrame,
    creates a DESeqDataSet with the specified design formula, optionally filters genes based
    on minimum cell count, and runs the DESeq2 normalization and dispersion estimation.

    Parameters:
    -----------
    adata : sc.AnnData
        The pseudobulk AnnData object containing aggregated gene counts.
        The .X attribute should contain the count data.
    design : str or list[str]
        The design formula for DESeq2 (e.g., '~Condition'). Can be a string or a list of strings
        for more complex designs.
    min_cells : int
        The minimum number of cells required for a gene to be kept. If greater than 0,
        genes not expressed in at least this many samples will be filtered out.
    independent_filter: bool
        Whether to run the deseq with it's internal filter or not. Defaults to true

    Returns:
    --------
    DeseqDataSet
        The DESeqDataSet object after running DESeq2 normalization and dispersion estimation.
    """

    counts = pd.DataFrame(data = adata.X, columns = adata.var_names)

    dds = DeseqDataSet(counts = counts, metadata=adata.obs, design=design)

    if min_cells > 0:
        sc.pp.filter_genes(dds, min_cells=min_cells)
    
    dds.deseq2(independent_filter=independent_filter)
    return dds