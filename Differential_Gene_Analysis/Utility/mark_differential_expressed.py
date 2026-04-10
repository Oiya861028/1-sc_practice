"""
Utility functions for marking differentially expressed genes in pseudobulk analysis results.
"""

import pandas as pd
from typing import Optional


def generate_differential_expressed_column(
    df: pd.DataFrame,
    index: Optional[int] = None,
    lfc_threshold: float = 0.5,
    padj_threshold: float = 0.05
) -> pd.DataFrame:
    """
    Add columns to mark differentially expressed genes based on log2 fold change and adjusted p-value thresholds.

    Args:
        df: DataFrame containing differential expression results with columns 'padj' and 'log2FoldChange'
        index: Optional index for processing specific rows (currently unused)
        lfc_threshold: Log2 fold change threshold for significance (default: 0.5)
        padj_threshold: Adjusted p-value threshold for significance (default: 0.05)

    Returns:
        DataFrame with additional columns:
        - 'Significance': Boolean mask indicating significant genes
        - 'diff_expr': String labels ('up', 'down', or 'NS' for not significant)

    Raises:
        KeyError: If required columns 'padj' or 'log2FoldChange' are missing from the DataFrame
    """
    # Validate required columns exist
    required_cols = ['padj', 'log2FoldChange']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise KeyError(f"Required columns missing from DataFrame: {missing_cols}")

    # Create significance mask based on thresholds
    significance_mask = (df["padj"] < padj_threshold) & (abs(df["log2FoldChange"]) > lfc_threshold)
    df = df.copy()  # Avoid modifying the original DataFrame
    df["Significance"] = significance_mask

    # Add categorical differential expression labels
    df["diff_expr"] = "NS"  # Default to not significant

    # Mark upregulated genes (significant and positive log2FC)
    up_mask = significance_mask & (df["log2FoldChange"] > lfc_threshold)
    df.loc[up_mask, "diff_expr"] = "up"

    # Mark downregulated genes (significant and negative log2FC)
    down_mask = significance_mask & (df["log2FoldChange"] < -lfc_threshold)
    df.loc[down_mask, "diff_expr"] = "down"

    return df