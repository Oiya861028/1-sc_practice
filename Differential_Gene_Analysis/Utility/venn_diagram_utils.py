"""
Utility functions for Venn diagram analysis of differential gene expression data.

This module provides functions for:
- Loading and processing filter result dataframes
- Extracting gene sets by expression type
- Creating labeled Venn diagrams
- Plotting filter size comparisons
"""

import pandas as pd
from pathlib import Path
import numpy as np
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
from typing import List, Tuple, Set, Optional


def create_labeled_venn3(
    sets: Tuple[Set[str], Set[str], Set[str]],
    labels: Tuple[str, str, str],
    title: str = "",
    figsize: Tuple[int, int] = (10, 8),
    label_regions: bool = True,
) -> None:
    """
    Create a Venn diagram with optional gene name labels for each region.

    Args:
        sets: Tuple of three gene sets to compare
        labels: Tuple of three labels for the sets
        title: Optional title for the plot
        figsize: Figure size as (width, height)
        label_regions: Whether to add gene names to each region
    """
    set1, set2, set3 = sets

    # Calculate gene intersections and unique sets
    center_intersection = set1 & set2 & set3  # All three sets
    set1_set2_only = (set1 & set2) - set3     # Set1 and Set2 only
    set1_set3_only = (set1 & set3) - set2     # Set1 and Set3 only
    set2_set3_only = (set2 & set3) - set1     # Set2 and Set3 only
    set1_only = set1 - set2 - set3            # Set1 only
    set2_only = set2 - set1 - set3            # Set2 only
    set3_only = set3 - set2 - set1            # Set3 only

    # Create figure
    plt.figure(figsize=figsize)
    if title:
        plt.title(title, fontsize=14, pad=20)

    # Create Venn diagram (hide default numbers)
    v = venn3(sets, labels, subset_label_formatter=lambda x: "")

    if label_regions:
        # Helper function to format gene lists
        def format_gene_list(genes: Set[str]) -> str:
            if not genes:
                return ""
            gene_list = sorted(list(genes))
            if len(gene_list) > 10:  # Limit to first 10 genes if too many
                gene_list = gene_list[:10] + [f"... +{len(genes)-10} more"]
            return '\n'.join(gene_list)

        # Label each region with gene names (only if not empty)
        region_labels = {
            '111': format_gene_list(center_intersection),
            '110': format_gene_list(set1_set2_only),
            '101': format_gene_list(set1_set3_only),
            '011': format_gene_list(set2_set3_only),
            '100': format_gene_list(set1_only),
            '010': format_gene_list(set2_only),
            '001': format_gene_list(set3_only)
        }

        # Apply labels only for non-empty regions
        for region_id, gene_text in region_labels.items():
            if gene_text.strip():  # Only label if there are genes
                v.get_label_by_id(region_id).set_text(gene_text)

    plt.tight_layout()
    plt.show()


def load_filter_dataframes(base_path: str | Path) -> List[pd.DataFrame]:
    """
    Load all CSV files from a directory containing filter results.

    Args:
        base_path: Path to directory containing filter result CSV files

    Returns:
        List of pandas DataFrames, one for each filter threshold
    """
    if type(base_path) == str:
        base_path = Path(base_path)
    if not Path(base_path).is_dir():
        return "Path is not a directory!"
    filter_files = sorted([f for f in base_path.iterdir() if f.is_file() and f.suffix == '.csv'])
    return [pd.read_csv(file_path) for file_path in filter_files]


def extract_gene_sets(filter_dfs: List[pd.DataFrame]) -> List[Set[str]]:
    """
    Extract sets of gene names from filter dataframes.

    Args:
        filter_dfs: List of filter result dataframes

    Returns:
        List of sets containing gene names for each filter
    """
    return [set(df.iloc[:, 1]) for df in filter_dfs]


def extract_differential_genes(
    filter_dfs: List[pd.DataFrame],
    expression_type: str
) -> List[Set[str]]:
    """
    Extract genes based on differential expression status.

    Args:
        filter_dfs: List of filter result dataframes
        expression_type: Type of expression to extract ('up', 'down', or 'all')

    Returns:
        List of sets containing gene names for each filter
    """
    if expression_type == 'all':
        return extract_gene_sets(filter_dfs)

    gene_sets = []
    for df in filter_dfs:
        if expression_type == 'up':
            mask = df['diff_expr'] == 'up'
        elif expression_type == 'down':
            mask = df['diff_expr'] == 'down'
        else:
            raise ValueError(f"Invalid expression_type: {expression_type}. Must be 'up', 'down', or 'all'")

        genes = set(df.loc[mask, df.columns[1]])
        gene_sets.append(genes)

    return gene_sets


def plot_filter_sizes(
    filter_sets_list: List[List[Set[str]]],
    labels: List[str],
    title: str = "Filter Size Comparison"
) -> None:
    """
    Plot scatter plots comparing the sizes of different filter sets.

    Args:
        filter_sets_list: List of lists of gene sets to compare
        labels: Labels for each filter set list
        title: Title for the plot
    """
    plt.figure(figsize=(10, 6))

    colors = ['red', 'blue', 'green', 'orange', 'purple']
    for i, (filter_sets, label) in enumerate(zip(filter_sets_list, labels)):
        sizes = [len(gene_set) for gene_set in filter_sets]
        plt.scatter(
            x=np.arange(len(sizes)),
            y=sizes,
            color=colors[i % len(colors)],
            alpha=0.7,
            label=label,
            s=50
        )

    plt.title(title)
    plt.xlabel("Filter Index")
    plt.ylabel("Number of Genes")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()