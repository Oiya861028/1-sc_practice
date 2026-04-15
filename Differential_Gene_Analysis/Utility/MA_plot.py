import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text


def custom_MA_plot(
        stats_df,
        ax,
        x="baseMean",
        y="log2FoldChange",
        n_top_genes=10,
        point_size=5,
        alpha_threshold=0.05,
        lfc_threshold=0.5,
        custom_degs_list=None,) -> None:
    """
    Create an MA plot (Mean vs Log Fold Change) with labeled top genes.

    Parameters:
    - stats_df: DataFrame with differential expression results
    - ax: matplotlib axis to plot on
    - x: column name for x-axis (mean expression)
    - y: column name for y-axis (log fold change)
    - n_top_genes: number of top genes to label from up/down regulated
    - point_size: size of scatter points
    - alpha_threshold: significance threshold for p-adjusted values
    - lfc_threshold: log fold change threshold for significance
    - custom_degs_list: optional list of genes to highlight
    """

    # Ensure we have the required columns
    if x not in stats_df.columns or y not in stats_df.columns:
        raise ValueError(f"DataFrame must contain columns '{x}' and '{y}'")

    # Create significance column if not present
    if 'Significance' not in stats_df.columns:
        if 'padj' in stats_df.columns:
            stats_df['Significance'] = (stats_df['padj'] < alpha_threshold) & (abs(stats_df[y]) > lfc_threshold)
        else:
            stats_df['Significance'] = abs(stats_df[y]) > lfc_threshold

    # Create regulation column if not present
    if 'diff_expr' not in stats_df.columns:
        stats_df['diff_expr'] = 'NS'
        stats_df.loc[stats_df['Significance'] & (stats_df[y] > lfc_threshold), 'diff_expr'] = 'up'
        stats_df.loc[stats_df['Significance'] & (stats_df[y] < -lfc_threshold), 'diff_expr'] = 'down'

    # Plot all points
    colors = {'NS': 'gray', 'up': 'red', 'down': 'blue'}
    for category, color in colors.items():
        subset = stats_df[stats_df['diff_expr'] == category]
        ax.scatter(subset[x], subset[y], c=color, s=point_size, alpha=0.6, label=category)

    # Select top genes to label
    up_genes = stats_df[(stats_df['diff_expr'] == 'up')].nlargest(n_top_genes, y)
    down_genes = stats_df[(stats_df['diff_expr'] == 'down')].nsmallest(n_top_genes, y)

    # Combine top genes
    top_genes = pd.concat([up_genes, down_genes])

    # Add custom genes if provided
    if custom_degs_list is not None:
        custom_subset = stats_df[stats_df.index.isin(custom_degs_list)]
        top_genes = pd.concat([top_genes, custom_subset]).drop_duplicates()

    # Label the genes
    texts = []
    for idx, row in top_genes.iterrows():
        gene_name = row.get('Unnamed: 0', str(idx))
        texts.append(ax.text(row[x], row[y], gene_name, fontsize=8, ha='center', va='center'))

    # Adjust text to avoid overlaps
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

    # Add threshold lines
    ax.axhline(y=lfc_threshold, color='red', linestyle='--', alpha=0.7, linewidth=1)
    ax.axhline(y=-lfc_threshold, color='red', linestyle='--', alpha=0.7, linewidth=1)

    # Set labels and title
    ax.set_xlabel('Mean Expression (baseMean)')
    ax.set_ylabel('Log2 Fold Change')
    ax.set_title('MA Plot')

    # Use log scale for x-axis if values span large range
    if stats_df[x].max() / stats_df[x].min() > 100:
        ax.set_xscale('log')

    # Add legend
    ax.legend(title='Regulation', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Add grid
    ax.grid(True, alpha=0.3)
