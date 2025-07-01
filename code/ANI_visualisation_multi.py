# Noma Metagenomics: Recovery and analysis of Treponema MAGs from noma samples - Figure 4C
# By: Angus M. O'Ferrall & Richard N. Goodman
# Shotgun metagenomic analysis of the oral microbiomes of children with noma reveals a novel disease-associated organism
# Michael Olaleye, Angus M O’Ferrall, Richard N Goodman, Deogracia W Kabila, Miriam Peters, Gregoire Falq, Joseph Samuel, Donal Doyle, Diana Gomez, Gbemisola Oloruntuyi, Shafiu Isah, Adeniyi S Adetunji, Elise N Farley, Nicholas J Evans, Mark Sherlock, Adam P Roberts, Mohana Amirtharajah, Stuart Ainsworth

# Load packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def create_multi_ani_heatmaps_with_rotation(files, output_file_png, output_file_svg, rotation=90, box_size=1, annotation_fontsize=20, axis_fontsize=20, colorbar_fontsize=20):
    """
    Creates a figure with multiple ANI heatmaps, ensuring proper rotation, constant box size, and consistent formatting.

    Parameters:
        files (list): List of file paths to the Excel files containing the ANI tables.
        output_file_png (str): Path to save the combined heatmap as a PNG file.
        output_file_svg (str): Path to save the combined heatmap as an SVG file.
        rotation (int): Rotation angle for the heatmap (0, 90, 180, 270). Default is 90.
        box_size (float): Fixed size of each square in the heatmap (inches).
        annotation_fontsize (int): Font size for annotations inside the heatmaps.
        axis_fontsize (int): Font size for axis labels and ticks.
        colorbar_fontsize (int): Font size for the color bar label and ticks.
    """
    # Load all ANI tables and determine the largest matrix dimensions
    ani_tables = []
    for file in files:
        table = pd.read_excel(file, index_col=0).apply(pd.to_numeric, errors='coerce').fillna(0)
        
        # Apply rotation to each matrix
        if rotation == 90:
            table = table.transpose()
        elif rotation == 180:
            table = table.iloc[::-1, ::-1]
        elif rotation == 270:
            table = table.transpose().iloc[::-1, ::-1]
        elif rotation not in [0, 90, 180, 270]:
            raise ValueError("Rotation must be 0, 90, 180, or 270 degrees.")
        
        ani_tables.append(table)

    # Find the largest matrix dimensions to calculate figure size
    max_rows = max([table.shape[0] for table in ani_tables])
    max_cols = max([table.shape[1] for table in ani_tables])

    # Calculate the figure dimensions based on the largest matrix
    fig_width = len(files) * max_cols * box_size  # Total width for all heatmaps in a single row
    fig_height = max_rows * box_size  # Total height based on the largest matrix

    # Create the figure and axes
    fig, axes = plt.subplots(1, len(files), figsize=(fig_width, fig_height), constrained_layout=True)
    if len(files) == 1:
        axes = [axes]  # Ensure axes is iterable even for one heatmap

    # Set a shared color scale for all heatmaps
    vmin, vmax = 95, 100
    cmap = "coolwarm"

    # Generate individual heatmaps
    for i, (ani_table, ax) in enumerate(zip(ani_tables, axes)):
        # Generate a mask for the upper triangle
        mask = np.triu(np.ones_like(ani_table, dtype=bool))
        if rotation == 90 or rotation == 270:
            mask = np.transpose(mask)
        elif rotation == 180:
            mask = np.flip(mask)

        heatmap = sns.heatmap(
            ani_table,
            annot=True,
            fmt=".1f",
            mask=mask,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            annot_kws={"size": annotation_fontsize},
            xticklabels=ani_table.columns,
            yticklabels=ani_table.index,
            ax=ax,
            cbar=i == len(files) - 1,  # Add colorbar only for the last heatmap
            cbar_kws={"label": "ANI (%)", "shrink": 0.8} if i == len(files) - 1 else None,
        )

        # Adjust color bar font size
        if i == len(files) - 1:  # Only for the last heatmap
            cbar = heatmap.collections[0].colorbar
            cbar.ax.tick_params(labelsize=colorbar_fontsize)  # Adjust tick font size
            cbar.set_label("ANI (%)", fontsize=colorbar_fontsize)  # Adjust label font size

        # Set axis font sizes and rotation
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=axis_fontsize, rotation=90, ha="center")
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=axis_fontsize)
        ax.set_title(f"Matrix {i + 1} (Rot: {rotation}°)", fontsize=axis_fontsize + 4)

    # Save the combined figure
    plt.savefig(output_file_png, dpi=600)
    plt.savefig(output_file_svg, dpi=600)
    plt.close()

    print(f"Combined heatmap saved to {output_file_png}")
    print(f"Combined heatmap saved to {output_file_svg}")

# Example Usage
files = ["../data/ANI_A.xlsx", "../data/ANI_B.xlsx", "../data/ANI_C.xlsx", "../data/ANI_D.xlsx", "../data/ANI_E.xlsx", "../data/ANI_F.xlsx"]
create_multi_ani_heatmaps_with_rotation(
    files=files,
    output_file_png="../imgs/Figure_4C.png",
    output_file_svg="../imgs/Figure_4C.svg",
    rotation=90,  # Rotate all heatmaps by 90 degrees
    box_size=1,  # Uniform box size
    colorbar_fontsize=16  # Bigger font for the color bar
)
