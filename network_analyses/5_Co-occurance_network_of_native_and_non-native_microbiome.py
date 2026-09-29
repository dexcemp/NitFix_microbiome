import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

#  1. LOAD DATA 
file_path = 'plant_bacteria_and_fungi_data_merged.csv'   # place CSV in the same folder as this script
df = pd.read_csv(file_path)

# 2. FILTER: rhizobial samples only 
legume_df = df[df['symbiont'] == 'rhizobial'].copy()

#  3. IDENTIFY MICROBIAL COLUMNS 
metadata_cols = [
    'sample-id', 'latitude', 'longitude', 'sample_type', 'host_species',
    'host_genus', 'preservation', 'state', 'tribe', 'symbiont',
    'native', 'habitat', 'sequencingbatch', 'cycle_number', 'annealtempdegrees'
]

all_cols = list(df.columns)
cycle_idx = all_cols.index('cycle_number')

bacteria_pool = [c for c in all_cols[:cycle_idx] if c not in metadata_cols]
fungi_pool    = [c for c in all_cols[cycle_idx:] if c not in metadata_cols]

#  4. NETWORK BUILDER 
def build_and_plot(subset_df, label, title, threshold=0.25):
    """Build and save a co-occurrence network for a given subset."""

    # Select top 30 bacteria and top 30 fungi by mean relative abundance
    bact_data = subset_df[bacteria_pool].apply(pd.to_numeric, errors='coerce')
    fung_data = subset_df[fungi_pool].apply(pd.to_numeric, errors='coerce')

    top_30_bact = bact_data.mean().sort_values(ascending=False).head(30)
    top_30_fung = fung_data.mean().sort_values(ascending=False).head(30)

    top_60_genera = top_30_bact.index.tolist() + top_30_fung.index.tolist()
    top_60_values = pd.concat([top_30_bact, top_30_fung])

    # Spearman correlation matrix across all 60 genera
    corr_data = subset_df[top_60_genera].apply(pd.to_numeric, errors='coerce')
    corr_matrix = corr_data.corr(method='spearman')

    #  Build graph 
    G = nx.Graph()

    for genus in top_60_genera:
        node_type = 'Bacteria' if genus in bacteria_pool else 'Fungi'
        G.add_node(genus, abundance=top_60_values[genus], type=node_type)

    for i in range(len(top_60_genera)):
        for j in range(i + 1, len(top_60_genera)):
            corr = corr_matrix.iloc[i, j]
            if abs(corr) >= threshold:
                G.add_edge(top_60_genera[i], top_60_genera[j], weight=corr)

    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    print(f"[{label}] n samples={len(subset_df)}, nodes={n_nodes}, edges={n_edges}")

    #  Layout & colours 
    pos = nx.circular_layout(G, scale=10)

    color_map   = {'Bacteria': '#1f77b4', 'Fungi': '#ff7f0e'}
    node_colors = [color_map[G.nodes[n]['type']] for n in G.nodes()]

    node_abundances = [G.nodes[n]['abundance'] for n in G.nodes()]
    min_a, max_a    = min(node_abundances), max(node_abundances)
    node_sizes = [
        ((a - min_a) / (max_a - min_a + 1e-9)) * 9000 + 2000
        for a in node_abundances
    ]

    edges_data  = list(G.edges(data=True))
    edge_colors = ['#d62728' if d['weight'] < 0 else '#2ca02c' for _, _, d in edges_data]
    edge_widths = [abs(d['weight']) * 12 for _, _, d in edges_data]

    #  Draw 
    fig, ax = plt.subplots(figsize=(45, 45))

    nx.draw_networkx_edges(
        G, pos,
        edge_color=edge_colors, width=edge_widths, alpha=0.25, ax=ax
    )
    nx.draw_networkx_nodes(
        G, pos,
        node_size=node_sizes, node_color=node_colors,
        alpha=0.85, edgecolors='black', linewidths=2.5, ax=ax
    )

    # Radial labels
    label_radius = 10.7
    for node in G.nodes():
        x, y      = pos[node]
        angle_deg = np.degrees(np.arctan2(y, x))
        if angle_deg > 90 or angle_deg < -90:
            rotation = angle_deg + 180
            ha       = 'right'
        else:
            rotation = angle_deg
            ha       = 'left'
        ax.text(
            x / 10 * label_radius, y / 10 * label_radius, node,
            rotation=rotation, rotation_mode='anchor',
            fontsize=26, fontweight='bold', va='center', ha=ha
        )

    #Legend: manually placed in lower-right in data coordinates 
    lx, ly   = 12.8, -11.8   # bottom-right anchor
    row_h    = 1.05           # row spacing (builds upward)
    dot_r    = 0.38           # circle radius
    line_len = 0.45           # half-length of line symbol (short)
    fs       = 26

    legend_items = [
        ('line',   '#d62728', f'Negative correlation (r \u2264 \u2212{threshold})'),
        ('line',   '#2ca02c', f'Positive correlation (r \u2265 {threshold})'),
        ('circle', '#ff7f0e', 'Fungi'),
        ('circle', '#1f77b4', 'Bacteria'),
    ]

    for i, (kind, color, lbl) in enumerate(legend_items):
        y_row = ly + i * row_h
        sym_x = lx - 1.5
        txt_x = lx - 0.7
        if kind == 'circle':
            ax.add_patch(plt.Circle((sym_x, y_row), dot_r, color=color,
                                    zorder=5, clip_on=False))
        else:
            ax.plot([sym_x - line_len, sym_x + line_len], [y_row, y_row],
                    color=color, lw=6, solid_capstyle='round',
                    zorder=5, clip_on=False)
        ax.text(txt_x, y_row, lbl, fontsize=fs, va='center', ha='left',
                fontweight='bold', color='#222222', clip_on=False)

    plt.axis('off')
    ax.set_aspect('equal')
    plt.xlim(-13, 13)
    plt.ylim(-12, 12)

    # Save 
    out_base = f'cooccurrence_network_{label}'
    plt.savefig(f'{out_base}.png', bbox_inches='tight', dpi=300)
    plt.savefig(f'{out_base}.svg', bbox_inches='tight', format='svg')
    plt.close()
    print(f"   saved {out_base}.png / .svg")


#  5. RUN FOR EACH GROUP 
native_df    = legume_df[legume_df['native'] == 'yes'].copy()
nonnative_df = legume_df[legume_df['native'] == 'no'].copy()

build_and_plot(
    native_df,
    label='native_plants',
    title='Bacteria–Fungi Co-occurrence Network\n(Native Plants, Rhizobial Samples)'
)

build_and_plot(
    nonnative_df,
    label='nonnative_plants',
    title='Bacteria–Fungi Co-occurrence Network\n(Non-Native Plants, Rhizobial Samples)'
)

print("Done")
