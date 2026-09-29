import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches

# 1. LOAD AND CLEAN DATA
df = pd.read_csv('merged_and_cleaned_16S_data.csv')
df_filtered = df[df['tribe'] != 'na'].copy()

# Identify genus columns
meta_cols = ['sample-id', 'latitude', 'longitude', 'sample_type', 'host_species', 
             'host_genus', 'preservation', 'state', 'tribe', 'symbiont', 
             'native', 'habitat', 'sequencingbatch']
genus_cols = [col for col in df_filtered.columns if col not in meta_cols]

# Identify Top 20 Genera
genus_sums = df_filtered[genus_cols].sum().sort_values(ascending=False)
top_20 = genus_sums.head(20).index.tolist()
tribes = sorted(df_filtered['tribe'].unique().tolist())

# Prepare Link Data
link_data = []
for t in tribes:
    tdf = df_filtered[df_filtered['tribe'] == t]
    for g in top_20:
        val = float(tdf[g].sum())
        if val > 0:
            link_data.append({'tribe': t, 'genus': g, 'value': val})
links_df = pd.DataFrame(link_data)

# 2. BALANCED VERTICAL SCALING
total_flow = links_df['value'].sum()
BAR_TOTAL_H = 0.75  # Proportion of height dedicated to data bars

def get_balanced_positions(nodes, group_col):
    pos = {}
    n = len(nodes)
    node_heights = []
    for node in nodes:
        node_sum = links_df[links_df[group_col] == node]['value'].sum()
        # Scale bar height relative to its abundance
        h = (node_sum / total_flow) * BAR_TOTAL_H
        h = max(h, 0.002) # Minimum height for visibility
        node_heights.append(h)
    
    # Calculate unique gap for this column to fill exactly [0, 1]
    actual_bar_h = sum(node_heights)
    gap = (1.0 - actual_bar_h) / (n - 1)
    
    curr_y = 1.0
    for i, node in enumerate(nodes):
        h = node_heights[i]
        pos[node] = (curr_y - h, curr_y)
        curr_y -= (h + gap)
    return pos

tribe_pos = get_balanced_positions(tribes, 'tribe')
genus_pos = get_balanced_positions(top_20, 'genus')

# 3. COLOR PALETTE
TRIBE_PALETTE = ['#E63946','#F4A261','#E9C46A','#2A9D8F','#264653',
                 '#A8DADC','#457B9D','#1D3557','#6D6875','#B5838D',
                 '#E8998D','#F7B267','#F4845F','#F25C54','#A23B72']
tribe_colors = {t: TRIBE_PALETTE[i % len(TRIBE_PALETTE)] for i, t in enumerate(tribes)}

# 4. PLOTTING
fig, ax = plt.subplots(figsize=(22, 28))
ax.set_xlim(-0.25, 1.25)
ax.set_ylim(-0.02, 1.05)
ax.axis('off')

tribe_curr_y = {t: pos[1] for t, pos in tribe_pos.items()}
genus_curr_y = {g: pos[1] for g, pos in genus_pos.items()}

# Draw Flow Ribbons
for _, row in links_df.sort_values(['tribe', 'value'], ascending=[True, False]).iterrows():
    t, g, val = row['tribe'], row['genus'], row['value']
    
    # Calculate height on source side (Tribe)
    t_total = links_df[links_df['tribe'] == t]['value'].sum()
    t_node_h = tribe_pos[t][1] - tribe_pos[t][0]
    h_tribe_side = (val / t_total) * t_node_h
    
    # Calculate height on target side (Genus)
    g_total = links_df[links_df['genus'] == g]['value'].sum()
    g_node_h = genus_pos[g][1] - genus_pos[g][0]
    h_genus_side = (val / g_total) * g_node_h
    
    y_s_top, y_s_bot = tribe_curr_y[t], tribe_curr_y[t] - h_tribe_side
    tribe_curr_y[t] = y_s_bot
    
    y_e_top, y_e_bot = genus_curr_y[g], genus_curr_y[g] - h_genus_side
    genus_curr_y[g] = y_e_bot
    
    # Draw Ribbon using Curve4
    x = [0.05, 0.4, 0.6, 0.95]
    verts = [(x[0], y_s_top), (x[1], y_s_top), (x[2], y_e_top), (x[3], y_e_top),
             (x[3], y_e_bot), (x[2], y_e_bot), (x[1], y_s_bot), (x[0], y_s_bot)]
    codes = [mpath.Path.MOVETO, mpath.Path.CURVE4, mpath.Path.CURVE4, mpath.Path.CURVE4,
             mpath.Path.LINETO, mpath.Path.CURVE4, mpath.Path.CURVE4, mpath.Path.CURVE4]
    
    path = mpath.Path(verts + [verts[0]], codes + [mpath.Path.CLOSEPOLY])
    ax.add_patch(mpatches.PathPatch(path, facecolor=tribe_colors[t], alpha=0.35, edgecolor='none'))

# Draw Bars and Labels
for t, (y_bot, y_top) in tribe_pos.items():
    ax.add_patch(mpatches.Rectangle((0, y_bot), 0.05, y_top-y_bot, facecolor=tribe_colors[t]))
    ax.text(-0.02, (y_bot+y_top)/2, t, ha='right', va='center', fontsize=20, fontweight='bold')

for g, (y_bot, y_top) in genus_pos.items():
    ax.add_patch(mpatches.Rectangle((0.95, y_bot), 0.05, y_top-y_bot, facecolor='#444444'))
    ax.text(1.02, (y_bot+y_top)/2, g, ha='left', va='center', fontsize=20, fontstyle='italic', fontweight='bold')

plt.savefig('sankey_figure_bacteria.png', dpi=300, bbox_inches='tight')
plt.savefig('sankey_figure_bacteria.svg', format='svg', bbox_inches='tight')
