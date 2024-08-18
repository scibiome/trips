import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl


def plot_grn_trips(G_all, all_tfs=[], domino_genes=[], figsize=(15, 15),
                   layout = "graphviz_layout", prog="sfdp",
                   tf_color="lightsteelblue", domino_color="sandybrown", both_color="lightgreen",
                   normal_color="lightgray", highlight_mode="tf_and_module", fontsize_legend=12,
                   nodesize=2200, fontsize=14, legend_loc="upper left", filename=None, dpi=300):

    """
    Plots the GRN using networkx
    """

    # plt.rcParams['figure.facecolor'] = 'black'
    mpl.rcParams['figure.dpi'] = dpi
    plt.style.use('default')

    G = G_all.copy()
    
    # Remove dashes in the node names, replace with underscores
    mapping = {k:k.replace("-","_") for k in G.nodes()}
    G = nx.relabel_nodes(G, mapping)

    if highlight_mode == "tf_and_module":
        common_genes = list(set(all_tfs).intersection(set(domino_genes)))

    node_colors = []
    for node in G.nodes():
        if node in common_genes:
            node_colors.append(both_color)
        elif node in all_tfs:
            node_colors.append(tf_color)
        elif node in domino_genes:
            node_colors.append(domino_color)
        else:
            node_colors.append(normal_color)
            
    fig, ax = plt.subplots(figsize=figsize)

    if layout == "spring_layout":
        pos = nx.spring_layout(G)
    
    elif layout == "random_layout":
        pos = nx.random_layout(G)

    elif layout == 'graphviz_layout':
        # https://networkx.org/documentation/stable/reference/generated/networkx.drawing.nx_pydot.graphviz_layout.html
        pos = nx.nx_pydot.graphviz_layout(G, prog=prog)

    nx.draw_networkx_nodes(G, pos, node_size=nodesize, ax=ax, node_color=node_colors)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=fontsize)
    nx.draw_networkx_edges(G, pos, edge_color="black", arrows=True, arrowsize=20, width=2, min_target_margin=25, ax=ax)

    from matplotlib.lines import Line2D
    
    custom_lines = [Line2D([0], [0], color=tf_color, lw=8),
                    Line2D([0], [0], color=domino_color, lw=8),
                    Line2D([0], [0], color=both_color, lw=8)]
    ax.legend(custom_lines, ['TF', 'DOMINO module', 'Both TF and in DOMINO module'], fontsize=fontsize_legend, loc=legend_loc)
    ax= plt.gca()
    ax.axis('off')

    if filename:
        if filename.endswith("svg"):
            plt.savefig(filename, format="svg", bbox_inches="tight")
        else:
            plt.savefig(filename, format="png", bbox_inches="tight")
    plt.show()