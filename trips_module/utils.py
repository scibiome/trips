import pandas as pd
import numpy as np
import networkx as nx
import os
import random
import subprocess
import glob

flatten = lambda l: [item for sublist in l for item in sublist]


def get_degs(file_degs, lfc_thresh=1.0, pval_thresh=0.05):
    """
    file_degs: txt file containing
                the differential expression analysis results
                should have columns ["Log_FoldChange","AdjPValue"]
    gene_module: list of genes in a module (e.g. DOMINO module)
                genes in the module but not DEGs will be assigned
                the mean scores of all DEGs as prize
    Returns
    degs: list of DEGs satisfying the cutoffs
    """

    df = pd.read_csv(file_degs, sep="\t")
    df_sub = df[df["AdjPValue"] < pval_thresh]
    if df_sub.empty:
        raise Exception("No genes satisfiying p-value cutoff")

    df_sub = df_sub[abs(df_sub["Log_FoldChange"]) > lfc_thresh]
    if df_sub.empty:
        raise Exception("No genes satisfiying logFoldChange cutoff")

    degs = list(df_sub["Gene_symbol"])
    return degs


def get_grn(file_grn, directed=False, n_top_edges=None):

    df_grn = pd.read_csv(file_grn, sep="\t")

    if n_top_edges:
        df_grn = df_grn.head(n_top_edges)
    if directed:
        G_grn = nx.from_pandas_edgelist(df_grn, source="node1", target="node2", create_using=nx.DiGraph)
    else:
        G_grn = nx.from_pandas_edgelist(df_grn, source="node1", target="node2")

    return G_grn


def get_ppi_net(file_ppi_net):

    """
    :param file_ppi_net: Should be a tab-separated file with two columns ["node1", "node2"]
    :return: the networkx graph (PPI network)
    """

    df_ppi = pd.read_csv(file_ppi_net, sep="\t")
    ppi_network = nx.from_pandas_edgelist(df_ppi, source="node1", target="node2")
    return ppi_network


def convert_to_networkx_graph(df):
    """
    Input data frame of TF-TF edges
    Must have columns ["tf1","tf2"]
    Get the largest connected component
    """

    df["weight"] = 1

    G = nx.from_pandas_edgelist(df, "tf1", "tf2", edge_attr=["weight"])
    GG = G.subgraph(max(nx.connected_components(G), key=len))

    return GG


def pd_edges_from_stp_file(file_stp):

    with open(file_stp, 'r') as f:

        edge_tuples = []
        lines = f.readlines()
        for line in lines:

            if line.startswith("E "):
                edge = line.split()
                edge_tuples.append((edge[1], edge[2]))

        return pd.DataFrame(edge_tuples, columns=["from", "to"])


def get_unique_nodes(df):
    """

    :param df: should have two columns
    :return: list of unique nodes from data frame
    """
    df.columns = ["from", "to"]
    list1 = df["from"]
    list2 = df["to"]
    all_nodes = list(list1) + list(list2)
    all_nodes = list(set(all_nodes))

    return all_nodes


def get_grn(file_grn, directed=True, n_top_edges=None, weighted=False):

    df_grn = pd.read_csv(file_grn, sep="\t")
    print(df_grn.head())
    if df_grn.shape[1] > 2:
        df_grn = df_grn.sort_values(by=["weight"], ascending=False)
    if n_top_edges:
        df_grn = df_grn.head(n_top_edges)

    if (directed and weighted):
        G_grn = nx.from_pandas_edgelist(df_grn, source="from", target="to", edge_attr='weight', create_using=nx.DiGraph)
    elif (directed and not weighted):
        G_grn = nx.from_pandas_edgelist(df_grn, source="from", target="to", create_using=nx.DiGraph)
    elif (not directed and weighted):
        G_grn = nx.from_pandas_edgelist(df_grn, source="from", target="to", edge_attr='weight')
    elif (not directed and not weighted):
        G_grn = nx.from_pandas_edgelist(df_grn, source="from", target="to")

    return G_grn


def load_ppi_biogrid_custom(file_biogrid_raw, filter_type="scored"):

    df_biogrid = pd.read_csv(file_biogrid_raw, sep="\t")
    if filter_type == "scored":
        df_biogrid = df_biogrid[df_biogrid["Score"] != "-"]
    df_biogrid = df_biogrid[["Official Symbol Interactor A", "Official Symbol Interactor B"]]
    df_biogrid.columns = ["gene1", "gene2"]
    ppi_network = nx.from_pandas_edgelist(df_biogrid, source="gene1", target="gene2")

    return ppi_network


def load_combined_solutions(target_folder, dataset, directed=False, verbose=False):

    """
    Combines the TRIPS result for a dataset
    target_folder:  Should contain the *_pcst.txt files after running TRIPS for the dataset
    dataset: keyword or dataset name that should match those in the output folder
    """

    files = glob.glob(os.path.join(target_folder, "{}_module_[0-9]*_pcst.txt".format(dataset)))
    indices = [int(f.split("module_")[1].split("_pcst")[0]) for f in files]
    if verbose:
        print("======================================")
        print("files", files)
        print("indices", indices)
        print("======================================")

    if directed:
        G_big = nx.DiGraph()
        if indices:
            if len(set(indices)) == 1:
                file_undirected = os.path.join(target_folder, "{}_module_{}_pcst.txt".format(dataset, indices[0]))
                if os.path.exists(file_undirected):
                    df = pd.read_csv(file_undirected, sep="\t")
                    G_directed = nx.from_pandas_edgelist(df, "source", "target", create_using=nx.DiGraph)
                    G_big = nx.compose(G_big, G_directed)
            else:
                max_index = max(indices)
                for index in range((max_index+1)):
                    file_undirected = os.path.join(target_folder, "{}_module_{}_pcst.txt".format(dataset, index))
                    if os.path.exists(file_undirected):
                        df = pd.read_csv(file_undirected, sep="\t")
                        G_directed = nx.from_pandas_edgelist(df, "source", "target", create_using=nx.DiGraph)
                        G_big = nx.compose(G_big, G_directed)
        else:
            print("No output files found.")

    else:
        G_big = nx.Graph()
        if indices:
            if len(set(indices)) == 1:
                file_undirected = os.path.join(target_folder, "{}_module_{}_pcst.txt".format(dataset, indices[0]))
                if os.path.exists(file_undirected):
                    df = pd.read_csv(file_undirected, sep="\t")
                    G_undirected = nx.from_pandas_edgelist(df, "source", "target", create_using=nx.Graph)
                    G_big = nx.compose(G_big, G_undirected)
            else:
                max_index = max(indices)
                for index in range((max_index+1)):
                    file_undirected = os.path.join(target_folder, "{}_module_{}_pcst.txt".format(dataset, index))
                    if os.path.exists(file_undirected):
                        df = pd.read_csv(file_undirected, sep="\t")
                        G_undirected = nx.from_pandas_edgelist(df, "source", "target", create_using=nx.Graph)
                        G_big = nx.compose(G_big, G_undirected)

        else:
            print("No output files found.")

    return G_big



def load_domino_solns(file_domino_modules, mapping):

    indiv_modules = []
    dict_domino = {}
    mapping_rev = {v: k for k, v in mapping.items()}

    result_genes = []
    if os.path.exists(file_domino_modules):
        with open(file_domino_modules, 'r') as results:
            for line in results:
                result_genes.append(line.strip())

        indiv_modules = []
        for module in result_genes:
            module_genes = list(module[1:-1].split(", "))
            module_genes = [mapping_rev[x] for x in module_genes if x in mapping_rev]
            indiv_modules.append(module_genes)

    return indiv_modules