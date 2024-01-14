
import os
import glob
import pandas as pd
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from trips_module.utils import *


def get_lcc(G):

    """
    Extract the largest connected component (LCC) from a graph
    Input
    G: (undirected) networkx graph
    """

    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    G0 = G.subgraph(Gcc[0])

    return G0

def get_unique_nodes(df):

    """
    df: a pandas edgelist dataframe with two columns
    """

    df.columns = ["from", "to"]
    list1 = df["from"]
    list2 = df["to"]
    all_nodes = list(list1) + list(list2)
    all_nodes = list(set(all_nodes))

    return all_nodes


def jaccard_similarity(list1, list2):
    
    if len(list1) == 0 or len(list2) == 0:
        return 0.0
    else:
        intersection = len(list(set(list1).intersection(list2)))
        print("intersection", intersection)
        union = (len(set(list1)) + len(set(list2))) - intersection
        return float(intersection) / union


def calculate_precision(list1, gold_standard):

    if len(list1) == 0:
        return 0.0
    else:
        com = list(set(list1).intersection(set(gold_standard)))
        return float(len(com) / len(list1))


def extract_cacts_results(file, n_top_reg=None):

    if os.path.exists(file):
        df_res = pd.read_csv(file, sep=" ")
        if n_top_reg:
            df_res = df_res.head(n_top_reg)
        all_covers = list(df_res["Name"])
        return all_covers
    else:
        return []
        
def extract_aracne_results(file, all_tfs=[], n_top_edges=10000):
    
    df = pd.read_csv(file, sep="\t")
    df = df.head(n_top_edges)
    G = nx.from_pandas_edgelist(df, "from","to")
    G_tfs = set(G.nodes()).intersection(set(all_tfs))
    return list(G_tfs)

def extract_mogamun_results(file):
    
    df = pd.read_csv(file, sep=" ")
    G = nx.from_pandas_edgelist(df, source="V1", target="V2")
    return list(G.nodes())
    

def extract_pcsf_results(file_pcsf):

    """
    Extracts the node list of the PCSF output
    """

    df = pd.read_csv(file_pcsf, sep="\t")
    df.columns = ["from", "to"]
    all_nodes = get_unique_nodes(df)

    return all_nodes


def extract_regenrich_results(folder_regenrich, dataset, sep="\t", n_top_reg=None):

    """
    Extracts the node list of the RegEnrich outputs
    """

    file = os.path.join(folder_regenrich, "{}_rankscore_GSEA.txt".format(dataset))
    if os.path.exists(file):
        df = pd.read_csv(file, sep=sep)
        if n_top_reg:
            df = df.head(n_top_reg)
        all_covers = list(df.reg)
        return all_covers
    else:
        return []


def extract_trips_results(folder, dataset, shuffling_method=None, directed=False, verbose=False):

    """
    Extracts the node list of the RegEnrich output
    """
    soln = []
    target_folder = os.path.join(folder, dataset)
    if shuffling_method:
        soln = load_combined_solutions_shuffled(target_folder, dataset, directed=directed, verbose=verbose,
                                                shuffling_method=shuffling_method)
    else:
        soln = load_combined_solutions(target_folder, dataset, directed=directed, verbose=verbose)

    return soln


def extract_diamond_pcst_results(folder_diam_pcst, dataset, directed=False, shuffling_method=None):

    soln = []

    if shuffling_method:

        file_soln = os.path.join(folder_diam_pcst, dataset, "{}_{}_pcst.txt".format(dataset, shuffling_method))
        if os.path.exists(file_soln):
            df = pd.read_csv(file_soln, sep="\t")
            if directed:
                G = nx.from_pandas_edgelist(df, "source" ,"target", create_using=nx.DiGraph)
            else:
                G = nx.from_pandas_edgelist(df, "source", "target", create_using=nx.Graph)
            soln = list(G.nodes())
        else:
            print("No. DIAMOnD+PCST files found.")

    else:
        file_soln = os.path.join(folder_diam_pcst, dataset, "{}_pcst.txt".format(dataset))
        if os.path.exists(file_soln):
            df = pd.read_csv(file_soln, sep="\t")
            if directed:
                G = nx.from_pandas_edgelist(df, "source" ,"target", create_using=nx.DiGraph)
            else:
                G = nx.from_pandas_edgelist(df, "source", "target", create_using=nx.Graph)
            soln = list(G.nodes())
        else:
            print("No. DIAMOnD+PCST files found.")

    return soln


def load_combined_solutions_shuffled(target_folder, dataset, directed=False, verbose=False, shuffling_method=None):

    """
    Combines the TRIPS/dapcst result for a dataset
    target_folder: contains the _pcst.txt files after running TRIPS for the dataset
    dataset: keyword or dataset name
    """

    files = glob.glob(os.path.join(target_folder, "*_module_[0-9]*_pcst.txt"))
    indices = [int(f.split("module_")[1].split("_pcst")[0]) for f in files]
    if verbose:
        print("======================================")
        print("files", files)
        print("indices", indices)

    if directed:
        G_big = nx.DiGraph()
        if indices:
            if len(set(indices)) == 1:
                file_undirected = os.path.join(target_folder, "{}_{}_module_{}_pcst.txt".format(dataset, shuffling_method, indices[0]))
                if os.path.exists(file_undirected):
                    df = pd.read_csv(file_undirected, sep="\t")
                    G_directed = nx.from_pandas_edgelist(df, "source", "target", create_using=nx.DiGraph)
                    G_big = nx.compose(G_big, G_directed)
            else:
                max_index = max(indices)
                for index in range((max_index+1)):
                    file_undirected = os.path.join(target_folder, "{}_{}_module_{}_pcst.txt".format(dataset, shuffling_method, index))
                    if os.path.exists(file_undirected):
                        df = pd.read_csv(file_undirected, sep="\t")
                        G_directed = nx.from_pandas_edgelist(df, "source", "target", create_using=nx.DiGraph)
                        G_big = nx.compose(G_big, G_directed)
        else:
            print("No. TRIPS output files found for workflow.")
    else:
        G_big = nx.Graph()
        if indices:
            if len(set(indices)) == 1:
                file_undirected = os.path.join(target_folder, "{}_{}_module_{}_pcst.txt".format(dataset, shuffling_method, indices[0]))
                if os.path.exists(file_undirected):
                    df = pd.read_csv(file_undirected, sep="\t")
                    G_undirected = nx.from_pandas_edgelist(df, "source", "target", create_using=nx.Graph)
                    G_big = nx.compose(G_big, G_undirected)
            else:
                max_index = max(indices)
                for index in range((max_index + 1)):
                    file_undirected = os.path.join(target_folder, "{}_{}_module_{}_pcst.txt".format(dataset, shuffling_method, index))
                    if os.path.exists(file_undirected):
                        df = pd.read_csv(file_undirected, sep="\t")
                        G_undirected = nx.from_pandas_edgelist(df, "source", "target", create_using=nx.Graph)
                        G_big = nx.compose(G_big, G_undirected)
        else:
            print("No. TRIPS output files found for workflow.")

    return G_big


def extract_degs_pcst_results(folder_degs_pcst, dataset, shuffling_method=None):

    soln = []
    if shuffling_method is not None:
        file_soln = os.path.join(folder_degs_pcst, dataset, "{}_{}_pcst.txt".format(dataset, shuffling_method))
    else:
        file_soln = os.path.join(folder_degs_pcst, dataset, "{}_pcst.txt".format(dataset))

    if os.path.exists(file_soln):
        df = pd.read_csv(file_soln, sep="\t")
        G = nx.from_pandas_edgelist(df, "source" ,"target", create_using=nx.DiGraph)
        soln = list(set(G.nodes()))

    return soln


def network_perturbation(folder, target_workflow, all_datasets, directed=False, title="", all_tfs=[], dict_dataset_map=None,
                        folder_degs=None, filename=None,
                        N=20, gold_standard=None, n_mech=None,
                        mode=None, pval_thresh=0.05, lfc_thresh = 1, verbose=False):

    """
    Generic method for evaluating robustness of a network-based method using perturbed backbone networks
    gold_standard: dictionary of disease DOIDs:[list of disease-associated genes]
    """

    if directed:
        all_shuffling_methods = ["REWIRED", "EDGE_SWITCHING", "EXPECTED_DEGREE", "UNIFORM"]
        print("Using directed networks.")
    else:
        all_shuffling_methods = ["REWIRED", "EXPECTED_DEGREE", "SCALE_FREE", "SHUFFLED"]
        print("Using undirected networks.")

    workflows = [target_workflow] + all_shuffling_methods

    dict_precisions = {k: [] for k in workflows}
    dict_recall = {k: [] for k in workflows}
    dict_soln_size = {k: [] for k in workflows}

    for dataset in all_datasets:

        if verbose:
                print("===============================")
                print("DATASET:", dataset)

        doid1 = dict_dataset_map[dataset]["DOID"]
        if doid1 == "80208":
            doid = "0080208"
        else:
            doid = doid1
        disease_genes = gold_standard[doid]
        disease_tfs = list(set(disease_genes).intersection(set(all_tfs)))

        for workflow in workflows:

            soln = []

            if target_workflow == "trips":
                if workflow == target_workflow:
                    soln = extract_trips_results(folder, dataset)
                else:
                    soln = extract_trips_results(folder, dataset, shuffling_method=workflow)

            elif target_workflow == "diamond_pcst":
                if workflow == target_workflow:
                    soln = extract_diamond_pcst_results(folder, dataset)
                else:
                    soln = extract_diamond_pcst_results(folder, dataset, shuffling_method=workflow)

            if target_workflow == "degs_pcst":
                if workflow == target_workflow:
                    soln = extract_degs_pcst_results(folder, dataset)
                else:
                    soln = extract_degs_pcst_results(folder, dataset, shuffling_method=workflow)

            elif target_workflow == "regenrich":

                if workflow == target_workflow:
                    folder_regenrich = os.path.join(folder, "REGENRICH")
                    soln = extract_regenrich_results(folder_regenrich, dataset, n_top_reg=N)
                else:
                    folder_regenrich = os.path.join(folder, workflow)
                    soln = extract_regenrich_results(folder_regenrich, dataset, n_top_reg=N)

            elif target_workflow == "pcsf":

                if workflow == target_workflow:
                    if n_mech:
                        folder_pcsf = os.path.join(folder, "PCSF", "{}_edges_{}.txt".format(dataset, n_mech))
                    else:
                        folder_pcsf = os.path.join(folder, "PCSF", "{}_edges.txt".format(dataset))
                    soln = extract_pcsf_results(folder_pcsf)
                else:
                    try:
                        if n_mech:
                            folder_pcsf = os.path.join(folder, workflow, "{}_edges_{}.txt".format(dataset, n_mech))
                        else:
                            folder_pcsf = os.path.join(folder, workflow, "{}_edges.txt".format(dataset))
                        soln = extract_pcsf_results(folder_pcsf)
                    except Exception as e:
                        pass
                        
            elif target_workflow == "mogamun":
                
                if workflow == target_workflow:
                    folder_sub = glob.glob(os.path.join(folder, "REAL", dataset, "Experiment_*".format(dataset)))
                    if len(folder_sub) == 1:
                        file = os.path.join(folder_sub[0], "A_ALL_FILTERED_INTERACTIONS_CYTOSCAPE.csv")
                        if os.path.exists(file):
                            soln = extract_mogamun_results(file)
                else:
                    folder_sub = glob.glob(os.path.join(folder, workflow, dataset, "Experiment_*".format(dataset)))
                    if len(folder_sub) == 1:
                        file = os.path.join(folder_sub[0], "A_ALL_FILTERED_INTERACTIONS_CYTOSCAPE.csv")
                        if os.path.exists(file):
                            soln = extract_mogamun_results(file)

            # We are interested only in the TFs
            all_covers = list(set(soln).intersection(set(all_tfs)))

            # ============UPDATE THE DICTIONARIES===============
            if all_covers:
                dict_precisions[workflow].append(calculate_precision(all_covers, disease_tfs))
            else:
                dict_precisions[workflow].append(float(0))

            common_genes = set(disease_tfs).intersection(set(all_covers))
            dict_recall[workflow].append(len(common_genes))
            dict_soln_size[workflow].append(len(all_covers))

        disease_tfs = [] # reset list
  
    return dict_precisions, dict_recall, dict_soln_size