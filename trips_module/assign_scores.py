
import pandas as pd
import numpy as np
import networkx as nx


def calculate_gene_scores_retaindegs(file_degs, G_grn, lfc_thresh=1.0, pval_thresh=0.05,
                          prize_multiplier=1, mu=0.05, gene_module=[],
                          score_others="mean", adjusted=False, custom_hub_penalty=None):
    """
    file_degs: txt file the differential expression analysis results
                should have columns ["Gene","log2fc","AdjPval"]
    G_grn: (undirected) networkx graph
    prize_multiplier: multiplier for the abs(logFC) for DEGs in the module
    gene_module: list of genes in a module (e.g. DOMINO module)
                genes in the module but not DEGs will be assigned
                the mean scores of all DEGs as prize
    custom_hub_penalty: dictionary of gene:degree (could be from literature bias, etc)

    Returns
    dict_scores: dictionary where key is gene and value is the calculated prize
                including the scores of genes in the module (in addition to DEGs)

    The prizes are proportional to the abs logFC determined from the file_degs file
    """

    assert len(gene_module) > 1, "The gene module should not be empty!"
    print("len(gene_module)", len(gene_module))
    print("=============================CALCULATING GENE SCORES==============================")

    df = pd.read_csv(file_degs, sep="\t")
    df_sub = df[df["AdjPValue"] < pval_thresh]
    if df_sub.empty:
        raise Exception("No genes satisfiying adjusted p-value cutoff")
    df_sub = df_sub[abs(df_sub["Log_FoldChange"]) > lfc_thresh]
    if df_sub.empty:
        raise Exception("No genes satisfiying logFC cutoff")

    df_sub["score"] = -np.log(df_sub["AdjPValue"]) * prize_multiplier
    dict_scores = dict(zip(df_sub["Gene_symbol"], df_sub["score"]))
    all_degs = list(dict_scores.keys())
    dict_scores2 = {k: v for k, v in dict_scores.items() if k in gene_module}  # DEGs in module

    if gene_module:
        if score_others == "mean":
            if len(dict_scores2) > 0:
                mean_score = np.mean(list(dict_scores2.values()))
                for gene in gene_module:
                    if gene not in dict_scores2:
                        dict_scores2[gene] = mean_score * prize_multiplier
        elif score_others == "min":
            if len(dict_scores2) > 0:
                min_score = min(list(dict_scores2.values()))
                for gene in gene_module:
                    if gene not in dict_scores2:
                        dict_scores2[gene] = min_score * prize_multiplier

    # Apply hub penalty
    if custom_hub_penalty:
        hub_penalty = dict(custom_hub_penalty)
    else:
        hub_penalty = dict(G_grn.degree())

    # Exclude original DEGs from penalty
    dict_scores_a = {k: (v - (mu * hub_penalty[k])) for k, v in dict_scores2.items() if ((k not in all_degs) and (k in hub_penalty))}
    dict_scores_b = {k:v for k, v in dict_scores2.items() if k in all_degs}
    dict_scores2 = {**dict_scores_a, **dict_scores_b}
    dict_scores2 = dict((k, 0) if v < 0 else (k, v) for k, v in dict_scores2.items())

    return dict_scores2


def calculate_gene_scores_degsonly(file_degs, G_grn, lfc_thresh=0.58, pval_thresh=0.05, mode="both",
                          prize_multiplier=1, mode_pval=True, adjusted=False):
    """
    file_degs: txt file (obtained from the GREIN interface) containing
                the differential expression analysis results
                should have columns ["Log_FoldChange","AdjPValue"]
    prize_multiplier: multiplier for the abs(logFC) for DEGs in the module
    gene_module: list of genes in a module (e.g. DOMINO module)
                genes in the module but not DEGs will be assigned
                the mean scores of all DEGs as prize

    Returns
    dict_scores: dictionary where key is gene and value is the calculated prize
                including the scores of genes in the module (in addition to DEGs)

    The prizes are proportional to the abs logFC determined from the file_degs file
    """

    df = pd.read_csv(file_degs, sep="\t")

    if adjusted:
        df_sub = df[df["AdjPValue"] < pval_thresh]
        if df_sub.empty:
            raise Exception("No genes satisfiying adjusted p-value cutoff")
    else:
        df_sub = df[df["PValue"] < pval_thresh]

        if df_sub.empty:
            raise Exception("No genes satisfiying p-value cutoff")

    if mode == "up":
        df_sub = df_sub[df_sub["Log_FoldChange"] > lfc_thresh]
    elif mode == "down":
        df_sub = df_sub[df_sub["Log_FoldChange"] < -lfc_thresh]
    elif mode == "both":
        df_sub = df_sub[abs(df_sub["Log_FoldChange"]) > lfc_thresh]
    if df_sub.empty:
        raise Exception("No genes satisfiying logFC cutoff")

    if mode_pval:
        df_sub["score"] = -np.log(df_sub["AdjPValue"]) * prize_multiplier
    else:
        df_sub["score"] = abs(df_sub["Log_FoldChange"]) * prize_multiplier

    dict_scores = dict(zip(df_sub["Gene_symbol"], df_sub["score"]))

    return dict_scores
