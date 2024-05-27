import argparse
import os

from trips_module.domino import DominoRunner
from trips_module.dapcstp import DAPCSTP
from trips_module.assign_scores import *
from trips_module.utils import *
# from trips_module.network_generators import  *


def validate_path(f):

    if not os.path.exists(f):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(f))
    return f


def parse_args():

    parser = argparse.ArgumentParser(description="Run TRIPS on gene expression data.",)

    parser.add_argument('--ppi', type=validate_path, required=True, help="Path to PPI network file. Should have two columns: source and target.")
    parser.add_argument('--grn', type=validate_path, required=True, help="Path to gene regulatory network file. Should have two columns: source and target.")
    parser.add_argument('--degs', type=validate_path, required=True, help="Path to the differential expression results. Should contain columns adjustedpvalue and logfc.")
    parser.add_argument('--path_domino', type=validate_path, required=True, help="Path to DOMINO.")
    parser.add_argument('--path_dapcstp', type=validate_path, required=True, help="Path to DAPCSTP.")
    parser.add_argument('--output_folder', type=validate_path, required=True, help="Specify the output folder.")
    parser.add_argument('--lfc_thresh', type=float, default=1.0, help="The logFoldChange threshold")
    parser.add_argument('--p', type=float, default=50, help="Set the pth percentile for calculating the edge cost.")
    parser.add_argument('--pval_thresh', type=float, default=0.05, help="Adjusted p-value threshold.")
    parser.add_argument('--prize_multiplier', type=float, default=1.0, help="Set the prize multipliers for DEGs or nodes in the gene module.")
    parser.add_argument('--mode', choices=["module","degs"], default="module", help="module will perform the TRIPS workflow while degs will perform the DEGs+PCST workflow.")
    parser.add_argument('--min_comm_size', type=int, default=0)

    args = parser.parse_args()
    return args


def main_trips(args):

    file_degs = args.degs
    file_grn = args.grn
    file_ppi = args.ppi
    path_to_dapcstp = args.path_to_dapcst
    path_to_domino = args.path_to_domino
    output_folder = args.output_folder

    lfc_thresh = args.lfc_thresh
    pval_thresh = args.pval_thresh
    prize_multiplier = args.prize_multiplier
    p = args.p
    mode = args.mode

    G_ppi = get_ppi_net(file_ppi)
    G_grn = get_grn(file_grn)

    run_trips_one(G_ppi, G_grn, file_degs, path_to_domino, path_to_dapcstp,
                  output_folder,
                  keyword="dataset1",
                  lfc_thresh=lfc_thresh, pval_thresh=pval_thresh,
                  edge_cost=p, prize_multiplier=prize_multiplier,
                  mapping=None,
                  min_comm_size=0,
                  directed_grn=False,
                  mode=mode,
                  custom_ppi=None)


def run_trips_one(ppi_network, G_grn, file_degs, path_to_domino, path_to_dapcst,
              main_output_folder,
              keyword="dataset1",
              lfc_thresh=1.0, pval_thresh=0.05,
              pp=50, prize_multiplier=1.0,
              mapping=None,
              min_comm_size=0,
              directed_grn=False,
              score_others="min",
              mode="module",
              domino_done=False,
              domino_modules=None,
              custom_hub_penalty=None, hybrid_score=None,
                         penalize_degs=False):

    print("Parameters:")
    # print("directed GRN: ", directed_grn)
    print("lfc threshold: ", lfc_thresh)
    print("pval threshold: ", pval_thresh)
    print("edge cost percentile: ", pp)

    if (domino_done == True) and (domino_modules is None):
        raise ValueError('Since domino_done is set to True, the location of DOMINO modules must be provided.')

    if not os.path.exists(main_output_folder):
        print("Creating main output folder...")
        os.makedirs(main_output_folder)

    final_output_folder = os.path.join(main_output_folder, "final_solutions_{}".format(keyword))
    if not os.path.exists(final_output_folder):
        os.mkdir(final_output_folder)

    # =============Get the PPI and GRN=============
    if directed_grn:
        assert(G_grn.is_directed()), "GRN should be directed."
    else:
        assert (G_grn.is_directed() == False), "GRN should be undirected."

    # =============Get the gene modules=============
    # degs = get_degs(file_degs, lfc_thresh=lfc_thresh, pval_thresh=pval_thresh)
    df_degs = pd.read_csv(file_degs, sep=",")
    df_degs = df_degs[df_degs["AdjPValue"] < pval_thresh]
    df_degs = df_degs[abs(df_degs["Log_FoldChange"]) > lfc_thresh]
    df_degs["score"] = -1 * (np.log10(df_degs["AdjPValue"]))
    dict_scores = dict(zip(df_degs["Gene_symbol"], df_degs["score"]))
    degs = list(dict_scores.keys())
    print("Found {} DEGs.".format(len(degs)))

    print("===========Running DOMINO===========")

    if domino_done:
        all_modules = load_domino_solns(file_domino_modules=domino_modules, mapping=mapping)
        print("Loaded {} active modules.".format(len(all_modules)))
    else:
        # Create DOMINO output folder
        domino_output_folder = os.path.join(main_output_folder, "domino_outputs")
        if not os.path.exists(domino_output_folder):
            os.mkdir(domino_output_folder)
        # Run DOMINO to get active modules on the PPI network
        dominor = DominoRunner(ppi_network, path_to_domino, domino_output_folder)
        all_modules = dominor.run_domino(degs, keyword=keyword, min_comm_size=min_comm_size)
        print("Found {} active modules.".format(len(all_modules)))

    # =============Find the DAPCST solutions=============
    print("===========Finding PCST solutions===========")
    dpcst = DAPCSTP(G_grn)

    # Calculate adaptive edge cost
    edge_cost = np.percentile(np.array(list(dict_scores.values())), pp)
    print("Edge cost: ", edge_cost)

    # Will run TRIPS
    if mode == "module":
        for index, domino_module in enumerate(all_modules):

            # dict_scores = calculate_gene_scores_retaindegs(file_degs, G_grn, lfc_thresh, pval_thresh, score_others=score_others,
            #                                            prize_multiplier=prize_multiplier, gene_module=domino_module)
            dict_scores = calculate_gene_scores(file_degs, G_grn, lfc_thresh, pval_thresh, score_others=score_others,
                                                prize_multiplier=prize_multiplier, gene_module=domino_module, 
                                                custom_hub_penalty=custom_hub_penalty, hybrid_score=hybrid_score, penalize_degs=penalize_degs)

            print("No. of genes in the domino module: ", len(domino_module))
            mod = "{}_module_{}".format(keyword, index)

            G_pcst = dpcst.run_dapcst_on_module(path_to_dapcst, degs, domino_module, dict_scores, edge_cost=edge_cost,
                                                 keyword=mod, final_output_folder=final_output_folder, directed=directed_grn)
            if G_pcst:
                # Write results to file
                file_out = os.path.join(final_output_folder, "{}_pcst.txt".format(mod))
                df = nx.to_pandas_edgelist(G_pcst)
                df.to_csv(file_out, sep="\t", index=False)
            else:
                print("No solution found for module ", index)

    # Will run the DEGs+PCST workflow
    elif mode == "degs":

        dict_scores = calculate_gene_scores_degsonly(file_degs, lfc_thresh, pval_thresh,
                                                       prize_multiplier=prize_multiplier)

        G_pcst = dpcst.run_dapcst_on_module(degs, path_to_domino, degs, dict_scores, edge_cost=edge_cost,
                                             keyword=mode, final_output_folder=final_output_folder)
        if G_pcst:
            # Write results to file
            file_out = os.path.join(final_output_folder, "{}_pcst.txt".format(keyword))
            df = nx.to_pandas_edgelist(G_pcst)
            df.to_csv(file_out, sep="\t", index=False)
        else:
            print("No solution found for DEGs.")

    print("Final solutions were written to folder: ", final_output_folder)

    # Load DAPCSTP solutions
    G_all = load_combined_solutions(final_output_folder, keyword, directed=False, verbose=False)

    # Save the DOMINO results
    file_output = os.path.join(final_output_folder, f"{keyword}_DOMINO_modules_one_per_line.txt")
    with open(file_output, "w") as f:
        for mod in all_modules:
            line = "\t".join(mod)
            f.write(line + "\n")

    # Save the final TRIPS output
    file_out = os.path.join(final_output_folder, f"{keyword}_trips_solution.txt")
    df_soln = nx.to_pandas_edgelist(G_all)
    df_soln.to_csv(file_out, sep="\t", index=False)
    
    return all_modules, G_all

if __name__ == "__main__":
    main_trips()