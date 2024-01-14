import pandas as pd
import networkx as nx
import pickle
import os
import subprocess
import collections

from trips_module.utils import *


class DAPCSTP:

    def __init__(self, G_grn):
        self.G_grn = G_grn

    def pd_edges_from_stp_file(self, file_stp):

        """
        Reads an stp file and converts to a pandas dataframe
        Returns df with columns ["from","to"]
        """

        with open(file_stp, 'r') as f:

            edge_tuples = []
            lines = f.readlines()
            for line in lines:
                if line.startswith("E "):
                    edge = line.split()
                    edge_tuples.append((edge[1], edge[2]))

            return pd.DataFrame(edge_tuples, columns=["from", "to"])


    def write_stp_from_nx_graph(self, G_new, dict_scores2, domino_module_ind, edge_cost, file_out, keyword):

        with open(file_out, 'w') as f:

            lines1 = ["33D32945 STP File, STP Format Version 1.0\n", "SECTION Comments", keyword, "END", "\n"]
            lines2 = ["SECTION Graph", "Nodes {}".format(G_new.number_of_nodes()),
                      "Edges {}".format(G_new.number_of_edges())]
            lines = lines1 + lines2
            for line in lines:
                f.write(f"{line}\n")

            # Write the edges
            for edge in G_new.edges(data=True):
                f.write(f"E {edge[0]} {edge[1]} {edge_cost}\n")
            f.write("END\n\n")

            # Write the terminals
            f.write("SECTION Terminals\n")
            f.write("Terminals {}\n".format(len(dict_scores2)))
            for terminal, score in dict_scores2.items():
                if terminal in domino_module_ind:
                    f.write(f"TP {terminal} {score}\n")
            f.write("END\n\n")
            f.write("EOF")
            
            
    def write_stp_from_nx_graph_directed(self, G_new, dict_scores2, domino_module_ind, edge_cost, file_out, keyword):

        with open(file_out, 'w') as f:

            lines1 = ["33D32945 STP File, STP Format Version 1.0\n", "SECTION Comments", keyword, "END", "\n"]
            lines2 = ["SECTION Graph", "Nodes {}".format(G_new.number_of_nodes()),
                      "Arcs {}".format(G_new.number_of_edges())]
            lines = lines1 + lines2
            for line in lines:
                f.write(f"{line}\n")

            # Write the edges
            for edge in G_new.edges(data=True):
                f.write(f"A {edge[0]} {edge[1]} {edge_cost}\n")
            f.write("END\n\n")

            # Write the terminals
            f.write("SECTION Terminals\n")
            f.write("Terminals {}\n".format(len(dict_scores2)))
            for terminal, score in dict_scores2.items():
                if terminal in domino_module_ind:
                    f.write(f"TP {terminal} {score}\n")
            f.write("END\n\n")
            f.write("EOF")


    def generate_dapcst_instance(self, all_degs, dict_scores, gene_module=[]):

        """
        all_degs: list of gene names of all DEGs
        dict_scores: dict of {gene name: node prize}
        """

        print("No. of nodes in the GRN:", self.G_grn.number_of_nodes())
        dict_geneids = {gene_name: gene_name for gene_name in self.G_grn.nodes()}
        nx.set_node_attributes(self.G_grn, dict_geneids, "GeneID")
        G_new = nx.convert_node_labels_to_integers(self.G_grn, 1)
        dict_node_names = nx.get_node_attributes(G_new, "GeneID")  # integer_id:gene_name
        dict_rev = {v: k for k, v in dict_node_names.items()}  # gene_name:integer_id

        degs_in_g = set(all_degs).intersection(set(self.G_grn.nodes()))
        print("No. of DEGs that are in the GRN: ", len(degs_in_g))

        if gene_module:
            gene_module_in_g = set(gene_module).intersection(set(self.G_grn.nodes()))
            dict_scores2 = {dict_rev[gene]: dict_scores[gene] for gene in gene_module_in_g if gene in dict_scores}
            print("No. of genes in the gene module in the GRN: ", len(dict_scores2))
        else:
            dict_scores2 = {dict_rev[gene]: dict_scores[gene] for gene in degs_in_g if gene in dict_scores}
            print("No. of DEGs from gene module in the GRN: ", len(dict_scores2))

        return G_new, dict_node_names, dict_scores2


    def run_dapcst_on_module(self, path_to_dapcstp, degs, domino_module, dict_scores, edge_cost, keyword,
                             final_output_folder,
                             directed=False):

        # Generate the dapcstp instance and stp file
        G_new, dict_node_names, dict_scores2 = self.generate_dapcst_instance(degs, dict_scores,
                                                                             gene_module=domino_module)
        file_stp = os.path.join(final_output_folder, "{}.stp".format(keyword))
        domino_module_in = [k for k, v in dict_node_names.items() if v in domino_module]

        try:

            if directed == True:
                self.write_stp_from_nx_graph_directed(G_new, dict_scores2, domino_module_in,
                                                               edge_cost, file_stp, keyword)
            else:
                self.write_stp_from_nx_graph(G_new, dict_scores2, domino_module_in, edge_cost, file_stp, keyword)

            # Run the PCST solver
            file_out = os.path.join(final_output_folder, "solution_{}.stp".format(keyword))
            command = f'cd {path_to_dapcstp}; ./dapcstp {file_stp} --type pcstp -o {file_out}'
            subprocess.call(command, shell=True)

            if os.path.exists(file_out):

                # Load the PCST solution
                df = self.pd_edges_from_stp_file(file_out).astype(int)
                if directed:
                    G = nx.from_pandas_edgelist(df, source="from", target="to", create_using=nx.DiGraph)
                else:
                    G = nx.from_pandas_edgelist(df, source="from", target="to")

                G = nx.relabel_nodes(G, dict_node_names)
                return G

            else:
                return None

        except Exception as e:
            pass

