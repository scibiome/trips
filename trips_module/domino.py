
import mygene
import pandas as pd
import networkx as nx
import subprocess


class DominoRunner:


    def __init__(self, ppi_network, path_to_domino, output_folder, mapping=None):

        """
        ppi_network:  nx graph of the PPI network to be used
        path_to_domino: the path to DOMINO
        output_folder: where to save the DOMINO results
        """

        self.path_to_domino = path_to_domino
        self.output_folder = output_folder
        
        if mapping:
            self.prepare_gene_mapping(ppi_network)
        else:
            self.mapping = mapping

        path_to_network = self.prepare_ppi_net_for_domino(ppi_network)
        self.path_to_network = path_to_network


    def prepare_gene_mapping(self, ppi_network):
    
        """
        Prepares the mapping from ensembl.gene ids (ENSGXXXX) to gene symbol using mygene
        """
    
        mg = mygene.MyGeneInfo()
        out = mg.querymany(list(ppi_network.nodes()), scopes="symbol", fields='ensembl.gene', species='human', verbose=False)
        mapping = dict()
        for line in out:
            try:
                res = line["ensembl"]
                if len(res) == 1:

                    mapping[line["query"]] = res["gene"]
                else:
                    mapping[line["query"]] = res[0]["gene"]

            except KeyError:
                mapping[line["query"]] = "pass"
                
        self.mapping = mapping


    def prepare_ppi_net_for_domino(self, ppi_network):

        # 1. Write GGI network in format required by your methodmapping[line["query"]] = res["gene"]
        path_to_network = f'{self.output_folder}/domino_ggi.cif'
        # gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')

        with open(path_to_network, 'w') as edge_list_file:
            edge_list_file.write(f'node_1\t combined_score \tnode_2\n')
            for u, v in ppi_network.edges():
                try:
                    edge_list_file.write(f'{self.mapping[u]}\t ppi \t{self.mapping[v]}\n')
                except:
                    pass
        path_to_network = f'{self.output_folder}/domino_ggi.cif'

        return path_to_network


    def run_domino(self, seed_genes, keyword="", min_comm_size=0):

        """Runs DOMINO.
        Parameters
        ----------
        seed_genes : list of str
            Seed genes (entries are gene IDs).
        min_comm_size : int
            Size filter used for finding the active modules
        Returns
        -------
        indiv_modules : list of list of str
            list of modules output by DOMINO.
        """

        assert (type(keyword) == str and len(keyword) > 0)

        outfile = f'{self.output_folder}/{keyword}_domino_slices.txt'
        command = f'cd {self.path_to_domino}; slicer --network_file {self.path_to_network} --output_file {outfile}'
        subprocess.call(command, shell = True)

        # 2. Write expression data and phenotype in the format required by your method
        seeds_ens = [self.mapping[x] for x in seed_genes if x in self.mapping]
        path_to_seeds = f'{self.output_folder}/{keyword}_domino_seeds.txt'
        path_to_seeds = path_to_seeds.replace(" ", "_")
        with open(path_to_seeds, 'w') as seeds_file:
            for g in seeds_ens:
                seeds_file.write(f'{g}\n')

        path_to_seeds = f'{self.output_folder}/{keyword}_domino_seeds.txt'
        path_to_seeds = path_to_seeds.replace(" ", "_")
        # print("path_to_seeds", path_to_seeds)
        # print("self.output_folder", self.output_folder)

        # 3. Insert the command to run your method, direct the output to path_to_output
        # path_to_output = f'../../temp/'
        command = f'cd {self.path_to_domino}; domino --active_genes_files {path_to_seeds} --network_file {self.path_to_network} --slices_file {outfile} --output_folder {self.output_folder} --visualization false'
        subprocess.call(command, shell = True)

        # 4. Process results such that they are formatted as a list of strings (entez IDs)
        # path_to_seeds = f'{output_folder}/{prefix}_domino_seeds.txt'

        result_genes = []
        with open(path_to_seeds[:-4] +"/modules.out", 'r') as results:
            # with open(os.path.join(output_folder, domino_test_domino_seeds, "modules.out"), 'r') as results:
            for line in results:
                result_genes.append(line.strip())

        mapping_rev = {v: k for k ,v in self.mapping.items()}
        indiv_modules = []
        for module in result_genes:
            module_genes = list(module[1:-1].split(", "))
            module_genes = [mapping_rev[x] for x in module_genes]
            indiv_modules.append(module_genes)

        # res =  flatten([x[1:-1].split(", ") for x in result_genes])
        # res = list(set(res))
        indiv_modules = [mod for mod in indiv_modules if len(mod) >= min_comm_size]

        # Delete temporary data.
        # subprocess.call(f'rm -r {output_folder}/{path_to_seeds[:-4]}', shell=True)
        # subprocess.call(f'rm ../temp/{prefix}_domino_*', shell=True)
        return indiv_modules