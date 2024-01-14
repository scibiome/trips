import networkx as nx
import numpy as np
import graph_tool.all as gt
import networkit as nk


from trips_module.utils import *


def convert_nx_to_nk_graph(G_enhanced, with_weight=True):

    # Convert networkx to networkit graph
    if with_weight:
        nxG = nk.nxadapter.nx2nk(G_enhanced, weightAttr="weight")
    else:
        nxG = nk.nxadapter.nx2nk(G_enhanced)

    nxGnodes = list(nxG.iterNodes())

    # Create mapping from node ids in nxG to node ids in nx
    nx_nodelist = list(G_enhanced.nodes())
    idmap = dict((id1, u) for (id1, u) in zip(nxGnodes, nx_nodelist))
    idmap_rev = dict({v: k for k, v in idmap.items()})

    return nxG, idmap, idmap_rev


def generate_RDPN(ggi_network, seed):

    """Generates random GGI network where nodes keep the exactly same degrees.
    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.
    Returns
    -------
    RDPN_network : nx.Graph
        Random GGI network where nodes keep their expected degrees.
    """

    d = nx.to_dict_of_lists(ggi_network)
    edges = [(i, j) for i in d for j in d[i]]
    GT = gt.Graph(directed=False)
    GT.add_vertex(sorted(ggi_network.nodes())[-1])
    GT.add_edge_list(edges)

    gt.random_rewire(GT,model = "constrained-configuration", n_iter = 100, edge_sweep = True)

    edges_new = list(GT.get_edges())
    edges_new = [tuple(x) for x in edges_new]
    rewired_network = nx.Graph()
    rewired_network.add_nodes_from(ggi_network.nodes())
    rewired_network.add_edges_from(edges_new)
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    nx.set_node_attributes(rewired_network, gene_ids, 'GeneID')

    return(rewired_network)


def generate_rewired_network(ggi_network, seed):

    """Generates random GGI network where nodes keep their expected degrees.
    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.
    Returns
    -------
    rewired_network : nx.Graph
        Random GGI network where nodes keep their expected degrees.
    """

    degree_view = ggi_network.degree()
    degree_sequence = [degree_view[node] for node in ggi_network.nodes()]
    rewired_network = nx.expected_degree_graph(degree_sequence, seed=seed, selfloops=False)
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    nx.set_node_attributes(rewired_network, gene_ids, 'GeneID')

    return rewired_network


def generate_shuffled_network(ggi_network, seed):

    """Generates random GGI network by shuffling the node labels (gene IDs).
    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.
    Returns
    -------
    shuffled_network : nx.Graph
        Random GGI network where the nodes labels are shuffled.
    """

    shuffled_network = nx.Graph(ggi_network)
    shuffled_gene_ids = list(nx.get_node_attributes(shuffled_network, 'GeneID').values())
    np.random.seed(seed=seed)
    np.random.shuffle(shuffled_gene_ids)
    shuffled_gene_ids = {node: shuffled_gene_ids[node] for node in shuffled_network.nodes()}
    nx.set_node_attributes(shuffled_network, shuffled_gene_ids, 'GeneID')

    return shuffled_network


def generate_scale_free_network(ggi_network, seed):

    """Generate random scale free network with as many nodes and approximately as many edges as the original network.
    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.
    Returns
    -------
    scale_free_network : nx.Graph
        Random scale free network with as many nodes and approximately as many edges as the original network.
    """

    num_nodes = ggi_network.number_of_nodes()
    num_edges = ggi_network.number_of_edges()
    m = np.max([1, round(num_nodes / 2.0 - np.sqrt((num_nodes * num_nodes) / 4.0 - num_edges))])
    scale_free_network = nx.barabasi_albert_graph(num_nodes, m, seed=seed)
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    nx.set_node_attributes(scale_free_network, gene_ids, 'GeneID')

    return scale_free_network


def generate_uniform_network(ggi_network, seed):

    """Generates random uniform network with as many nodes and edges as the original network.
    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.
    Returns
    -------
    uniform_network : nx.Graph
        Random uniform network with as many nodes and edges as the original network.
    """

    num_nodes = ggi_network.number_of_nodes()
    num_edges = ggi_network.number_of_edges()
    uniform_network = nx.gnm_random_graph(num_nodes, num_edges, seed=seed)
    gene_ids = nx.get_node_attributes(ggi_network, 'GeneID')
    nx.set_node_attributes(uniform_network, gene_ids, 'GeneID')
    return uniform_network


def generate_network2(ggi_network, shuffling_method, seed):

    """Generates random network based on the original network.
    Parameters
    ----------
    ggi_network : nx.Graph
        Original GGI network.
    seed : convertible to np.uint32 or None
        Seed used by the generator.
    network_generator_selector : utils.NetworkGeneratorSelector
        Specifies which generator should be used.
    Returns
    -------
    random_network : nx.Graph
        Random network that was generated based on the original GGI network.
    """

    if shuffling_method == "EXPECTED_DEGREE":
        return generate_rewired_network(ggi_network, seed)
    elif shuffling_method == "SHUFFLED":
        return generate_shuffled_network(ggi_network, seed)
    elif shuffling_method == "SCALE_FREE":
        return generate_scale_free_network(ggi_network, seed)
    elif shuffling_method == "UNIFORM":
        return generate_uniform_network(ggi_network, seed)


def degree_preserving_randomize(G_nx, seed=42, verbose=False):

    """
    Perform degree-preserving randomization
    of a networkx graph using networkit functions
    Input: networkx graph
    Output: networkx graph
    Note: does not consider any weight attributes
    """

    # Set seed
    nk.setSeed(seed, True)

    # Convert back to networkx graph and relabel nodes
    G_tf, idmap, idmap_rev = convert_nx_to_nk_graph(G_nx, with_weight=False)
    if verbose:
        print()
        print("No. of nodes of shuffled network: ", G_tf.numberOfNodes())
        print("No. of edges of shuffled network: ", G_tf.numberOfEdges())

    # Run degree-preserving randomization
    dps = nk.randomization.DegreePreservingShuffle(G_tf)

    # Run algorithm
    dps.run()

    G_tf_shuffled = dps.getGraph()
    # Verify
    for u in range(G_tf.upperNodeIdBound()):
        assert (G_tf.degree(u) == G_tf_shuffled.degree(u))

    # Relabel nodes
    g_shuffled = nk.nxadapter.nk2nx(G_tf_shuffled)
    if verbose:
        print("No. of nodes: ", g_shuffled.number_of_nodes())
        print("No. of edges: ", g_shuffled.number_of_edges())
    g_new = nx.relabel_nodes(g_shuffled, idmap)

    return g_new


def shuffle_nx_graph(G, shuffling_method="REWIRED", seed=42):

    if shuffling_method == "REWIRED":
        G_shuffled = degree_preserving_randomize(G, seed=seed)

    else:

        # add GeneID attribute
        dict_geneids = {gene_name: gene_name for gene_name in G.nodes()}
        nx.set_node_attributes(G, dict_geneids, "GeneID")
        G_new = nx.convert_node_labels_to_integers(G)
        # Actual shuffling happens here
        G_shuffled = generate_network2(G_new, shuffling_method, seed)
        dict_rename = nx.get_node_attributes(G_shuffled, "GeneID")
        G_shuffled = nx.relabel_nodes(G_shuffled, dict_rename, "GeneID")

    return G_shuffled

