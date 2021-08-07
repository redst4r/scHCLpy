import networkx as nx
import obonet

"""
Code to deal with cell ontology celltypes.
In particular, it provides
- G_CL, a graph of ONLY cell-ontology celltypes, linked via the 'is_a', 'part_of' edges. (Note that CL has these "develops from" edges, which tend to shortcircuit)
- G_uberon, a graph of uberon anatomies
"""



graph = obonet.read_obo('http://purl.obolibrary.org/obo/cl/cl.obo')

# the uberon only subgraph
G_uberon = nx.subgraph(graph, [n for n in graph.nodes if n.startswith('UBERON')]).copy()
rm_edges = [(a, b, edgetype) for (a, b, edgetype) in G_uberon.edges if edgetype not in['is_a', 'part_of']]
G_uberon.remove_edges_from(rm_edges)

# the cell ontology subgraph
G_CL = nx.subgraph(graph, [n for n in graph.nodes if n.startswith('CL')]).copy()
rm_edges = [(a,b, edgetype) for (a,b,edgetype) in G_CL.edges if edgetype not in['is_a', 'part_of']]
G_CL.remove_edges_from(rm_edges)


def node2name(n, G):
    """
    turn a CL id into a celltype name
    """
    return G.nodes[n]['name']

def get_uberon_sub_anatomy(uberon_id):
    """
    any anatomy thats more specialized then the query term, ie 'pancreas' -> 'islet of Langerhans'
    """
    # subset the graph to uberon, and only is_a and part_of relations.
    # the more specialized cell types are ancestors
    uberon_ids = [_ for _ in nx.ancestors(G_uberon, uberon_id) if _.startswith('UBERON:')]
    return uberon_ids


def get_subcelltypes(CL_id):
    """
    get all subtypes of a cell type (i.e. more specialized cells)
    """
    return nx.ancestors(G_CL, CL_id)


def get_celltypes_of_anatomy(uberon_id):
    """
    get celltypes linked directly to the uberonterm and that celltypes speci
    """
   # CL terms directyly linked to the UBERON id
    candidates = [_ for _ in graph.predecessors(uberon_id) if _.startswith('CL:')]

    # any more specialized cells than the candidates also are part of the anatomy
    celltypes = []
    for c in candidates:
        celltypes.extend(get_subcelltypes(c))
    celltypes.extend(candidates)
    return celltypes

def get_celltypes_of_anatomy_and_subtypes(top_uberon):
    """
    get all the celltypes that are downstream of a given anatomy term
    """
    uberon_subtypes = get_uberon_sub_anatomy(top_uberon)
    uberon_subtypes.append(top_uberon)  # add the node itself in case it has any CL terms linked
    CLs = []
    for uid in uberon_subtypes:
        c = get_celltypes_of_anatomy(uid)
        CLs.extend(c)
    return list(set(CLs))
