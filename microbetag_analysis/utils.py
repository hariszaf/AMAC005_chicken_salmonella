from variables import *
from stats import *
from visualization import *

from dataclasses import dataclass

# To parse network
#    pip install ndex2
from ndex2.cx2 import RawCX2NetworkFactory


import pandas as pd


def load_cx2(filepath):
    """
    Load a .cx2 network with ndex2
    """
    factory = RawCX2NetworkFactory()
    net = factory.get_cx2network(filepath)
    return net


def parse_microbetag_edges(cx2, env_set, metabolites_set):
    """
    Gets a ndex2 object as retruned from loading a .cx2 microbetag annotated network
    and returns

    Arguments:
    - cx2 (ndex2.cx2.CX2Network):
    - env_set (Set):
    - metabolites_set (Set):

    Returns:
    - edge_id_2_number_of_complements (Dict)
    - edge_types (Dict):
    - taxa_pairs_with_positive_weight (List):
    - edge_id_2_cross_feeding_compounds (Dict):
    """

    edge_id_sign = {}
    edge_id_2_number_of_complements = {}
    edge_id_2_cross_feeding_compounds = {}
    edge_types = {
        "taxon_to_taxon"  : {"edges": [], "counts": 0},  # done
        "taxon_to_metabo" : {"edges": [], "counts": 0},
        "taxon_to_env"    : {"edges": [], "counts": 0},
        "metabo_to_metabo": {"edges": [], "counts": 0},  # done
        "env_to_env"      : {"edges": [], "counts": 0},  # done
        "metabo_to_env"   : {"edges": [], "counts": 0},
    }

    source_target_edge_ids = {}
    taxa_pairs_with_positive_weight = []

    for _, edge in cx2.get_edges().items():

        edge_id = edge["id"]

        if edge["v"]["interaction type"] == COMPL_INT_TYPE:

            for k in edge["v"]:
                # NOTE For PATHWAY COMPLEMENTS
                # if "compl::" in k:
                #     edge_id_2_number_of_complements[edge_id] = len(edge["v"][k])
                uniq_cross_feeding_compounds = set()
                if k.startswith("seedCompl::"):
                    edge_id_2_number_of_complements[edge_id] = len(edge["v"][k])
                    for compl in edge["v"][k]:
                        for i in compl.split("^")[2].split(";"):
                            uniq_cross_feeding_compounds.add(i)
                    edge_id_2_cross_feeding_compounds[edge_id] = (
                        uniq_cross_feeding_compounds
                    )
        else:

            if "microbetag::weight" in edge["v"]:

                if edge["v"]["microbetag::weight"] > 0:
                    edge_id_sign[edge_id] = 1
                    taxa_pairs_with_positive_weight.append((edge["s"], edge["t"]))
                else:
                    edge_id_sign[edge_id] = -1

        source_target_edge_ids[(edge["s"], edge["t"])] = edge_id

        source, target = edge["s"], edge["t"]

        node_source = cx2.get_nodes()[source]
        source_name = node_source["v"]["name"]

        node_target = cx2.get_nodes()[target]
        target_name = node_target["v"]["name"]

        if (
            source_name not in env_set | metabolites_set
            and target_name not in env_set | metabolites_set
        ):
            edge_types["taxon_to_taxon"]["counts"] += 1
            edge_types["taxon_to_taxon"]["edges"].append(edge_id)

        elif source_name in env_set and target_name in env_set:
            edge_types["env_to_env"]["counts"] += 1
            edge_types["env_to_env"]["edges"].append(edge_id)

        elif source_name in metabolites_set and target_name in metabolites_set:
            edge_types["metabo_to_metabo"]["counts"] += 1
            edge_types["metabo_to_metabo"]["edges"].append(edge_id)

        elif {
            source_name,
            target_name,
        } <= env_set | metabolites_set and source_name in env_set ^ metabolites_set:
            edge_types["metabo_to_env"]["counts"] += 1
            edge_types["metabo_to_env"]["edges"].append(edge_id)

        elif {source_name, target_name} & env_set and not {
            source_name,
            target_name,
        } & metabolites_set:
            edge_types["taxon_to_env"]["counts"] += 1
            edge_types["taxon_to_env"]["edges"].append(edge_id)

        elif {source_name, target_name} & metabolites_set and not {
            source_name,
            target_name,
        } & env_set:
            edge_types["taxon_to_metabo"]["counts"] += 1
            edge_types["taxon_to_metabo"]["edges"].append(edge_id)

        else:
            print("PROBLEM:", source_name, target_name)

    return (
        edge_id_2_number_of_complements,
        edge_types,
        taxa_pairs_with_positive_weight,
        edge_id_2_cross_feeding_compounds,
    )


def get_compls_and_compounds_in_positive_associated_taxa(
    cx2, positive_pairs, edge_id_compl_number, compounds
):

    number_of_complements_in_positive_associated_taxa = 0
    unique_compounds = set()
    for _, edge in cx2.get_edges().items():

        s_t = (edge["s"], edge["t"])
        t_s = (edge["t"], edge["s"])
        edge_id = edge["id"]

        if s_t in positive_pairs or t_s in positive_pairs:
            if edge["v"]["interaction type"] == COMPL_INT_TYPE:
                for k in edge["v"]:
                    if k.startswith("seedCompl::"):
                        number_of_complements_in_positive_associated_taxa += (
                            edge_id_compl_number[edge_id]
                        )
                        for i in compounds[edge_id]:
                            unique_compounds.add(i)

    return number_of_complements_in_positive_associated_taxa, unique_compounds


@dataclass
class MggParser:
    """A `dataclass` for the parsed microbetag annotated network once loaded with `ndex2`; instances of this class will store outputs of the `parse_microbetag_egdes()`."""

    import ndex2

    num_total_compls           : dict
    edge_types                 : dict
    pos_pairs_node_ids         : list  # of tuples
    edge_id_to_unique_compounds: dict
    cx2                        : ndex2.cx2.CX2Network

    def __init__(self, cx2, env_set, metabolites_set) -> None:

        (
            self.num_total_compls,
            self.edge_types,
            self.pos_pairs_node_ids,
            self.edge_id_to_unique_compounds,
        ) = parse_microbetag_edges(cx2, env_set, metabolites_set)
        self.cx2 = cx2


def remove_cols_with_null(df):

    # Store original columns
    original_columns = df.columns

    # Drop columns with all NaNs
    df_clean = df.dropna(axis=1, how='any')

    # Find which columns were dropped
    dropped_columns = original_columns.difference(df_clean.columns)

    return df_clean, dropped_columns


def check_cooccurrence(s, t, cx2):
    """
    IN USE!
    """
    all_edges = cx2.get_edges()

    # Check for either direction (source -> target or target -> source)
    return any(
        ((edge["s"] == s and edge["t"] == t) or (edge["s"] == t and edge["t"] == s))
        and edge["v"]["interaction type"] == COOCCURENCE
        for edge in all_edges.values()
    )


def init_node_for_counting(netw_node):
    """
    IN USE!
    """
    node                     = {}
    node["id"]               = netw_node["id"]
    node[NEIGHBORS_NUMBER]   = 0
    node["neighbors"]        = set()
    node[COMPLEMENTS_NUMBER] = 0
    node["name"]             = netw_node["v"]["name"]
    node["species"]          = netw_node["v"]["taxonomy::species"]
    node["family"]           = netw_node["v"]["taxonomy::family"]
    node["order"]            = netw_node["v"]["taxonomy::order"]  # typo fixed in microbetag ...

    return node


def process_cooccurrence_and_regression(
    parsed_network, condition, metabolites, env_set
):
    """
    Plot number of a node's neighbors vs its seed complements per neighbor
    applying a Weighted Least Squares (WLS) model
    """

    # Initialize dictionary to track node interactions
    number_of_neighbors_pos = {}

    # Process edges for co-occurrence
    for edge in parsed_network.cx2.get_edges().values():
        if edge["v"]["interaction type"] == COMPL_INT_TYPE:
            id1, id2 = edge["s"], edge["t"]
            if check_cooccurrence(id1, id2, parsed_network.cx2):
                n1 = parsed_network.cx2.get_node(id1)
                n2 = parsed_network.cx2.get_node(id2)

                # Process interactions for each node pair
                for node, neighbor in [(n1, n2), (n2, n1)]:
                    if node["v"]["name"] in metabolites or node["v"]["name"] in env_set:
                        continue

                    if node["id"] not in number_of_neighbors_pos:
                        number_of_neighbors_pos[node["id"]] = init_node_for_counting(
                            node
                        )

                    number_of_neighbors_pos[node["id"]]["neighbors"].add(neighbor["id"])

                    # Count seed complements
                    for case, items in edge["v"].items():
                        if case.startswith("seedCompl"):
                            number_of_neighbors_pos[node["id"]][
                                COMPLEMENTS_NUMBER
                            ] += len(items)

    # Calculate neighbors count
    for counts in number_of_neighbors_pos.values():
        counts[NEIGHBORS_NUMBER] = len(counts["neighbors"])

    # Create dataframe
    df              = pd.DataFrame.from_dict(number_of_neighbors_pos, orient="index")
    df["condition"] = condition
    df[RATIO]       = df[COMPLEMENTS_NUMBER] / df[NEIGHBORS_NUMBER]

    return df


def get_node_ids_from_sp_name(sp_name, cx2):
    # NOTE (Haris Zafeiropoulos, 2025-04-02): NOT IN USE
    return [
        node["id"]
        for _, node in cx2.get_nodes().items()
        if node["v"]["name"] == sp_name
    ]


def get_sp_name_from_node_id(node_id, cx2):
    # NOTE (Haris Zafeiropoulos, 2025-04-02):  NOT IN USE
    for _, node in cx2.get_nodes().items():
        if node["id"] == node_id:
            return node["v"]["name"]


# ------------


def compute_weight_scores(cx_net, parsed_net):
    """
    Computes positive and negative weight scores for a given network.

    Parameters:
        cx_net: Network object containing edges.
        parsed_net: Parsed network dictionary with edge types.

    Returns:
        pos_weight_scores, neg_weight_scores: Dictionaries containing computed scores.
    """

    # Create a mapping of taxa pairs to edge IDs
    taxa_pair_2_edge_id = {
        (edge["s"], edge["t"]): edge_id
        for edge_id in cx_net.get_edges()
        if (edge := cx_net.get_edge(edge_id))["v"]["interaction type"]
        != COMPL_INT_TYPE
    }

    pos_weight_scores, neg_weight_scores = {}, {}

    for edge_id in parsed_net.edge_types["taxon_to_taxon"]["edges"]:
        edge = cx_net.get_edge(edge_id)

        if edge["v"]["interaction type"] in [COOCCURENCE, COEXCLUSION]:
            continue

        try:
            comp, coop, s, t = (
                edge["v"]["seed::competition"],
                edge["v"]["seed::cooperation"],
                edge["s"],
                edge["t"],
            )
        except Exception as e:
            print(f"WARNING: No seed scores for edge: {edge}")
            print(f"This would lead to ERROR: {e}")
            continue

        eid = taxa_pair_2_edge_id.get((s, t), taxa_pair_2_edge_id.get((t, s)))

        flashweave_score = cx_net.get_edge(eid)["v"]["microbetag::weight"]
        target           = pos_weight_scores if flashweave_score > 0 else neg_weight_scores

        target[edge_id] = {
            "cooperation": coop,
            "competition": comp,
            COOCCURENCE if flashweave_score > 0 else COEXCLUSION: flashweave_score,
        }

    return {"pos": pos_weight_scores, "neg": neg_weight_scores}


def drop_sample(abd_file=None, colnames=[], rownames=[], meta_file=None, abd_sep=',', meta_sep='\t'):

    if abd_file:
        abd = pd.read_csv(abd_file, sep=abd_sep)
        if len(colnames) > 0:
            abd = abd.drop(columns=colnames, errors="ignore")
            save_with_suffix(abd, abd_file, delim=abd_sep)

    if meta_file:
        meta = pd.read_csv(meta_file, sep=meta_sep)

        if len(colnames) > 0:
            meta = meta.drop(columns=colnames, errors="ignore")

        if len(rownames) > 0:
            meta = meta[~meta["sample"].isin(rownames)]
            save_with_suffix(meta, meta_file, delim=meta_sep)


def save_with_suffix(df, file_path, suffix="_cropped", delim=","):
    from pathlib import Path

    file_path = Path(file_path)
    new_path  = file_path.with_name(file_path.stem + suffix + file_path.suffix)
    df.to_csv(new_path, index=False, sep=delim)
    return new_path
