from dataclasses import dataclass

# To parse network
#    pip install ndex2
from ndex2.cx2 import RawCX2NetworkFactory

import venn

import statsmodels.api as sm

import numpy as np
import pandas as pd

from plotnine import ggplot, aes, geom_point, theme_minimal, labs, ggsave, geom_line, annotate



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
        "taxon_to_taxon": {"edges": [], "counts": 0},  # done
        "taxon_to_metabo": {"edges": [], "counts": 0},
        "taxon_to_env": {"edges": [], "counts": 0},
        "metabo_to_metabo": {"edges": [], "counts": 0}, # done
        "env_to_env": {"edges": [], "counts": 0},       # done
        "metabo_to_env": {"edges": [], "counts": 0}
    }

    source_target_edge_ids = {}
    taxa_pairs_with_positive_weight = []

    for _, edge in cx2.get_edges().items():

        edge_id = edge["id"]

        if "completes/competes with" == edge["v"]["interaction type"]:
            for k in edge["v"]:
                # NOTE For PATHWAY COMPLEMENTS
                # if "compl::" in k:
                #     edge_id_2_number_of_complements[edge_id] = len(edge["v"][k])
                uniq_cross_feeding_compounds = set()
                if k.startswith("seedCompl::"):
                    edge_id_2_number_of_complements[edge_id] = len(edge["v"][k])
                    for compl in edge["v"][k]:
                        for i in compl.split("^")[2].split(";"): uniq_cross_feeding_compounds.add(i)
                    edge_id_2_cross_feeding_compounds[edge_id] = uniq_cross_feeding_compounds
        else:
            if 'microbetag::weight' in edge["v"]:
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

        if source_name not in env_set | metabolites_set and target_name not in env_set | metabolites_set:
            edge_types["taxon_to_taxon"]["counts"] += 1
            edge_types["taxon_to_taxon"]["edges"].append(edge_id)

        elif source_name in env_set and target_name in env_set:
            edge_types["env_to_env"]["counts"] += 1
            edge_types["env_to_env"]["edges"].append(edge_id)

        elif source_name in metabolites_set and target_name in metabolites_set:
            edge_types["metabo_to_metabo"]["counts"] += 1
            edge_types["metabo_to_metabo"]["edges"].append(edge_id)

        elif {source_name, target_name} <= env_set | metabolites_set and source_name in env_set ^ metabolites_set:
            edge_types["metabo_to_env"]["counts"] += 1
            edge_types["metabo_to_env"]["edges"].append(edge_id)

        elif {source_name, target_name} & env_set and not {source_name, target_name} & metabolites_set:
            edge_types["taxon_to_env"]["counts"] += 1
            edge_types["taxon_to_env"]["edges"].append(edge_id)

        elif {source_name, target_name} & metabolites_set and not {source_name, target_name} & env_set:
            edge_types["taxon_to_metabo"]["counts"] += 1
            edge_types["taxon_to_metabo"]["edges"].append(edge_id)

        else:
            print("PROBLEM:", source_name, target_name)

    return edge_id_2_number_of_complements, edge_types, taxa_pairs_with_positive_weight, edge_id_2_cross_feeding_compounds


def get_compls_and_compounds_in_positive_associated_taxa(cx2, positive_pairs, edge_id_compl_number, compounds):

    number_of_complements_in_positive_associated_taxa = 0
    unique_compounds = set()
    for _, edge in cx2.get_edges().items():

        s_t = (edge["s"], edge["t"])
        t_s = (edge["t"], edge["s"])
        edge_id = edge["id"]

        if s_t in positive_pairs or t_s in positive_pairs:
            if edge["v"]["interaction type"] == "completes/competes with":
                for k in edge["v"]:
                    if k.startswith("seedCompl::"):
                        number_of_complements_in_positive_associated_taxa += edge_id_compl_number[edge_id]
                        for i in compounds[edge_id]:
                            unique_compounds.add(i)

    return number_of_complements_in_positive_associated_taxa, unique_compounds



@dataclass
class MggParser():
    """A `dataclass` for the parsed microbetag annotated network once loaded with `ndex2`; instances of this class will store outputs of the `parse_microbetag_egdes()`."""

    import ndex2

    num_total_compls: dict
    edge_types: dict
    pos_pairs_node_ids: list  # of tuples
    edge_id_to_unique_compounds: dict
    cx2: ndex2.cx2.CX2Network

    def __init__(self, cx2, env_set, metabolites_set) -> None:

        (
            self.num_total_compls,
            self.edge_types,
            self.pos_pairs_node_ids,
            self.edge_id_to_unique_compounds
         ) = parse_microbetag_edges(cx2, env_set, metabolites_set)
        self.cx2 = cx2


def check_cooccurrence(s, t, cx2):
    """
    IN USE!
    """
    all_edges = cx2.get_edges()

    # Check for either direction (source -> target or target -> source)
    return any(
        (
         (edge["s"] == s and edge["t"] == t)
         or
         (edge["s"] == t and edge["t"] == s)
        )
        and edge["v"]['interaction type'] == "cooccurrence"
        for edge in all_edges.values()
    )


def init_node_for_counting(netw_node):
    """
    IN USE!
    """
    node  = {}
    node["id"] = netw_node["id"]
    node["neighbors_num"] = 0
    node["neighbors"] = set()
    node["seed_compls_num"] = 0
    node["name"] = netw_node["v"]["name"]
    node["species"] = netw_node["v"]["taxonomy::species"]
    node["family"] = netw_node["v"]["taxonomy::family"]
    node["order"] = netw_node["v"]["taxonomy::oder"]      # typo fixed in microbetag ...

    return node


def process_cooccurrence_and_regression(day_net, day, metabolites, env_set, model="WLS"):
    """
    Plot number of a node's neighbors vs its seed complements per neighbor
    applying a Weighted Least Squares (WLS) model
    """

    # Initialize dictionary to track node interactions
    number_of_neighbors_pos = {}

    # Process edges for co-occurrence
    for edge in day_net.cx2.get_edges().values():
        if edge['v']['interaction type'] == 'completes/competes with':
            id1, id2 = edge["s"], edge["t"]
            if check_cooccurrence(id1, id2, day_net.cx2):
                n1, n2 = day_net.cx2.get_node(id1), day_net.cx2.get_node(id2)

                # Process interactions for each node pair
                for node, neighbor in [(n1, n2), (n2, n1)]:
                    if node["v"]["name"] in metabolites or node["v"]["name"] in env_set:
                        continue

                    if node["id"] not in number_of_neighbors_pos:
                        number_of_neighbors_pos[node["id"]] = init_node_for_counting(node)

                    number_of_neighbors_pos[node["id"]]["neighbors"].add(neighbor["id"])

                    # Count seed complements
                    for case, items in edge["v"].items():
                        if case.startswith("seedCompl"):
                            number_of_neighbors_pos[node["id"]]["seed_compls_num"] += len(items)

    # Calculate neighbors count
    for counts in number_of_neighbors_pos.values():
        counts["neighbors_num"] = len(counts["neighbors"])

    # Create dataframe
    df = pd.DataFrame.from_dict(number_of_neighbors_pos, orient="index")
    df["day"] = day
    df["compl_ratio"] = df["seed_compls_num"] / df["neighbors_num"]

    if model == "WLS":
        regression_df, r_squared, p_value = wls(df)

    elif model == "OLS":
        regression_df, r_squared, p_value = ols(df)

    # Generate the plot
    plot = (
        ggplot(df, aes(x='neighbors_num', y='seed_compls_num', color='order')) +
        geom_point(size=3) +
        geom_line(regression_df, aes(x="neighbors_num", y="compl_ratio"), color="red") +
        theme_minimal() +
        labs(title=f'{day} at the family level',
            x='Neighbors Number',
            y='Seed Complements Number') +
        annotate("text", x=df["neighbors_num"].max() * 0.8, y=df["compl_ratio"].max() * 0.9,
                 label=f'R² = {r_squared:.3f}\nP-value = {p_value:.3g}', size=10, color="black")
    )

    return plot, df



def wls(df):
    """
                    Weighted Least Squares (WLS) Regression

    """

    # Fit Weighted Least Squares (WLS) model
    weights = 1 / df['neighbors_num'].value_counts().reindex(df['neighbors_num']).values
    X = sm.add_constant(df["neighbors_num"])
    y = df["compl_ratio"]
    model = sm.WLS(y, X, weights=weights).fit()

    # Create regression line for plot
    x_range = np.linspace(df["neighbors_num"].min(), df["neighbors_num"].max(), 100)
    X_pred = sm.add_constant(x_range)
    y_pred = model.predict(X_pred)
    wls_regression_df = pd.DataFrame({"neighbors_num": x_range, "compl_ratio": y_pred})

    # Compute R² and p-value
    r_squared = model.rsquared
    p_value = model.pvalues["neighbors_num"]

    return wls_regression_df, r_squared, p_value


def ols(df):
    """
                                OLS Regression

    """

    # Fit Ordinary Least Squares (OLS) model
    X = sm.add_constant(df["neighbors_num"])
    y = df["compl_ratio"]
    model = sm.OLS(y, X).fit()

    # Create regression line for plot
    x_range = np.linspace(df["neighbors_num"].min(), df["neighbors_num"].max(), 100)
    X_pred = sm.add_constant(x_range)
    y_pred = model.predict(X_pred)
    ols_regression_df = pd.DataFrame({"neighbors_num": x_range, "compl_ratio": y_pred})

    # Compute R² and p-value
    r_squared = model.rsquared
    p_value = model.pvalues["neighbors_num"]

    return ols_regression_df, r_squared, p_value


def get_node_ids_from_sp_name(sp_name, cx2):
    # NOTE (Haris Zafeiropoulos, 2025-04-02): NOT IN USE
    return [node["id"] for _, node in cx2.get_nodes().items() if node["v"]["name"] == sp_name]


def get_sp_name_from_node_id(node_id, cx2):
    # NOTE (Haris Zafeiropoulos, 2025-04-02):  NOT IN USE
    for _, node in cx2.get_nodes().items():
        if node["id"] == node_id:
            return node["v"]["name"]


def fit_polynomial_regression(df, degree=2):
    # Create polynomial features
    X = np.vander(df["neighbors_num"], N=degree + 1, increasing=True)

    # Fit the polynomial regression model
    y = df["compl_ratio"]
    model = sm.OLS(y, X).fit()

    # Create regression curve
    x_range = np.linspace(df["neighbors_num"].min(), df["neighbors_num"].max(), 100)
    X_pred = np.vander(x_range, N=degree + 1, increasing=True)
    y_pred = model.predict(X_pred)

    poly_regression_df = pd.DataFrame({"neighbors_num": x_range, "compl_ratio": y_pred})

    # Compute R² and p-values
    r_squared = model.rsquared
    p_values = model.pvalues  # p-values for each coefficient

    return poly_regression_df, r_squared, p_values