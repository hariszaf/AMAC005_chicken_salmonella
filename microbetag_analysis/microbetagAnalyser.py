#!/usr/bin/env python3
import os
import json
import pandas as pd
# pip install ndex2
from ndex2.cx2 import RawCX2NetworkFactory

COOCCURENCE        = "co-occurrence"
COEXCLUSION        = "co-exclusion"
COMPETITION        = "competition"
COOPERATION        = "cooperation"
COMPL_INT_TYPE     = "complementarity"
SEED_COMP          = "::".join(["seed", COMPETITION])
SEED_COOP          = "::".join(["seed", COOPERATION])

class MicrobetagNetsAnalysis:

    """
    Usage example:
    >>> anal = MicrobetagNetsAnalysis("microbetag_nets/prev20/per_day/")  # path to cx2 files or directory
    >>> interactors, interscetion = anal.get_stats_df()
    >>> interscetion.head()
    >>> interactors.head()
    """

    def __init__(self, path_to_cx2, taxonomic_level='species', competition_threshold=0.5, cooperation_threshold=0.2):

        self.cx2_nets        = self.load_cx2_nets(path_to_cx2)
        self.taxonomic_level = taxonomic_level
        self.comp_threshold  = competition_threshold
        self.coop_threshold  = cooperation_threshold

        self.intersection    = None
        self.top_interactors = None
        self._interactors    = None

    def load_cx2_nets(self, path_to_cx2):

        cx2_nets = set()

        if os.path.isfile(path_to_cx2):
            print(f"Loading {path_to_cx2}...")
            cx2_nets.add(MicrobetagNetsAnalysis.load_cx2(path_to_cx2))   # or self.load_cx2()

        elif os.path.isdir(path_to_cx2):

            for filename in os.listdir(path_to_cx2):
                full_path = os.path.join(path_to_cx2, filename)

                if os.path.isfile(full_path):
                    print(f"Loading {full_path}...")
                    cx2_nets.add(MicrobetagNetsAnalysis.load_cx2(full_path))   # or self.load_cx2()
        return cx2_nets

    def get_stats(self):
        """
        Interactors
        ----------

        Node             Total edges   Cooccurrences   Mutual exclusions   Complementarities   Competitions
        --------------------------------------------------------------------------------------------------
        Bacteroides           18             14                 2                   1                 1
        Prevotella            15             12                 1                   1                 1

        Intersections
        ------------

        Edge                       Number across networks    Interaction types
        ----------------------------------------------------------------------
        Bacteroides<->Prevotella         12                  [COOCCURRENCE]
        Escherichia<->Acinetobacter       9                  [COOCCURRENCE, COMPETITION]
        """

        interactors   = {}
        intersections = {}

        for filename, cx2 in self.cx2_nets:

            per_net_interactors = {}

            edges = cx2.get_edges()

            for _, edge in edges.items():

                interaction = MicrobetagNetsAnalysis.get_interaction_types(
                    edge, self.comp_threshold, self.coop_threshold
                )

                if interaction is not None:
                    s_node = cx2.get_node(edge["s"])
                    t_node = cx2.get_node(edge["t"])
                    s_type, s_taxa = MicrobetagNetsAnalysis.get_higher_taxonomic_level(s_node, self.taxonomic_level)
                    t_type, t_taxa = MicrobetagNetsAnalysis.get_higher_taxonomic_level(t_node, self.taxonomic_level)

                    # TOP INTERACTORS
                    for type, taxon in [(s_type, s_taxa), (t_type, t_taxa)]:

                        if taxon not in per_net_interactors:

                            per_net_interactors = MicrobetagNetsAnalysis.init_node_counts(
                                per_net_interactors, taxon, type
                            )

                        per_net_interactors[taxon]["total_edges"] += 1

                        if isinstance(interaction, tuple):
                            for inter in interaction:
                                per_net_interactors[taxon][inter] += 1
                        else:
                            per_net_interactors[taxon][interaction] += 1

                    # INTERSECTIONS
                    key     = f"{s_taxa}-{t_taxa}"
                    antikey = f"{t_taxa}-{s_taxa}"
                    if key not in intersections:
                        if antikey in intersections:
                            key = antikey
                        else:
                            intersections = MicrobetagNetsAnalysis.init_edge_counts(intersections, key)
                    if interaction == COOCCURENCE or interaction == COEXCLUSION:
                        intersections[key]["occurrences_across_networks"] += 1
                    curr = intersections[key]["interaction_types"]
                    if isinstance(interaction, tuple):
                        curr.update(interaction)
                    else:
                        curr.add(interaction)

            interactors[filename] = per_net_interactors

        top_interactors = {}
        for net in interactors:
            for taxon, counts in interactors[net].items():
                if taxon not in top_interactors:
                    top_interactors[taxon] = counts
                else:
                    for key in counts.keys() - {"type"}:
                        top_interactors[taxon][key] += counts[key]

        self._interactors    = interactors
        self.top_interactors = top_interactors
        self.intersection    = intersections

    def get_stats_df(self):

        if self.top_interactors is None or self.intersection is None:
            self.get_stats()

        else:
            self.intersection    = None
            self.top_interactors = None
            self._interactors    = None
            print(f"Recomputing stats for taxonomic level: {self.taxonomic_level}...")
            self.get_stats()

        intersection = pd.DataFrame.from_dict(self.intersection).T.sort_values(by="occurrences_across_networks", ascending=False)
        interactors  = pd.DataFrame.from_dict(self.top_interactors).T.sort_values(by="total_edges", ascending=False)

        return interactors, intersection

    @staticmethod
    def init_node_counts(dict, new_key, type):
        dict[new_key] = {
            "type"       : type,
            "total_edges": 0,
            COOCCURENCE  : 0,
            COEXCLUSION  : 0,
            COMPETITION  : 0,
            COOPERATION  : 0
        }
        return dict

    @staticmethod
    def init_edge_counts(dict, new_key):
        dict[new_key] = {
            "occurrences_across_networks": 0,
            "interaction_types"          : set()
        }
        return dict

    @staticmethod
    def load_cx2(path_to_cx2):

        with open(path_to_cx2, 'r') as cx2_file:
            cx2_dict = json.load(cx2_file)

        factory = RawCX2NetworkFactory()
        return path_to_cx2.split("/")[-1].split(".cx2")[0], factory.get_cx2network(cx2_dict)

    @staticmethod
    def get_interaction_types(edge, comp_threshold=0.5, coop_threshold=0.2):

        interaction = edge["v"]['interaction type']

        if interaction == COOCCURENCE:
            return COOCCURENCE
        elif interaction == COEXCLUSION:
            return COEXCLUSION
        else:
            try:
                comp = edge["v"][SEED_COMP]
                coop = edge["v"][SEED_COOP]
            except KeyError as e:
                comp, coop = None, None
                pass

            if comp is not None and coop is not None:

                if comp > comp_threshold and coop < coop_threshold:
                    return COMPETITION

                elif comp > coop_threshold and coop > coop_threshold:
                    return COMPETITION, COOPERATION

                elif comp < comp_threshold and coop > coop_threshold:
                    return COOPERATION
                else:
                    return None

    @staticmethod
    def get_higher_taxonomic_level(node, taxonomic_level='species'):

        levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',]

        if 'microbetag::ncbi-tax-level' not in node["v"]:
            metavar = node["v"]["name"]
            return "metavar", metavar

        elif node["v"]['microbetag::ncbi-tax-level'] == taxonomic_level:
            return "taxon", node["v"]["microbetag::taxon"].split("__")[1]

        else:

            taxa     = [x.split("__")[1] for x in node["v"]['microbetag::taxonomy'].split(";")]
            tax_dict = dict(zip(levels, taxa))

            current = node["v"]['microbetag::ncbi-tax-level']

            current_level = levels[levels.index(current)]
            desired_level = levels[levels.index(taxonomic_level)]

            if levels.index(current_level) < levels.index(desired_level):
                return "taxon", tax_dict.get(current_level, None)
            else:
                return "taxon", tax_dict.get(taxonomic_level, None)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Microbetag Networks Analysis")
    parser.add_argument("--path", type=str, required=True, help="Path to the directory containing the cx2 files. Can be either a folder or a single cx2 file.")
    parser.add_argument("--taxonomic_level", type=str, default='species', help="Taxonomic level to analyze (default: species).")
    parser.add_argument("--competition_threshold", type=float, default=0.5, help="Threshold for competition interaction (default: 0.5).")
    parser.add_argument("--cooperation_threshold", type=float, default=0.2, help="Threshold for cooperation interaction (default: 0.2).")
    parser.add_argument("--output_interactors", type=str, default="interactors.csv", help="Output CSV file for top interactors (default: interactors.csv).")
    parser.add_argument("--output_intersections", type=str, default="intersections.csv", help="Output CSV file for intersections (default: intersections.csv).")
    parser.add_argument("--format", type=str, choices=['csv', 'html'], default='csv', help="Output file format: 'csv' or 'html' (default: csv).")
    args = parser.parse_args()

    anal = MicrobetagNetsAnalysis(
        args.path,
        taxonomic_level       = args.taxonomic_level,
        competition_threshold = args.competition_threshold,
        cooperation_threshold = args.cooperation_threshold
    )
    interactors, intersection = anal.get_stats_df()

    if args.format == 'csv':
        interactors.to_csv(args.output_interactors.split(".")[0] + ".csv")
        intersection.to_csv(args.output_intersections.split(".")[0] + ".csv")
    elif args.format == 'html':
        interactors.to_html(args.output_interactors.split(".")[0] + ".html")
        intersection.to_html(args.output_intersections.split(".")[0] + ".html")
    else:
        print("Unsupported format. Please use 'csv' or 'html'.")
