#!/usr/bin/env python3
import os
import sys
import json
import pandas as pd
from glob import glob
from ndex2.cx2 import RawCX2NetworkFactory

COOCCURRENCE        = "co-occurrence"
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

        self.cx2_nets        = MicrobetagNetsAnalysis.load_cx2_nets(path_to_cx2)
        self.taxonomic_level = taxonomic_level
        self.comp_threshold  = competition_threshold
        self.coop_threshold  = cooperation_threshold

        self.intersection    = None
        self.top_interactors = None
        self._interactors    = None

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

                s_node = cx2.get_node(edge["s"])
                t_node = cx2.get_node(edge["t"])
                s_type, s_taxa = MicrobetagNetsAnalysis.get_higher_taxonomic_level(s_node, self.taxonomic_level)
                t_type, t_taxa = MicrobetagNetsAnalysis.get_higher_taxonomic_level(t_node, self.taxonomic_level)

                # TOP INTERACTORS
                for type, taxon in [(s_type, s_taxa), (t_type, t_taxa)]:

                    if taxon not in per_net_interactors:

                        MicrobetagNetsAnalysis.init_node_counts(per_net_interactors, taxon, type)

                    if interaction is not None:

                        if isinstance(interaction, tuple):
                            for inter in interaction:
                                per_net_interactors[taxon][inter] += 1
                                per_net_interactors[taxon]["total_edges"] += 1
                        else:
                            per_net_interactors[taxon][interaction] += 1
                            per_net_interactors[taxon]["total_edges"] += 1

                # INTERSECTIONS
                key     = f"{s_taxa}-{t_taxa}"
                antikey = f"{t_taxa}-{s_taxa}"
                if key not in intersections:
                    if antikey in intersections:
                        key = antikey
                    else:
                        MicrobetagNetsAnalysis.init_edge_counts(intersections, key)

                if interaction == COOCCURRENCE or interaction == COEXCLUSION:
                    intersections[key]["occurrences_across_networks"] += 1

                curr_types = intersections[key]["interaction_types"]
                curr_nets  = intersections[key]["networks_with_the_edge"]

                curr_nets.add(filename)

                if interaction is not None:
                    if isinstance(interaction, tuple):
                        curr_types.update(interaction)
                    else:
                        curr_types.add(interaction)

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
    def load_cx2_nets(paths):

        cx2_nets = set()

        # Normalize to list
        if isinstance(paths, (str, os.PathLike)):
            paths = [paths]
        elif not isinstance(paths, (list, tuple, set)):
            raise ValueError(f"Invalid type for paths: {type(paths)}")

        # Expand all paths and collect .cx2 files
        all_files = []
        for p in paths:
            for f in glob(p, recursive=True):  # handles wildcards & directories
                if os.path.isfile(f) and f.endswith(".cx2"):
                    all_files.append(f)
                elif os.path.isdir(f):
                    # recursively add all .cx2 files in directory
                    all_files.extend(
                        os.path.join(root, file)
                        for root, _, files in os.walk(f)
                        for file in files
                        if file.endswith(".cx2")
                    )

        # Load all files
        for file_path in all_files:
            print(f"Loading {file_path}...")
            cx2_nets.add(MicrobetagNetsAnalysis.load_cx2(file_path))  # or self.load_cx2()

        return cx2_nets

    @staticmethod
    def init_node_counts(dict, new_key, type):
        dict[new_key] = {
            "type"       : type,
            "total_edges": 0,
            COOCCURRENCE  : 0,
            COEXCLUSION  : 0,
            COMPETITION  : 0,
            COOPERATION  : 0
        }
        return dict

    @staticmethod
    def init_edge_counts(dict, new_key):
        dict[new_key] = {
            "occurrences_across_networks": 0,
            "interaction_types"          : set(),
            "networks_with_the_edge"     : set()
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

        if interaction == COOCCURRENCE:
            return COOCCURRENCE
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

        elif node["v"]['microbetag::ncbi-tax-level'] not in levels:
            metavar = node["v"]["name"]
            return "metavar", metavar

        elif node["v"]['microbetag::ncbi-tax-level'] == taxonomic_level:

            taxon_full = node["v"].get("microbetag::taxon", "")
            # partition() always returns three parts: (before, sep, after):
            _, sep, after = taxon_full.partition("__")
            taxon = after if sep else taxon_full

            taxon = parse_taxonomy_underscores(taxon_full)

            return "taxon", taxon

        else:

            # taxa     = [x.split("__")[1] for x in node["v"]['microbetag::taxonomy'].split(";")]
            taxa     = [parse_taxonomy_underscores(x) for x in node["v"]['microbetag::taxonomy'].split(";")]
            tax_dict = dict(zip(levels, taxa))

            current = node["v"]['microbetag::ncbi-tax-level']

            current_level = levels[levels.index(current)]
            desired_level = levels[levels.index(taxonomic_level)]

            if levels.index(current_level) < levels.index(desired_level):
                return "taxon", tax_dict.get(current_level, None)
            else:
                return "taxon", tax_dict.get(taxonomic_level, None)

def parse_taxonomy_underscores(taxonomy):
    _, sep, after = taxonomy.partition("__")
    taxon = after if sep else taxonomy
    return taxon


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Microbetag Networks Analysis")
    # # NOTE (Haris Zafeiropoulos, 2025-12-03): nargs means one or more arguments can be passed after -p
    parser.add_argument("-p", "--path", type=str, nargs="+", required=True, help="Path to the directory containing the cx2 files. Can be either a folder or a single cx2 file.")

    parser.add_argument("-l", "--taxonomic-level", type=str, default='species', help="Taxonomic level to analyze (default: species).")
    parser.add_argument("-comp", "--competition_threshold", type=float, default=0.5, help="Threshold for competition interaction (default: 0.5).")
    parser.add_argument("-coop", "--cooperation_threshold", type=float, default=0.2, help="Threshold for cooperation interaction (default: 0.2).")
    parser.add_argument("-r", "--output_interactors", type=str, default="interactors.csv", help="Output CSV file for top interactors (default: interactors.csv).")
    parser.add_argument("-s", "--output_intersections", type=str, default="intersections.csv", help="Output CSV file for intersections (default: intersections.csv).")
    parser.add_argument("-f", "--format", type=str, choices=['csv', 'html'], default='csv', help="Output file format: 'csv' or 'html' (default: csv).")
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
