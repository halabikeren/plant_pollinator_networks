from ete3 import Tree
import click
import os
import pandas as pd
import numpy as np
from collections import Counter

def process_tree(tree_path: str) -> Tree:
    tree = Tree(tree_path, format=1)
    for leaf in tree.get_leaves():
        leaf.name = leaf.name.replace("_"," ").lower()
    return tree

@click.command()
@click.option(
    "--input_path",
    help="path to a network",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--output_dir",
    help="path to the phylogenetic features files",
    type=click.Path(exists= False),
    required=True,
)
@click.option(
    "--tree_path",
    help="path to a tree",
    type=click.Path(exists=True),
    required=False,
    default="/groups/itay_mayrose/halabikeren/plant_pollinator_networks/trees/ALLOTB_expanded_by_unresolved_names.tre",
)
@click.option(
    "--name_resolution_path",
    help="path to name resolution of plant names",
    type=click.Path(exists=True),
    required=False,
    default="/groups/itay_mayrose/halabikeren/plant_pollinator_networks/name_resolution/resolved_plant_names.csv",
)
@click.option(
    "--use_name_resolution",
    help="indicator of weather the resolved names should be projected on the tree or the unresolved names",
    type=bool,
    required=False,
    default=False,
)
def compute_phylogenetic_features(input_path: str,
                                  output_dir: str,
                                  tree_path: str,
                                  name_resolution_path: str,
                                  use_name_resolution: str):
   
    os.makedirs(output_dir, exist_ok=True)
    network = pd.read_csv(input_path)
    unresolved_plant_names = network.Plant.str.lower().tolist()
    plant_features = pd.DataFrame(columns=["species", "resolved_name", "phylogenetic_uniqueness"])
    plant_features["species"] = unresolved_plant_names
    name_field = "species"

    if use_name_resolution:
        name_resolution = pd.read_csv(name_resolution_path)
        name_resolution.original_name = name_resolution.original_name.str.lower()
        name_resolution.resolved_name = name_resolution.resolved_name.str.lower()
        plant_features.set_index("species", inplace=True)
        plant_features["resolved_name"].fillna(value=name_resolution.set_index("original_name")["resolved_name"].to_dict(), inplace=True)
        plant_features.reset_index(inplace=True)
        name_field = "resolved_name"

    network_tree = process_tree(tree_path)
    plant_names = plant_features[name_field].dropna().tolist()
    if len(plant_names) == 0:
        raise ValueError(f"no plant names are available and so phylogenetic features cannot be computed")
    tree_names = set(network_tree.get_leaf_names())
    plant_names_in_tree = [name for name in plant_names if name in tree_names]
    if len(plant_names_in_tree) == 0:
        raise ValueError(f"no resolved plant names are present in the tree so phylogenetic features cannot be computed")
    network_tree.prune(plant_names_in_tree, preserve_branch_length=True)
    pd_val = np.sum([node.dist for node in network_tree.traverse()])
    network_features = pd.DataFrame({"network": [os.path.basename(input_path)],
                                     "num_network_species": len(plant_names),
                                     "num_network_species_in_tree": len(plant_names_in_tree),
                                     "coverage": len(plant_names_in_tree)/len(plant_names),
                                     "pd": [pd_val]})
    network_features.to_csv(f"{output_dir}/network_features.csv", index=False)

    # phylogenetic uniqueness per plant species -  distance from other plant species in the network across the tree
    leaf_names = set(network_tree.get_leaf_names())
    plant_features["phylogenetic_uniqueness"] = plant_features[name_field].apply(lambda sp1: np.mean([network_tree.get_distance(sp1, sp2) for sp2 in leaf_names if sp1 != sp2]) if sp1 in leaf_names else np.nan)
    plant_features["genus"] = plant_features[name_field].apply(lambda name: name.split(" ")[0])
    genera_counter = Counter(plant_features["genus"].tolist())
    plant_features["has_relative_in_network"] = plant_features["genus"].apply(lambda genus: genera_counter[genus] > 1)
    plant_features.to_csv(f"{output_dir}/plant_features.csv", index=False)
    
if __name__ == '__main__':
    compute_phylogenetic_features()