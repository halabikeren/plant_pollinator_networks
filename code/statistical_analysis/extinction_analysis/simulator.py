import os.path
import random
from enum import Enum
import pandas as pd
import numpy as np


class ExtinctionOrder(Enum):
    random = 0
    polyploids_first = 1
    diploids_first = 2


class Simulator:
    classification_data: dict
    ext_order: ExtinctionOrder
    rewiring_flag: bool
    rewiring_probability: float

    def __init__(self,
                 classification_path: str,
                 ext_order: ExtinctionOrder = ExtinctionOrder.random,
                 rewiring_flag: bool = True,
                 rewiring_probability: float = 0):
        classification_data = pd.read_csv(classification_path)
        ploidy_level_col = "is_polyploid_by_resolved"
        if "Plant" in classification_data.columns:
            classification_data["original_name"] = classification_data.Plant.str.lower()
        else:
            classification_data.original_name = classification_data.original_name.str.lower()
        classification_data[ploidy_level_col].fillna(-1, inplace=True)
        self.classification_data = classification_data.set_index("original_name")[ploidy_level_col].to_dict()
        self.ext_order = ext_order
        self.rewiring_flag = rewiring_flag
        self.rewiring_probability = rewiring_probability

    def cascade(self,
                start: str,
                network: pd.DataFrame,
                primary_extinction_candidates: list) -> tuple[pd.DataFrame, pd.DataFrame]:
        already_extinct = network.columns[np.where(network.sum(axis=0) == 0)]
        network.loc[start] = 0
        cascade_data = pd.DataFrame({"cascade_iteration": [1],
                                     "extinction_type": ["primary"],
                                     "extinction_level": ["plant"],
                                     "extinct_taxon": [start]})
        total_extinct = network.columns[np.where(network.sum(axis=0) == 0)]
        second_extinct = [sp for sp in total_extinct if sp not in already_extinct]
        remaining_plants = list(set((network.sum(axis=1) > 0).index) & set(primary_extinction_candidates))
        if self.rewiring_flag and len(remaining_plants) > 0:
            network[second_extinct] = 0
            rewiring_candidates = random.sample(second_extinct, int(len(second_extinct) * self.rewiring_probability))
            for rewiring_candidate in rewiring_candidates:
                switch_taxon = random.sample(remaining_plants, 1)
                network.loc[switch_taxon, rewiring_candidate] = 1
                second_extinct.remove(rewiring_candidate)
        cascade_step = pd.DataFrame({"extinct_taxon": second_extinct})
        cascade_step["cascade_iteration"] = 2
        cascade_step["extinction_type"] = "secondary"
        cascade_step["extinction_level"] = "pollinator"
        cascade_step["extinct_taxon"] = second_extinct
        cascade_data = pd.concat([cascade_data, cascade_step])
        return cascade_data, network

    def simulate(self,
                 network_path: str) -> pd.DataFrame:
        network = pd.read_csv(network_path)
        if "Plant" not in network.columns:
            network = network.rename(columns={"Unnamed: 0": "Plant"})
        if "Unnamed: 0" in network.columns:
            network = network.drop(["Unnamed: 0"], axis=1)
        network.Plant = network.Plant.str.lower()
        network.set_index("Plant", inplace=True)
        plants = network.index.str.lower().tolist()
        primary_extinction_order = plants
        random.shuffle(primary_extinction_order)
        polyploids = [p for p in plants if self.classification_data.get(p, -1) == 1]
        random.shuffle(polyploids)
        diploids = [p for p in plants if self.classification_data.get(p, -1) == 0]
        random.shuffle(diploids)
        unknown = list(set(plants) - set(polyploids) - set(diploids))
        random.shuffle(unknown)
        if self.ext_order == ExtinctionOrder.polyploids_first:
            primary_extinction_order = polyploids + diploids + unknown
        if self.ext_order == ExtinctionOrder.diploids_first:
            primary_extinction_order = diploids + polyploids + unknown
        extinction_data = []
        i = 0
        assert (len(primary_extinction_order) == network.shape[0])
        network_copy = network.copy()
        while network_copy.sum().sum() > 0:
            if len(primary_extinction_order) == 0:
                print(
                    f"no species left for primary extinction, but there are somehow still {network.sum().sum()} interactions in the network {network_path}")
                return pd.DataFrame(columns=["cascade_iteration", "extinction_type", "extinction_level",
                                             "extinct_taxon", "primary_iteration", "simulation_index"])
            candidate = primary_extinction_order.pop(0)
            cascade_out, network_copy = self.cascade(start=candidate,
                                                     network=network_copy,
                                                     primary_extinction_candidates=primary_extinction_order)
            i += 1
            cascade_out["primary_iteration"] = i
            extinction_data.append(cascade_out)
        extinction_data = pd.concat(extinction_data)
        return extinction_data

    def write_simulations(self,
                          network_path: str,
                          output_path: str,
                          nsim: int = 100):

        if not os.path.exists(output_path):
            full_extinction_data = []
            for i in range(nsim):
                extinction_data = self.simulate(network_path=network_path)
                extinction_data["simulation_index"] = i
                full_extinction_data.append(extinction_data)
            full_extinction_data = pd.concat(full_extinction_data)
            full_extinction_data.to_csv(output_path)

if __name__ == '__main__':
    simulator = Simulator(classification_path="/groups/itay_mayrose/halabikeren/plant_pollinator_networks/data/ploidy_classification/plant_classification.csv",
                          ext_order = ExtinctionOrder.diploids_first,
                          rewiring_flag = True,
                          rewiring_probability = 0)
    simulator.write_simulations(network_path="/groups/itay_mayrose/halabikeren/plant_pollinator_networks/data/networks/all/weighted/470.csv",
                                output_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/data/statistical_analysis/extinction_analysis/rewiring_prob_0/weighted/diploids_first/470.csv",
                                nsim=5)