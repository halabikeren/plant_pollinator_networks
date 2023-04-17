import pandas as pd
import os
import numpy as np
import random
from enum import Enum


class ExtinctionOrder(Enum):
    random: 0
    polyploids_first: 1
    diploids_first: 2


def generate_network_extinction_simulations(network_path: str,
                                            network_species_metedata_path: str,
                                            output_path: str,
                                            polyploid_range: tuple = (0.9, 0.9),
                                            diploid_range: tuple = (0.9, 0.9),
                                            unknown_range: tuple = (0.9, 0.9),
                                            pollinator_range: tuple = (0.9, 0.9),
                                            nsim=5):
    if not os.path.exists(output_path):
        network = pd.read_csv(network_path, encoding='unicode_escape')
        network_metedata = pd.read_csv(network_species_metedata_path).query(
            f"network == '{os.path.basename(network_path)}'")
        network.rename(columns={"plant/pollinator": "Plant"}, inplace=True)
        if "Plant" not in network.columns:
            network.rename(columns={network.columns.tolist()[0]: "Plant"}, inplace=True)
        processed_network = network.melt(id_vars=["Plant"],
                                              var_name="Pollinator",
                                              value_name="Value")
        processed_network['Plant'] = processed_network['Plant'].str.lower()  # lowercase
        processed_network['Pollinator'] = processed_network['Pollinator'].str.lower()  # lowercase
        processed_network['Value'] = pd.to_numeric(processed_network['Value'], errors='coerce')
        processed_network = processed_network[processed_network["Value"] > 0]
        processed_network = processed_network.merge(network_metedata[["original_name", "is_polyploid"]],
                                                    how='left',
                                                    left_on='Plant',
                                                    right_on='original_name')
        processed_network['is_polyploid'].fillna(-1, inplace=True)
        extinction_dfs = []
        for i in range(nsim):
            extinction_df = extinction_analysis(classified_network=processed_network,
                                                polyploid_range=polyploid_range,
                                                diploid_range=diploid_range,
                                                unknown_range=unknown_range,
                                                pollinator_range=pollinator_range)
            extinction_df["simulation_index"] = i
            extinction_dfs.append(extinction_df)
        all_extinction_data = pd.concat(extinction_dfs)
        all_extinction_data.to_csv(output_path)


def All_networks_extinction_analysis(networks_data: dict,
                                     output_dir: str,
                                     range_of_rates: tuple,
                                     num_sim: int = 100):
    network_folder = networks_data['networks']
    species_data_path = networks_data['data']
    networks_files = {f: os.path.join(network_folder, f) for f in os.listdir(network_folder) if
                      os.path.isfile(os.path.join(network_folder, f))}
    for i, network in enumerate(networks_files):
        network_path = networks_files[network]
        if not "reference" in network_path and network_path.endswith(".csv"):
            generate_network_extinction_simulations(network_path,
                                                    network_species_metedata_path=species_data_path,
                                                    polyploid_range=range_of_rates,
                                                    diploid_range=range_of_rates,
                                                    unknown_range=range_of_rates,
                                                    pollinator_range=range_of_rates,
                                                    nsim=num_sim,
                                                    output_path=f"{output_dir}/extinction_simulations_network_{i}.csv")
    return 0


def update_edgelist_totals(ep_network_edgelist):
    ep_network_edgelist.loc[:, 'total_per_plant'] = ep_network_edgelist.groupby(['Plant'])['Value'].transform(
        lambda x: x.sum())
    ep_network_edgelist.loc[:, 'total_per_pollinator'] = ep_network_edgelist.groupby(['Pollinator'])['Value'].transform(
        lambda x: x.sum())
    ep_network_edgelist.loc[:, 'd_plant_poli'] = ep_network_edgelist['Value'] / ep_network_edgelist['total_per_plant']
    ep_network_edgelist.loc[:, 'd_poli_plant'] = ep_network_edgelist['Value'] / ep_network_edgelist[
        'total_per_pollinator']
    return ep_network_edgelist


def get_config(name: str,
               ext_order: ExtinctionOrder):
    if name == 'Plant':
        return {'base_extinction_col': 'Plant', 'primary_extinction_col': 'Pollinator', 'd_col': 'd_poli_plant',
                'rate_col': "pollinator_r"}
    else:
        return {'base_extinction_col': 'Pollinator', 'primary_extinction_col': 'Plant', 'd_col': 'd_plant_poli',
                'rate_col': "plant_r"}


def get_plant_rates(is_polyploid: bool,
                    unknown_ext_rate: float,
                    polyploid_ext_rate: float,
                    diploid_ext_rate: float):
    if is_polyploid == -1:
        return unknown_ext_rate
    elif is_polyploid == 1:
        return polyploid_ext_rate
    elif is_polyploid == 0:
        return diploid_ext_rate
    else:
        raise Exception("Polyploidity should be -1,1 or 0")


def iteration_metrics(ep_network_edgelist, primary_extinct, base_extinct, config, degree, baseline):
    res = {
        'base_extinct_level': config['base_extinction_col'],
        'total_pollinator_in_base_extinctions': len(
            ep_network_edgelist[ep_network_edgelist['Pollinator'].isin(base_extinct)]['Pollinator'].unique()),
        'total_plants_in_base_extinctions': len(
            ep_network_edgelist[ep_network_edgelist['Plant'].isin(base_extinct)]['Plant'].unique()),
        'total_plants_in_base_extinctions_diploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid == 0][ep_network_edgelist['Plant'].isin(base_extinct)][
                'Plant'].unique()),
        'total_plants_in_base_extinctions_polyploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid == 1][ep_network_edgelist['Plant'].isin(base_extinct)][
                'Plant'].unique()),

        'total_plants_in_base_extinctions_unknown': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid == -1][
                ep_network_edgelist['Plant'].isin(base_extinct)]['Plant'].unique()),

        'connectence': ep_network_edgelist.shape[0] / ((len(ep_network_edgelist['Plant'].unique()) * len(
            ep_network_edgelist['Pollinator'].unique())) + 0.001),
        'total_pollinator_in_primary_extinctions': len(
            ep_network_edgelist[ep_network_edgelist['Pollinator'].isin(primary_extinct)]['Pollinator'].unique()),
        'total_plants_in_primary_extinctions': len(
            ep_network_edgelist[ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),
        'total_plants_in_primary_extinctions_diploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid == 0][
                ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),
        'total_plants_in_primary_extinctions_polyploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid == 1][
                ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),
        'total_plants_in_primary_extinctions_unknown': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid == -1][
                ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),

        'total_pollinators_prior_extinction': len(
            ep_network_edgelist['Pollinator'].unique()),
        'total_plants_prior_extinction': len(
            ep_network_edgelist['Plant'].unique()),
        'total_plants_prior_extinction_diploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid == 0]['Plant'].unique()),
        'total_plants_prior_extinction_polyploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid == 1]['Plant'].unique()),
        'total_plants_prior_extinction_unknown': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid == -1]['Plant'].unique()),

    }
    res.update({'degree': degree})
    res.update(baseline)
    return res


def extinction_analysis(classified_network: pd.DataFrame,
                        polyploid_range: tuple,
                        diploid_range: tuple,
                        unknown_range: tuple,
                        pollinator_range: tuple,
                        extinction_order: ExtinctionOrder = ExtinctionOrder.random):
    polyploid_ext_rate = np.random.uniform(low=polyploid_range[0], high=polyploid_range[1])
    diploid_ext_rate = np.random.uniform(low=diploid_range[0], high=diploid_range[1])
    unknown_ext_rate = np.random.uniform(low=unknown_range[0], high=unknown_range[1])
    pollinator_ext_rate = np.random.uniform(low=pollinator_range[0], high=pollinator_range[1])
    classified_network["pollinator_r"] = pollinator_ext_rate
    classified_network["plant_r"] = classified_network["is_polyploid"].apply(
        lambda is_poly: get_plant_rates(is_poly, unknown_ext_rate, polyploid_ext_rate, diploid_ext_rate))
    classified_network = update_edgelist_totals(classified_network)

    Base = 'Plant'
    Primary = 'Pollinator'
    config = get_config(Base)
    base_extinct = random.sample(list(classified_network[config['base_extinction_col']].unique()), 1)
    # print(base_extinct)
    iters_df = pd.DataFrame()
    degree = 0
    # print(ep_network_edgelist[['Plant', 'Pollinator']])
    baseline = {'baseline_n_pollinators': len(classified_network['Pollinator'].unique()),
                'baseline_n_plant': len(classified_network['Plant'].unique()),
                'baseline_n_polyploids': len(classified_network.query("is_polyploid == 1")['Plant'].unique()),
                'baseline_n_diploids': len(classified_network.query("is_polyploid == 0")['Plant'].unique()),
                'baseline_n_unknowns': len(classified_network.query("is_polyploid == -1")['Plant'].unique()),
                'baseline_connectence': len(classified_network) / ((len(classified_network['Plant'].unique()) * len(
                    classified_network['Pollinator'].unique())) + 0.001)}
    while True:
        # print(f"degree: {degree}")
        at_extinction_risk = \
            classified_network.loc[classified_network[config['base_extinction_col']].isin(base_extinct)][
                ['is_polyploid', config['d_col'], config['base_extinction_col'], config['primary_extinction_col'],
                 config['rate_col']]]
        if len(at_extinction_risk.index) == 0:
            primary_extinct = []
            break
        prob = np.random.uniform()
        at_extinction_risk['primary_extinct_prob'] = at_extinction_risk.apply(lambda row:
                                                                              row[config['d_col']] * row[
                                                                                  config['rate_col']], axis=1)
        primary_extinct = at_extinction_risk[at_extinction_risk.primary_extinct_prob >= prob][
            config['primary_extinction_col']].unique()
        if len(primary_extinct) == 0:
            primary_extinct = []
            break
        iter_metrics = iteration_metrics(classified_network, primary_extinct, base_extinct, config, degree, baseline)

        iters_df = iters_df.append(iter_metrics, ignore_index=True)
        # print(f'Primary extinct:\n{primary_extinct}')
        classified_network = classified_network.loc[
            ~classified_network[config['base_extinction_col']].isin(base_extinct)]  # Removing the base extinction
        # print( ep_network_edgelist[['Plant','Pollinator']])
        classified_network = update_edgelist_totals(classified_network)
        # print(iter_metrics)

        base_extinct = primary_extinct

        pass
        degree += 1
        Base, Primary = Primary, Base  # update
        config = get_config(Base)

    iter_metrics = iteration_metrics(classified_network, primary_extinct, base_extinct, config, degree, baseline)

    iters_df = iters_df.append(iter_metrics, ignore_index=True)
    # print(iters_df)
    return iters_df

# def main():
#     # working_dir = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/extinction_analysis/"
#     # weighted_features_data = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/features/plant_features/weighted/features_with_classification.csv"
#     # weighted_networks_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/networks/all/weighted/"
#     # binary_features_data = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/features/plant_features/binary/features_with_classification.csv"
#     # binary_networks_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/networks/all/binary/"
#     # binarized_features_data = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/features/plant_features/binarized_weighted/features_with_classification.csv"
#     # binarized_networks_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/networks/all/binarized_weighted/"
#     # networks = {'weighted': {'networks': weighted_networks_path, 'data': weighted_features_data},
#     #             'binary': {'networks': binary_networks_path, 'data': binary_features_data},
#     #             'binarized': {'networks': binarized_networks_path, 'data': binarized_features_data}
#     #             }
#     #
#     #
#     # for _type in ['binary']: #networks:
#     #     output_dir = f"/groups/itay_mayrose/halabikeren/plant_pollinator_networks/extinction_analysis/{_type}/"
#     #     os.makedirs(output_dir, exist_ok=True)
#     #     type_df = All_networks_extinction_analysis(networks[_type], output_dir=output_dir)
#     #     type_df_path = f"/groups/itay_mayrose/halabikeren/plant_pollinator_networks/extinction_analysis/all_{_type}_networks_simulation_data.csv"
#     #     type_df.to_csv(type_df_path)
#
#     network_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/networks/all/weighted/118.csv"
#     network_species_metedata_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/features/plant_features/weighted/features_with_classification.csv"
#     output_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/extinction_analysis/test.csv"
#     generate_network_extinction_simulations(network_path, network_species_metedata_path, output_path, iter=10)
#
# if __name__ == '__main__':
#     main()
