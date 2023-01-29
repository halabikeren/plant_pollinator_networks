import pandas as pd
# import networkx as nx
import os
import pickle
import numpy as np
import random
import matplotlib.pyplot as plt
# from networkx.algorithms import bipartite
# from networkx_analysis.constants import *



def generate_network_extinction_simulations(network_path, network_species_metedata_path, output_path,  polyploid_range = [0.9,0.9], diploid_range = [0.9,0.9], unknown_range = [0.9,0.9], pollinator_range =[0.9,0.9], iter  = 5):
    ep_network = pd.read_csv(network_path, encoding='unicode_escape')
    network_species_metedata = pd.read_csv(network_species_metedata_path).query(f"network == '{os.path.basename(network_path)}'")
    ep_network.rename(columns={"plant/pollinator": "Plant"}, inplace=True)
    if "Plant" not in ep_network.columns:
        ep_network.rename(columns={ep_network.columns.tolist()[0]: "Plant"}, inplace=True)

    ep_network_edgelist = ep_network.melt(id_vars=["Plant"],
                                          var_name="Pollinator",
                                          value_name="Value")
    ep_network_edgelist['Plant'] = ep_network_edgelist['Plant'].str.lower()  # lowercase
    ep_network_edgelist['Pollinator'] = ep_network_edgelist['Pollinator'].str.lower()  # lowercase
    ep_network_edgelist['Value'] = pd.to_numeric(ep_network_edgelist['Value'], errors='coerce')
    ep_network_edgelist = ep_network_edgelist[ep_network_edgelist["Value"] > 0]
    ep_network_edgelist = ep_network_edgelist.merge(network_species_metedata, how='left', left_on='Plant',
                                                    right_on='original_name')
    ep_network_edgelist['is_polyploid'] = ep_network_edgelist['is_polyploid'].fillna(-1)
    extinction_dfs = []
    for i in range(iter):
        extinction_df = extinction_analysis(ep_network_edgelist, pollinator_range= pollinator_range, polyploid_range = polyploid_range, diploid_range = diploid_range, unknown_range= unknown_range)
        extinction_df["simluaiotn_index"] = i
        extinction_dfs.append(extinction_df)
        if (len(extinction_df.index)>2):
            print(f"succeeded in simulating more than 2 iterations at extinction simulation {i}")
    pd.concat(extinction_dfs).to_csv(output_path)

    print(f"failed to simulate extinction after {iter} attemps")




def All_networks_extinction_analysis(networks_data, output_dir, range_of_rates: list[int], num_sim: int = 100):
    network_folder = networks_data['networks']
    species_data_path = networks_data['data']
    networks_files = {f: os.path.join(network_folder, f) for f in os.listdir(network_folder) if
                      os.path.isfile(os.path.join(network_folder, f))}
    for i,network in enumerate(networks_files):
        network_path = networks_files[network]
        if not "reference" in network_path and network_path.endswith(".csv"):
             generate_network_extinction_simulations(network_path,
                                                     network_species_metedata_path=species_data_path,
                                                     polyploid_range=range_of_rates,
                                                     diploid_range=range_of_rates,
                                                     unknown_range=range_of_rates,
                                                     pollinator_range=range_of_rates,
                                                     iter=num_sim,
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


def get_config(name):
    if name == 'Plant':
        return {'base_extinction_col': 'Plant', 'primary_extinction_col': 'Pollinator', 'd_col': 'd_poli_plant',
                'rate_col': "pollinator_r"}
    else:
        return {'base_extinction_col': 'Pollinator', 'primary_extinction_col': 'Plant', 'd_col': 'd_plant_poli',
                'rate_col': "plant_r"}


def get_plant_rates(is_polyploid, unknown_ext_rate, polyploid_ext_rate, diploid_ext_rate):
    if is_polyploid == -1:
        return unknown_ext_rate
    elif is_polyploid == 1:
        return polyploid_ext_rate
    elif is_polyploid == 0:
        return diploid_ext_rate
    else:
        raise Exception("Polyploidity should be -1,1 or 0")


def iteration_metrics(ep_network_edgelist, primary_extinct, base_extinct, config, degree, baseline):
    res =   {
        'base_extinct_level': config['base_extinction_col'],
        'total_pollinator_in_base_extinctions':len(ep_network_edgelist[ep_network_edgelist['Pollinator'].isin(base_extinct)]['Pollinator'].unique()),
        'total_plants_in_base_extinctions':len(ep_network_edgelist[ep_network_edgelist['Plant'].isin(base_extinct)]['Plant'].unique()),
        'total_plants_in_base_extinctions_diploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==0][ep_network_edgelist['Plant'].isin(base_extinct)]['Plant'].unique()),
        'total_plants_in_base_extinctions_polyploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==1][ep_network_edgelist['Plant'].isin(base_extinct)]['Plant'].unique()),

        'total_plants_in_base_extinctions_unknown': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==-1][ep_network_edgelist['Plant'].isin(base_extinct)]['Plant'].unique()),

        'connectence': ep_network_edgelist.shape[0]/((len(ep_network_edgelist['Plant'].unique())*len(ep_network_edgelist['Pollinator'].unique()))+0.001),
        'total_pollinator_in_primary_extinctions':len(ep_network_edgelist[ep_network_edgelist['Pollinator'].isin(primary_extinct)]['Pollinator'].unique()),
        'total_plants_in_primary_extinctions':len(ep_network_edgelist[ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),
        'total_plants_in_primary_extinctions_diploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==0][ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),
        'total_plants_in_primary_extinctions_polyploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==1][ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),
        'total_plants_in_primary_extinctions_unknown': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==-1][ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),

        'total_pollinators_prior_extinction': len(
            ep_network_edgelist['Pollinator'].unique()),
        'total_plants_prior_extinction': len(
           ep_network_edgelist['Plant'].unique()),
        'total_plants_prior_extinction_diploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==0]['Plant'].unique()),
        'total_plants_prior_extinction_polyploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==1]['Plant'].unique()),
        'total_plants_prior_extinction_unknown': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==-1]['Plant'].unique()),


    }
    res.update({'degree': degree})
    res.update(baseline)
    return res



def extinction_analysis(ep_network_edgelist, polyploid_range, diploid_range, unknown_range, pollinator_range):
    polyploid_ext_rate = np.random.uniform(low=polyploid_range[0], high=polyploid_range[1])
    diploid_ext_rate = np.random.uniform(low=diploid_range[0], high=diploid_range[1])
    unknown_ext_rate = np.random.uniform(low=unknown_range[0], high=unknown_range[1])
    pollinator_ext_rate = np.random.uniform(low=pollinator_range[0], high=pollinator_range[1])
    ep_network_edgelist["pollinator_r"] = pollinator_ext_rate
    ep_network_edgelist["plant_r"] = ep_network_edgelist["is_polyploid"].apply(
        lambda x: get_plant_rates(x, unknown_ext_rate, polyploid_ext_rate, diploid_ext_rate))
    ep_network_edgelist = update_edgelist_totals(ep_network_edgelist)

    Base = 'Plant'
    Primary = 'Pollinator'
    config = get_config(Base)
    base_extinct = random.sample(list(ep_network_edgelist[config['base_extinction_col']].unique()), 1)
    #print(base_extinct)
    iters_df = pd.DataFrame()
    degree = 0
    #print(ep_network_edgelist[['Plant', 'Pollinator']])
    baseline = {'baseline_n_pollinators': len(ep_network_edgelist['Pollinator'].unique()),
                'baseline_n_plant': len(ep_network_edgelist['Plant'].unique()),
                'baseline_n_polyploids': len(ep_network_edgelist.query("is_polyploid == 1")['Plant'].unique()),
                'baseline_n_diploids': len(ep_network_edgelist.query("is_polyploid == 0")['Plant'].unique()),
                'baseline_n_unknowns': len(ep_network_edgelist.query("is_polyploid == -1")['Plant'].unique()),
                'baseline_connectence': len(ep_network_edgelist)/((len(ep_network_edgelist['Plant'].unique())*len(ep_network_edgelist['Pollinator'].unique()))+0.001)}
    while True:
        #print(f"degree: {degree}")
        at_extinction_risk = \
        ep_network_edgelist.loc[ep_network_edgelist[config['base_extinction_col']].isin(base_extinct)][
            ['is_polyploid', config['d_col'], config['base_extinction_col'], config['primary_extinction_col'], config['rate_col']]]
        if len(at_extinction_risk.index) == 0:
            primary_extinct = []
            break
        prob = np.random.uniform()
        at_extinction_risk['primary_extinct_prob'] = at_extinction_risk.apply(lambda row:
                                                                              row[config['d_col']] * row[
                                                                                  config['rate_col']], axis=1)
        primary_extinct = at_extinction_risk[at_extinction_risk.primary_extinct_prob >= prob][config['primary_extinction_col']].unique()
        if len(primary_extinct) == 0:
            primary_extinct = []
            break
        iter_metrics = iteration_metrics(ep_network_edgelist, primary_extinct, base_extinct, config,degree, baseline)

        iters_df = iters_df.append(iter_metrics, ignore_index=True)
        #print(f'Primary extinct:\n{primary_extinct}')
        ep_network_edgelist = ep_network_edgelist.loc[
            ~ep_network_edgelist[config['base_extinction_col']].isin(base_extinct)]  # Removing the base extinction
        #print( ep_network_edgelist[['Plant','Pollinator']])
        ep_network_edgelist = update_edgelist_totals(ep_network_edgelist)
        #print(iter_metrics)

        base_extinct = primary_extinct

        pass
        degree += 1
        Base, Primary = Primary, Base  # update
        config = get_config(Base)


    iter_metrics = iteration_metrics(ep_network_edgelist, primary_extinct, base_extinct, config, degree, baseline)

    iters_df = iters_df.append(iter_metrics, ignore_index=True)
    # print(iters_df)
    return iters_df

def main():
    # working_dir = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/extinction_analysis/"
    # weighted_features_data = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/features/plant_features/weighted/features_with_classification.csv"
    # weighted_networks_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/networks/all/weighted/"
    # binary_features_data = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/features/plant_features/binary/features_with_classification.csv"
    # binary_networks_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/networks/all/binary/"
    # binarized_features_data = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/features/plant_features/binarized_weighted/features_with_classification.csv"
    # binarized_networks_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/networks/all/binarized_weighted/"
    # networks = {'weighted': {'networks': weighted_networks_path, 'data': weighted_features_data},
    #             'binary': {'networks': binary_networks_path, 'data': binary_features_data},
    #             'binarized': {'networks': binarized_networks_path, 'data': binarized_features_data}
    #             }
    #
    #
    # for _type in ['binary']: #networks:
    #     output_dir = f"/groups/itay_mayrose/halabikeren/plant_pollinator_networks/extinction_analysis/{_type}/"
    #     os.makedirs(output_dir, exist_ok=True)
    #     type_df = All_networks_extinction_analysis(networks[_type], output_dir=output_dir)
    #     type_df_path = f"/groups/itay_mayrose/halabikeren/plant_pollinator_networks/extinction_analysis/all_{_type}_networks_simulation_data.csv"
    #     type_df.to_csv(type_df_path)

    network_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/networks/all/weighted/118.csv"
    network_species_metedata_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/features/plant_features/weighted/features_with_classification.csv"
    output_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/extinction_analysis/test.csv"
    generate_network_extinction_simulations(network_path, network_species_metedata_path, output_path, iter=10)

if __name__ == '__main__':
    main()
