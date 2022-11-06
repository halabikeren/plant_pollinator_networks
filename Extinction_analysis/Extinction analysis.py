import pandas as pd
import networkx as nx
import os
import pickle
import numpy as np
import random
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
from networkx_analysis.constants import *



def generate_network_extinction_simulations(network_path, network_species_metedata,network_name,  polyploid_range = [0.9,0.9], diploid_range = [0.9,0.9], unknown_range = [0.9,0.9], pollinator_range =[0.9,0.9]):
    ep_network = pd.read_csv(network_path, encoding='unicode_escape')
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
    extinction_df = extinction_analysis(ep_network_edgelist, pollinator_range= pollinator_range, polyploid_range = polyploid_range, diploid_range = diploid_range, unknown_range= unknown_range)
    if(len(extinction_df.index)>2):
        extinction_df.to_csv(f"/Users/noa/Workspace/plant_pollinator_networks/Extinction_analysis/Extinction_results/{network_name}.csv")

    # sb = bipartite.spectral_bipartivity(P, nodes=plant_nodes, weight='weight')
    # draw_network(G)


def All_networks_extinction_analysis(networks_data):
    all_networks_simulation_data = []
    network_folder = networks_data['networks']
    species_data = networks_data['data']
    networks_files = {f: os.path.join(network_folder, f) for f in os.listdir(network_folder) if
                      os.path.isfile(os.path.join(network_folder, f))}
    for i,network in enumerate(networks_files):
        network_path = networks_files[network]
        network_species_metedata = species_data.loc[species_data.network == network]
        if not "reference" in network_path and network_path.endswith(".csv"):
            curr_network_simulation_data = generate_network_extinction_simulations(network_path,
                                                                                   network_species_metedata, network_name=i)
            if curr_network_simulation_data is not None:
                all_networks_simulation_data.append(curr_network_simulation_data)
    return pd.concat(all_networks_simulation_data)


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
        'total_pollinator_base_extinctions':len(ep_network_edgelist[ep_network_edgelist['Pollinator'].isin(base_extinct)]['Pollinator'].unique()),
        'total_plants_base_extinctions':len(ep_network_edgelist[ep_network_edgelist['Plant'].isin(base_extinct)]['Plant'].unique()),
        'total_plants_base_extinctions_diploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==0][ep_network_edgelist['Plant'].isin(base_extinct)]['Plant'].unique()),
        'total_plants_base_extinctions_polyploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==1][ep_network_edgelist['Plant'].isin(base_extinct)]['Plant'].unique()),

        'total_plants_base_extinctions_unknown': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==-1][ep_network_edgelist['Plant'].isin(base_extinct)]['Plant'].unique()),

        'connectence': len(ep_network_edgelist)/((len(ep_network_edgelist['Plant'].unique())*len(ep_network_edgelist['Pollinator'].unique()))+0.001),
        'total_pollinator_primary_extinctions':len(ep_network_edgelist[ep_network_edgelist['Pollinator'].isin(primary_extinct)]['Pollinator'].unique()),
        'total_plants_primary_extinctions':len(ep_network_edgelist[ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),
        'total_plants_primary_extinctions_diploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==0][ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),
        'total_plants_primary_extinctions_polyploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==1][ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),
        'total_plants_primary_extinctions_unknown': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==-1][ep_network_edgelist['Plant'].isin(primary_extinct)]['Plant'].unique()),

        'total_plants_prior_update': len(
           ep_network_edgelist['Plant'].unique()),
        'total_plants_prior_update_diploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==0]['Plant'].unique()),
        'total_plants_prior_update_polyploid': len(
            ep_network_edgelist[ep_network_edgelist.is_polyploid==1]['Plant'].unique()),
        'total_plants_prior_update_unknown': len(
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
    baseline = {'baseline_n_pollinators': len(ep_network_edgelist['Pollinator'].unique()), 'baseline_n_plant': len(ep_network_edgelist['Plant'].unique()), 'baseline_connectence': len(ep_network_edgelist)/((len(ep_network_edgelist['Plant'].unique())*len(ep_network_edgelist['Pollinator'].unique()))+0.001)}
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
    print(iters_df)
    return iters_df

def main():
    working_dir = '/Users/noa/Workspace/plant_pollinator_networks_data'
    weighted_features_data = pd.read_csv(
        f'{working_dir}/weighted_plants/plant_features_with_classification.csv')
    weighted_networks_path = f'{working_dir}/networks/all/weighted'
    binary_features_data = pd.read_csv(
        f'{working_dir}/binary_plants/plant_features_with_classification.csv')
    binary_networks_path = f'{working_dir}/networks/all/binary'
    binarized_features_data = pd.read_csv(
        f'{working_dir}/binarized_plants/plant_features_with_classification.csv')
    binarized_networks_path = f'{working_dir}/networks/all/binary'
    networks = {'weighted': {'networks': weighted_networks_path, 'data': weighted_features_data},
                'binary': {'networks': binary_networks_path, 'data': binary_features_data}
                }


    for type in networks:
        type_df = All_networks_extinction_analysis(networks[type])
        type_df_path = f'all_{type}_networks_simulation_data.csv'
        type_df.to_csv(type_df_path)


if __name__ == '__main__':
    main()
