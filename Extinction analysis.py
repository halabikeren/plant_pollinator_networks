import pandas as pd
import networkx as nx
import os
import pickle
import numpy as np
import random
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
from networkx_analysis.constants import *






def generate_networkx_graph(ep_network_edgelist,network_species_metedata):
    G = nx.Graph()
    G.add_nodes_from(ep_network_edgelist['Plant'], bipartite=0)

    nodes_classification = {node_name: label for node_name, label in
                            zip(network_species_metedata['original_name'],
                                network_species_metedata['is_polyploid'].fillna(-1))}
    nx.set_node_attributes(G, values=nodes_classification, name='polyploidity')
    G.add_nodes_from(ep_network_edgelist['Pollinator'], bipartite=1)
    G.add_weighted_edges_from(
        [(row['Plant'], row['Pollinator'], row['Value']) for idx, row in ep_network_edgelist.iterrows()],
        weight='Value')


def generate_network_extinction_simulations(network_path, network_species_metedata, n_iter = 1):
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
    ep_network_edgelist = ep_network_edgelist.merge(network_species_metedata, how = 'left', left_on = 'Plant',right_on='original_name')
    ep_network_edgelist['is_polyploid']=ep_network_edgelist['is_polyploid'].fillna(-1)
    generate_networkx_graph(ep_network_edgelist, network_species_metedata)
    simulation_results = pd.DataFrame()
    for i in range(n_iter):
        #simulation_results = simulation_results.append(extinction_analysis(ep_network_edgelist, min_r=0, max_r=0.3), ignore_index= True)
        #simulation_results = simulation_results.append(extinction_analysis(ep_network_edgelist, min_r=0.3, max_r=0.6),ignore_index= True)
        simulation_results = simulation_results.append(extinction_analysis(ep_network_edgelist, min_r=0.6, max_r=1),ignore_index= True)
    return simulation_results

    # sb = bipartite.spectral_bipartivity(P, nodes=plant_nodes, weight='weight')
    # draw_network(G)




def All_networks_extinction_analysis(networks_data):
    all_networks_simulation_data = []
    network_folder  = networks_data['networks']
    species_data=  networks_data['data']
    networks_files = {f : os.path.join(network_folder, f) for f in os.listdir(network_folder) if
                      os.path.isfile(os.path.join(network_folder, f))}
    for network in networks_files:
        network_path = networks_files[network]
        network_species_metedata = species_data.loc[species_data.network==network]
        if not "reference" in network_path and network_path.endswith(".csv"):
            curr_network_simulation_data= generate_network_extinction_simulations(network_path, network_species_metedata)
            if curr_network_simulation_data is not None:
                all_networks_simulation_data.append(curr_network_simulation_data)
    return pd.concat(all_networks_simulation_data)

def update_edgelist(ep_network_edgelist):
    ep_network_edgelist.loc[:,'total_per_plant'] = ep_network_edgelist.groupby(['Plant'])['Value'].transform(
        lambda x: x.sum())
    ep_network_edgelist.loc[:,'total_per_pollinator'] = ep_network_edgelist.groupby(['Pollinator'])['Value'].transform(
        lambda x: x.sum())
    ep_network_edgelist.loc[:,'d_plant_poli'] = ep_network_edgelist['Value'] / ep_network_edgelist['total_per_plant']
    ep_network_edgelist.loc[:,'d_poli_plant'] = ep_network_edgelist['Value'] / ep_network_edgelist['total_per_pollinator']
    return ep_network_edgelist

def switch(curr_config):
    if curr_config['base_extinction_col']=='Pollinator':
        return {'base_extinction_col': 'Plant','primary_extinction_col':'Pollinator','d_col':'d_poli_plant'}
    else:
        return {'base_extinction_col': 'Pollinator','primary_extinction_col':'Plant','d_col':'d_plant_poli'}


def extinction_analysis(ep_network_edgelist, min_r, max_r):
    r = np.random.uniform(low = min_r, high=max_r)
    ep_network_edgelist = update_edgelist(ep_network_edgelist)
    config = {'base_extinction_col': 'Plant','primary_extinction_col':'Pollinator','d_col':'d_poli_plant'}
    base_extinct = random.sample(list(ep_network_edgelist[config['base_extinction_col']].unique()), 1)
    degree=0
    total_extinctions = 0
    total_poli_extinctions = 0
    total_non_poli_extinctions = 0
    total_unknown_extinctions = 0
    while len(base_extinct)>0:
        degree+=1
        network = ep_network_edgelist[[config['d_col'],config['base_extinction_col'],config['primary_extinction_col']]]
        #print(f"network:\n{network}")
        #print(f"iter={degree}")
        #print(f"base_extinct:\n {base_extinct}")
        at_extinction_risk = ep_network_edgelist.loc[ep_network_edgelist[config['base_extinction_col']].isin(base_extinct)][['is_polyploid',config['d_col'],config['base_extinction_col'],config['primary_extinction_col']]]
        if len(at_extinction_risk)==0:
            break
        prob = np.random.uniform()
        at_extinction_risk['primary_extinct_prob'] = at_extinction_risk.apply(lambda row:
                                     row[config['d_col']]*r, axis=1)
        #print(f"At extinction risk:\n {at_extinction_risk[[config['primary_extinction_col'],'primary_extinct_prob']]}")
        primary_extinct = at_extinction_risk.loc[at_extinction_risk.primary_extinct_prob>=prob][[config['primary_extinction_col'],'is_polyploid']]
        total_extinctions = total_extinctions+ len(primary_extinct.index)
        if config['primary_extinction_col']=='Plant':
            total_poli_extinctions += len(primary_extinct.loc[primary_extinct['is_polyploid']==1])
            total_non_poli_extinctions += len(primary_extinct.loc[primary_extinct['is_polyploid'] == 0])
            total_unknown_extinctions += len(primary_extinct.loc[primary_extinct['is_polyploid'] == -1])
        #print(f"primary_extinct:\n {primary_extinct}")
        ep_network_edgelist = ep_network_edgelist.loc[~ep_network_edgelist[config['base_extinction_col']].isin(base_extinct)] # Removing the base extinction
        ep_network_edgelist = update_edgelist(ep_network_edgelist)
        base_extinct = primary_extinct[config['primary_extinction_col']] #update
        config = switch(config)
    res = {'degree': degree, 'total_extinctions': total_extinctions,
            'total_unknown_extinctions': total_unknown_extinctions,
            'total_poli_extinctions': total_poli_extinctions,
            'total_not_poli_extinctions': total_non_poli_extinctions,
            'r_range': f'{min_r}_{max_r}'}
    print(res)
    return res


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
    networks = {'weighted': {'networks':weighted_networks_path, 'data': weighted_features_data},
                'binary': {'networks':binary_networks_path, 'data': binary_features_data}
                }

    for type in networks:
        type_df = All_networks_extinction_analysis(networks[type])
        type_df_path=  f'all_{type}_networks_simulation_data.csv'
        type_df.to_csv(type_df_path)




if __name__ == '__main__':
    main()
