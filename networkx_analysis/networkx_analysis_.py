import pandas as pd
import networkx as nx
import os
import pickle
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite


def draw_network(G):
    node_colors = []
    for key in list(G.nodes.keys()):
        if G.nodes[key].get('is_polyploid') == 1 and G.nodes[key]['bipartite'] == 0:  # polyploids
            node_colors.append('b')
        elif G.nodes[key].get('is_polyploid') == 0 and G.nodes[key]['bipartite'] == 0:
            node_colors.append('g')
        elif G.nodes[key].get('is_polyploid') == -1 and G.nodes[key]['bipartite'] == 0:
            node_colors.append('y')
        else:
            node_colors.append('r')

    nx.draw_networkx(G,
                     pos=nx.kamada_kawai_layout(G, weight='Value'),
                     node_size=200,
                     font_size=5,
                     node_color=node_colors,
                     # width=[G.edges[e]['Value'] for e in G.edges],
                     with_labels=True)
    plt.show()



def enrich_nodes_with_graph_features(G, enriched_node_names):
        plant_nodes = {n for n, d in G.nodes(data=True) if d["bipartite"] == 0}
        pollinator_nodes = {n for n, d in G.nodes(data=True) if d["bipartite"] == 1}
        cc = bipartite.closeness_centrality(G, nodes=plant_nodes)
        enriched_node_names['cc'] = enriched_node_names['Node'].map(cc)
        dc = bipartite.degree_centrality(G, nodes=plant_nodes)
        enriched_node_names['dc'] = enriched_node_names['Node'].map(dc)
        bc = bipartite.betweenness_centrality(G, nodes=plant_nodes)
        enriched_node_names['dc'] = enriched_node_names['Node'].map(bc)
        clusr = bipartite.clustering(G, nodes=plant_nodes)
        enriched_node_names['clusr'] = enriched_node_names['Node'].map(clusr)
        redundancies = {}
        for node in G.nodes():
            try:
                redundancy = bipartite.node_redundancy(G, nodes=[node])
                redundancies.update(redundancy)
            except Exception:
                redundancy = {node: -1}
                redundancies.update(redundancy)
        enriched_node_names['redundancy'] = enriched_node_names['Node'].map(redundancies)
        P = bipartite.projected_graph(G, plant_nodes)
        ec = nx.eigenvector_centrality(P, max_iter=1000)
        enriched_node_names['pr_ec'] = enriched_node_names['Node'].map(ec)
        enriched_node_names['network_n_edges'] = G.size
        enriched_node_names['network_n_nodes'] = G.number_of_nodes()
        enriched_node_names['number_of_pollinators'] = len(pollinator_nodes)
        enriched_node_names['number_of_plants'] = len(plant_nodes)



def Generate_single_network(network_path, datasets_dict):
    trimmed_network_path = network_path.replace('/Users/noa/Workspace/plant_pollinator_networks_data','')
    network_features = datasets_dict['network_features'][datasets_dict['network_features']['orig_path']==trimmed_network_path]
    species_features = datasets_dict['species_features'][datasets_dict['species_features']['orig_path']==trimmed_network_path]
    species_features['species'] = species_features['species'].str.lower()
    ep_network = pd.read_csv(network_path, encoding='unicode_escape')
    ep_network.rename(columns={"plant/pollinator": "Plant"}, inplace=True)
    if "Plant" not in ep_network.columns:
        ep_network.rename(columns={ep_network.columns.tolist()[0]: "Plant"}, inplace=True)

    ep_network_edgelist = ep_network.melt(id_vars=["Plant"],
                                          var_name="Pollinator",
                                          value_name="Value")
    ep_network_edgelist['Plant'] = ep_network_edgelist['Plant'].str.lower()  # lowercase
    ep_network_edgelist['Pollinator'] = ep_network_edgelist['Pollinator'].str.lower()  # lowercase
    ep_network_edgelist["Value"] = pd.to_numeric(ep_network_edgelist['Value'], errors='coerce')
    ep_network_edgelist = ep_network_edgelist[ep_network_edgelist["Value"] > 0]
    G = nx.Graph()
    G.add_nodes_from(ep_network_edgelist['Plant'], bipartite=0)
    G_nodes_df = pd.DataFrame({'Node': G.nodes()})
    G_nodes_df['orig_path'] = trimmed_network_path
    G_nodes_df = G_nodes_df.merge(network_features,how='left', on = 'orig_path', suffixes=[None,"_redundant"])
    G_nodes_df = G_nodes_df.merge(species_features, how='left', left_on='Node', right_on='species',suffixes=[None,"_redundant"])
    G_nodes_df = G_nodes_df.merge(datasets_dict['ploidity_classifier'], how='left', left_on='Node', right_on='Species',suffixes=[None,"_redundant"])
    pct_missing = G_nodes_df['is_polyploid'].isnull().sum() / len(G_nodes_df)
    G_nodes_df['pct_missing_data'] = pct_missing
    if G_nodes_df['is_polyploid'].isnull().sum() / len(G_nodes_df) > 0.9:
        return None, pd.DataFrame()

    nodes_classification = {node_name: label for node_name, label in
                            zip(G_nodes_df['Species'],G_nodes_df['is_polyploid'].fillna(-1))}
    nx.set_node_attributes(G, values=nodes_classification, name='is_polyploid')
    G.add_nodes_from(ep_network_edgelist['Pollinator'], bipartite=1)
    G.add_weighted_edges_from(
        [(row['Plant'], row['Pollinator'], row['Value']) for idx, row in ep_network_edgelist.iterrows()],
        weight='Value')
    enrich_nodes_with_graph_features(G, G_nodes_df)
    # sb = bipartite.spectral_bipartivity(P, nodes=plant_nodes, weight='weight')
    # draw_network(G)
    return G, G_nodes_df


def extract_folder_networks(network_folder, datasets_dict):
    networks = []
    networks_files = [os.path.join(network_folder, f) for f in os.listdir(network_folder) if
                      os.path.isfile(os.path.join(network_folder, f))]
    for i, network_path in enumerate(networks_files):
        if not "reference" in network_path and network_path.endswith(".csv"):
            curr_network_obj, metadata = Generate_single_network(network_path, datasets_dict)
            metadata["network_ind"] = i
            if curr_network_obj is not None:
                networks.append({'network_obj': curr_network_obj, 'metadata': metadata})
    return networks


def extract_and_filter_networks(network_sources_folders, datasets_dict):
    relevant_networks = {}
    for source in network_sources_folders:
        if source != "all":
            binary_networks_folder = os.path.join(network_sources_folders[source], 'binary')
            binary_networks = extract_folder_networks(binary_networks_folder, datasets_dict)
            weighted_networks_folder = os.path.join(network_sources_folders[source], 'weighted')
            weighted_networks = extract_folder_networks(weighted_networks_folder, datasets_dict)
            relevant_networks[source] = {'binary': binary_networks, 'weighted': weighted_networks}
    return relevant_networks


def unify_metadata(relevant_networks):
    dfs = []
    for source in relevant_networks:
        binary_dfs = [network['metadata'] for network in relevant_networks[source]['binary']]
        if len(binary_dfs) > 0:
            binary_df = pd.concat(binary_dfs, sort=False)
            binary_df['network_type'] = 'binary'
        else:
            binary_df = pd.DataFrame()
        weighted_dfs = [network['metadata'] for network in relevant_networks[source]['weighted']]
        if len(weighted_dfs) > 0:
            weighted_df = pd.concat(weighted_dfs, sort=False)
            weighted_df['network_type'] = 'weighted'
        else:
            weighted_df = pd.DataFrame()
        all_source_networks = pd.concat([weighted_df, binary_df], sort=False)
        all_source_networks['network_source'] = source
        dfs.append(all_source_networks)
    return pd.concat(dfs, sort=False)


def generate_unified_features_datasets(network_metadata, weighted_species_features, binary_species_features, weighted_network_features, binary_network_features):
    network_metadata['orig_path'] = network_metadata['orig_path'].str.replace('/groups/itay_mayrose/halabikeren/plant_pollinator_networks','')
    network_metadata['network_numbered_file'] = network_metadata['network_index'].astype(str)+".csv"
    weighted_network_features['network_type'] = 'weighted'
    binary_network_features['network_type'] = 'binary'
    all_network_data = pd.concat([weighted_network_features,binary_network_features])
    network_features = network_metadata.merge(all_network_data,how = 'left', left_on= ['network_numbered_file', 'network_type'], right_on=['network','network_type'])
    weighted_species_features['network_type'] = 'weighted'
    binary_species_features['network_type'] = 'binary'
    all_species_data = pd.concat([weighted_species_features, binary_species_features])
    all_species_features = network_metadata.merge(all_species_data,how = 'left', left_on= ['network_numbered_file', 'network_type'], right_on=['network','network_type'])
    return network_features, all_species_features



def main():
    working_dir = '/Users/noa/Workspace/plant_pollinator_networks_data'
    binary_network_features = pd.read_csv(f'{working_dir}/features/network_features/binary/network_features.csv')
    weighted_network_features = pd.read_csv(f'{working_dir}/features/network_features/weighted/network_features.csv')
    network_metadata = pd.read_csv(f'{working_dir}/networks/networks_metadata.csv')
    classifier = pd.read_csv(f'{working_dir}/classification/classified_unresolved_plant_names.csv')
    weighted_species_features = pd.read_csv(f'{working_dir}/features/plant_features/weighted/plant_species_features.csv')
    binary_species_features = pd.read_csv(f'{working_dir}/features/plant_features/binary/plant_features.csv')
    network_features, species_features = generate_unified_features_datasets(network_metadata, weighted_species_features, binary_species_features,  #binary_species_features_df,
                                                                            weighted_network_features, binary_network_features)

    datasets_dict = {'network_features': network_features,'species_features': species_features,'ploidity_classifier': classifier}

    networks_folder_path = f"{working_dir}/networks"
    network_sources_folders = {o: os.path.join(networks_folder_path, o) for o in os.listdir(networks_folder_path) if
                               os.path.isdir(os.path.join(networks_folder_path, o)) and "." not in o}
    all_networks_path = 'all_networks'
    if os.path.exists(all_networks_path) and 1 == 2:
        with open(all_networks_path, 'rb') as ALL_NETWORKS:
            relevant_networks = pickle.load(ALL_NETWORKS)
    else:
        relevant_networks = extract_and_filter_networks(network_sources_folders, datasets_dict)
        with open(all_networks_path, 'wb') as ALL_NETWORKS:
            pickle.dump(relevant_networks, ALL_NETWORKS)
    print(relevant_networks)
    print(sum([len(relevant_networks[k]['binary']) for k in relevant_networks]))
    print(sum([len(relevant_networks[k]['weighted']) for k in relevant_networks]))
    all_metadata = unify_metadata(relevant_networks)
    all_metadata.to_csv('networks_metadata', index=False)
    print(all_metadata.head(10))


if __name__ == '__main__':
    main()
