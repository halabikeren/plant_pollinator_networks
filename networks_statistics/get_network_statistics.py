
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

def extract_single_network_edgelist(network_path):
    print(network_path)
    ep_network = pd.read_csv(network_path, encoding= 'unicode_escape')
    ep_network.rename(columns={"plant/pollinator": "Plant"}, inplace=True)
    if "Plant" not in ep_network.columns:
        ep_network.rename(columns={ep_network.columns.tolist()[0]: "Plant"}, inplace=True)

    ep_network_edgelist = ep_network.melt(id_vars=["Plant"],
                                          var_name="Pollinator",
                                          value_name="Value")
    ep_network_edgelist["Value"] = pd.to_numeric(ep_network_edgelist['Value'], errors='coerce')
    ep_network_edgelist = ep_network_edgelist[ep_network_edgelist["Value"] > 0]



    return ep_network_edgelist


def extract_folder_networks(network_folder,source,network_type):
    networks_edgelist = pd.DataFrame()
    networks_files = [os.path.join(network_folder, f) for f in os.listdir(network_folder) if
                             os.path.isfile(os.path.join(network_folder, f))]
    for i,network_path in enumerate(networks_files):
        if not "reference" in network_path and network_path.endswith(".csv"):
            curr_network_edgelist = extract_single_network_edgelist(network_path)
            curr_network_edgelist["source"] = source
            curr_network_edgelist["network_type"] = network_type
            curr_network_edgelist["network_path"] = network_path
            networks_edgelist = pd.concat([networks_edgelist,curr_network_edgelist], ignore_index= True)
    return networks_edgelist


def main():
    networks_folder_path = "/groups/itay_mayrose/halabikeren/plant_pollinator_networks/networks/"
    name_resolution_file = os.path.join(networks_folder_path,'processed_resolved_plant_names.csv')
    name_resolution = pd.read_csv(name_resolution_file)
    networks_edgelists_path = "All_networks_edgelists.csv"
    if not os.path.exists(networks_edgelists_path):
        network_sources_folders = {o:os.path.join(networks_folder_path, o) for o in os.listdir(networks_folder_path) if os.path.isdir(os.path.join(networks_folder_path, o)) and "." not in o}
        print(network_sources_folders)
        all_networks_edgelists = pd.DataFrame()
        for source in network_sources_folders:
                binary_networks_folder = os.path.join(network_sources_folders[source],'binary')
                binary_network_edgelists = extract_folder_networks(binary_networks_folder,source,  network_type = "binary")
                weighted_networks_folder = os.path.join(network_sources_folders[source],'weighted')
                weighted_network_edgelists = extract_folder_networks(weighted_networks_folder,source,  network_type = "weighted")
                all_networks_edgelists =  pd.concat([all_networks_edgelists,binary_network_edgelists,weighted_network_edgelists], ignore_index= True)
        print(all_networks_edgelists)
        all_networks_edgelists.to_csv(networks_edgelists_path)
    else:
        all_networks_edgelists = pd.read_csv(networks_edgelists_path)

    all_networks_edgelists["Plant"] = all_networks_edgelists["Plant"].str.lower().str.replace("_"," ")

    all_networks_edgelists = pd.merge(all_networks_edgelists, name_resolution, how='left', left_on=['Plant'],
                                   right_on=['original_name'])
    all_networks_edgelists = all_networks_edgelists[["source","network_type","network_path","Plant", "matched_name", "genus", "family", "Pollinator", "Value"]]

    all_networks_edgelists = all_networks_edgelists.astype(
         {'matched_name': 'string', 'genus': 'string', 'family': 'string', 'Pollinator': 'string', 'Plant':'string','network_path':'string', 'source':'string','network_type': 'string'})
    all_networks_edgelists[["matched_name", "genus", "family", "Pollinator","Plant","network_path","network_type","source"]] = all_networks_edgelists[
         ["matched_name", "genus", "family", "Pollinator","Plant","network_path","network_type","source"]].fillna('unknown')

    all_networks_edgelists.to_csv(f"{os.getcwd()}/All_networks_edgelists_named.csv")
    #unique_plants = (all_networks_edgelists[["Plant","matched_name"]]).drop_duplicates()
    #print(f"Percentage of known plants: {np.mean(all_networks_edgelists['matched_name']!='unknown')}")


        # network_level_d = {'number_of_plants': len(np.unique(ep_network_edgelist["matched_name"])), 'number_of_pollinators': len(np.unique(ep_network_edgelist["Pollinator"])),
        #                     'total_number_of_visits': (np.sum(ep_network_edgelist["Value"])),
        #                    'total_number_of_interactions': len(ep_network_edgelist.index),
        #                    'pct_unknown_names' : np.mean(ep_network_edgelist["matched_name"]=="unknown") }
        #

if __name__ == '__main__':
    main()






