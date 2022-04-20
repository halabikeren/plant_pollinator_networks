import pandas as pd
import os


def main():
    rmangal_networks = pd.read_csv("all_rmangal_networks_new.csv")
    mangal_networks_folder = os.path.join(os.getcwd(), "all_rmangal_networks")
    if not os.path.exists(mangal_networks_folder):
        os.mkdir(mangal_networks_folder)
    mangal_weighted_networks = os.path.join(mangal_networks_folder, "all_rmangal_weighted")
    os.mkdir(mangal_weighted_networks)
    mangal_binary_networks = os.path.join(mangal_networks_folder, "all_rmangal_binary")
    os.mkdir(mangal_binary_networks)

    for network_ind in rmangal_networks["network_id"].unique():
        network_data = rmangal_networks[rmangal_networks["network_id"] == network_ind]
        if (network_data["total_interactions"]).isin([0, 1]).all():
            folder = mangal_binary_networks
        else:
            folder = mangal_weighted_networks
        pivoted_network = network_data.pivot(index="node_to_name", columns="node_from_name",
                                             values="total_interactions")
        pivoted_network.to_csv(os.path.join(folder, f'{network_ind}.csv'))

    rmangal_networks[["network_id", "doi"]].drop_duplicates().to_csv('mangal_refs.csv')


if __name__ == "__main__":
    main()
