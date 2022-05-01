import pandas as pd
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import re
import numpy as np
import os
import shutil



def author_name_wol(str):
        return re.findall(r'(\w+)\b',str)[0]



def author_name_stouffer(str):
        match = re.search(r' ([A-z]+) (&|et al|\()',str)
        if match:
            return match.group(1)
        else:
            return 'not found'


def year_regex(str):
    matches = re.findall(r"\b((?:19|20)\d{2})\b",str) #r"[()A-z .&,]+\d{4}"
    if len(matches)>0:
        return matches[0]
    else:
        return 'unknown'



def main():
    wol_weighted_references_path = '/Users/noa/Workspace/networks/web_of_life/weighted/references.csv'
    wol_weighted_references = pd.read_csv(wol_weighted_references_path)[["Reference","ID"]]
    wol_weighted_references["source"] = "weighted"
    wol_binary_references_path = '/Users/noa/Workspace/networks/web_of_life/binary/references.csv'
    wol_binary_references = pd.read_csv(wol_binary_references_path)[["Reference","ID"]].drop_duplicates()
    wol_binary_references["source"] = "binary"
    all_wol_references = pd.concat([wol_binary_references,wol_weighted_references]).drop_duplicates()
    stouffer_papers_path = 'raw_CSVs/stouffer_network_papers.csv'
    stouffer_papers = pd.read_csv(stouffer_papers_path)
    stouffer_papers["year"]=stouffer_papers["Reference"].apply(year_regex)
    stouffer_papers["names"]=stouffer_papers["Reference"].apply(author_name_stouffer)
    all_wol_references["year"]=all_wol_references["Reference"].apply(year_regex)
    all_wol_references["names"]=all_wol_references["Reference"].apply(author_name_wol)


    stouffer_papers.to_csv("stouffer_papers_new.csv")
    all_wol_references.to_csv("all_wol_papers_new.csv")

    res = pd.merge(stouffer_papers,all_wol_references,how = 'left', left_on=['year','names'],right_on=['year','names'])
    res.to_csv("stouffer_wol.csv")
    res["ID"].drop_duplicates().to_csv("relevant_wol_networks.csv")
    exracted_IDS = list(pd.unique(res["ID"]))
    extra_nwtowrks = ["M_PL_039","M_PL_036","M_PL_026","M_PL_045","M_PL_042","M_PL_051","M_PL_052"]
    all_stouffer_networks = exracted_IDS+extra_nwtowrks
    stouffer_networks_folder = "stouffer_networks"
    if not os.path.exists(stouffer_networks_folder):
        os.mkdir(stouffer_networks_folder)
    for network in all_stouffer_networks:
        network_weighted = f'/Users/noa/Workspace/networks/web_of_life/weighted/{network}.csv'
        network_binary = f'/Users/noa/Workspace/networks/web_of_life/binary/{network}.csv'
        if os.path.exists(network_weighted):
            shutil.copyfile(network_weighted, f'{stouffer_networks_folder}/{network}.csv')
        elif os.path.exists(network_binary):
            shutil.copyfile(network_binary, f'{stouffer_networks_folder}/{network}.csv')



if __name__ == '__main__':
    main()




