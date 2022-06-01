import os
import re
from enum import Enum
from typing import Optional, Tuple
import logging
import numpy as np
import pandas as pd
from collections import Counter

logger = logging.getLogger(__name__)

class SelectionMethod(str, Enum):
    MAJORITY = "majority"
    BEST = "best"


class ResolvedNamesCurator:

    unresolved_names: pd.DataFrame
    resolved_names: list[pd.DataFrame]
    selection_method: SelectionMethod = SelectionMethod.BEST
    reference_resolved_names: Optional[pd.DataFrame]
    reference_unresolved_names: Optional[pd.DataFrame]
    work_dir: str = os.getcwd()

    def __init__(self,
                 unresolved_names_path: str,
                 resolved_names_paths: list[str],
                 selection_method: SelectionMethod,
                 reference_resolved_names: Optional[pd.DataFrame] = None,
                 reference_unresolved_names: Optional[pd.DataFrame] = None,
                 work_dir: Optional[str] = None,
                 is_processed: bool = False):
        logger.info(f"parsing unresolved names")
        self.unresolved_names = pd.read_csv(unresolved_names_path)
        logger.info(f"parsing {len(resolved_names_paths)} resolved names alternatives")
        if not is_processed:
            self.resolved_names = [self.get_resolved_names(path) for path in resolved_names_paths]
        else:
            self.resolved_names = [pd.read_csv(path) for path in resolved_names_paths]
        self.selection_method = selection_method
        self.reference_resolved_names = reference_resolved_names
        self.reference_unresolved_names = reference_unresolved_names
        if work_dir is not None:
            self.work_dir = work_dir

    @staticmethod
    def parse_name_mapping(resolved_names: pd.DataFrame):
        resolved_names.drop([col for col in resolved_names.columns if col.startswith("Unnamed")], inplace=True, axis=1)
        resolved_names.dropna(subset=["original_name", "resolved_name"], how="any", inplace=True)
        resolved_names["original_name"] = resolved_names["original_name"].str.lower()
        resolved_names["matched_name"] = resolved_names["matched_name"].str.lower()
        resolved_names["resolved_name"] = resolved_names["resolved_name"].str.lower()
        if "Old.status" in resolved_names.columns: # historic correction
            resolved_names.drop("orig_status", axis=1, inplace=True)
            resolved_names.rename(columns={"Old.status": "orig_status"}, inplace=True)

    @staticmethod
    def map_synonym_to_accepted(resolved_names: pd.DataFrame):
        if "name_string" in resolved_names.columns:
            resolved_names["resolved_name"].fillna(resolved_names["name_string"].to_dict(), inplace=True)
            name_string_to_matched_name = \
            resolved_names.loc[resolved_names.orig_status == "ACCEPTED"].set_index("name_string")[
                "resolved_name"].to_dict()
            resolved_names["resolved_name"] = resolved_names[["name_string", "resolved_name"]].apply(
                lambda row: name_string_to_matched_name[
                    row.matched_name] if row.matched_name in name_string_to_matched_name else row.matched_name, axis=1)
        resolved_names.matched_name = resolved_names.matched_name.str.lower()

    @staticmethod
    def parse_taxonomic_data(resolved_names: pd.DataFrame):
        if "taxon_rank" in resolved_names.columns:
            resolved_names.taxon_rank = resolved_names.taxon_rank.apply(lambda s: s.lower() if type(s) == str else s)
            resolved_names.taxon_rank = resolved_names.taxon_rank.replace(
                {'spec.': 'species', 'gen.': 'genus', 'f.': 'form', 'var.': 'variety', 'subsp.': 'subspecies'})
            if "classification_path_ranks" in resolved_names.columns:
                resolved_names.classification_path_ranks = resolved_names.classification_path_ranks.apply(
                    lambda x: x.replace('spec.', 'species').replace('gen.', 'genus').replace('f.', 'form').replace('var.',
                                                                                                                   'variety') if pd.notna(
                        x) else np.nan)
                # collapse taxonomic levels lower than species to a species level, namely: 'subspecies', 'variety', 'subvariety', 'form'
                resolved_names.loc[
                    resolved_names.taxon_rank.isin(['subspecies', 'variety', 'subvariety', 'form']), "matched_name"] = \
                    resolved_names.loc[
                        resolved_names.taxon_rank.isin(['subspecies', 'variety', 'subvariety', 'form'])].apply(
                        lambda row: ResolvedNamesCurator.get_taxonomic_classification(row, rank="species"), axis=1)
                resolved_names.loc[resolved_names.taxon_rank.isin(
                    ['subspecies', 'variety', 'subvariety', 'form']), "taxon_rank"] = "species"
                resolved_names["genus"] = resolved_names.apply(
                    lambda row: ResolvedNamesCurator.get_taxonomic_classification(row, rank="genus"),
                    axis=1)
                resolved_names["family"] = resolved_names.apply(
                    lambda row: ResolvedNamesCurator.get_taxonomic_classification(row, rank="family"), axis=1)

    @staticmethod
    def get_taxonomic_classification(row: pd.Series, rank: str):
        if pd.isna(row.classification_path_ranks) or pd.isna(row.classification_path):
            return np.nan
        ordered_ranks = row.classification_path_ranks.lower().split("|")
        ordered_names = row.classification_path.lower().split("|")
        rank_index, rank_value = np.nan, np.nan
        try:
            rank_index = ordered_ranks.index(rank)
            rank_value = ordered_names[rank_index]
        except Exception as e:
            if rank == 'species':
                rank_value_components = ordered_names[-1].split(" ")
                rank_value = " ".join(rank_value_components[:2])
            else:
                pass
        return rank_value

    def get_resolved_names(self, path: str) -> pd.DataFrame:
        processed_data_path = f"{self.work_dir}/processed_{os.path.basename(path)}"
        if not os.path.exists(processed_data_path):
            resolved_names = pd.read_csv(path)
            ResolvedNamesCurator.parse_name_mapping(resolved_names=resolved_names)
            ResolvedNamesCurator.map_synonym_to_accepted(resolved_names=resolved_names)
            ResolvedNamesCurator.parse_taxonomic_data(resolved_names=resolved_names)
            resolved_names.drop_duplicates(["original_name", "resolved_name"], inplace=True)
            resolved_names.to_csv(processed_data_path, index=False)
        else:
            resolved_names = pd.read_csv(processed_data_path)
        return resolved_names

    @staticmethod
    def is_contradicted(record: pd.Series) -> bool:
        matched_name_fields = [col for col in record.index if "_matched_name" in col]
        matched_names = record[matched_name_fields].unique()
        return len(matched_names) > 1

    def get_majority_names(self):
        method_to_matches = [self.resolved_names[i][["original_name", "matched_name", "orig_status", "genus", "family"]].rename(
            columns={"matched_name": f"{i}_matched_name",
                     "orig_status": f"{i}_orig_status",
                     "genus": f"{i}_genus",
                     "family": f"{i}_family"}) for i in range(len(self.resolved_names))]
        methods_to_consider = list(method_to_matches.keys())
        joined_data = method_to_matches[methods_to_consider[0]]
        for method in methods_to_consider[1:]:
            joined_data = joined_data.merge(method_to_matches[method], on="original_name", how="outer")
        joined_data = joined_data.loc[:, ~joined_data.columns.duplicated()]
        logger.info(
            f"# records in joined data before removing records without matched names = {joined_data.shape[0]:,}")
        joined_data.dropna(subset=[col for col in joined_data.columns if col.endswith("resolved_name")], how="all",
                           axis=0, inplace=True)
        joined_data = joined_data.loc[
            joined_data["original_name"].str.lower().isin((self.unresolved_names['species_name'].str.lower().unique()))]
        logger.info(
            f"# records in joined data after removing records without resolved names = {joined_data.shape[0]:,}")
        joined_data["is_contradicted"] = joined_data.apply(ResolvedNamesCurator.is_contradicted, axis=1)
        logger.info(f"# contradicted records = {joined_data.loc[joined_data.is_contradicted].shape[0]:,}")

        joined_data[["resolved_name", "matched_name", "genus", "family", "taxon_rank",
                     "majority_support"]] = joined_data.apply(ResolvedNamesCurator.get_majority, axis=1,
                                                              result_type="expand")
        logger.info(f"majority computed")
        joined_data[["is_reported_as_synonym", "most_supported_synonym_resolved_name",
                     "most_supported_synonym_matched_name", "most_supported_synonym_genus",
                     "most_supported_synonym_family", "most_supported_synonym_rank"]] = joined_data.apply(
            lambda x: ResolvedNamesCurator.get_majority_accepted_name(x), axis=1, result_type="expand")
        logger.info(f"superior synonym computed")
        joined_data[["is_classified_at_lower_than_genus", "most_supported_non_genus_resolved_name",
                     "most_supported_non_genus_matched_name", "most_supported_non_genus_genus",
                     "most_supported_non_genus_family", "most_supported_non_genus_taxon_rank"]] = joined_data.apply(
            lambda x: ResolvedNamesCurator.get_lower_than_genus_rank(x), axis=1, result_type="expand")
        logger.info(f"superior taxonomic classification computed")
        reported_as_synonym = joined_data.loc[joined_data.is_reported_as_synonym & (
                joined_data.resolved_name != joined_data.most_supported_synonym_resolved_name)][
            ["original_name"] + [col for col in joined_data.columns if col.startswith("most_supported_synonym_")]]
        reported_as_synonym.rename(
            columns={col: col.replace("most_supported_synonym_", "") for col in reported_as_synonym.columns},
            inplace=True)
        reported_at_lower_than_genus = \
            joined_data.loc[joined_data.is_classified_at_lower_than_genus & (joined_data.taxon_rank == 'genus')][
                ["original_name"] + [col for col in joined_data.columns if col.startswith("most_supported_non_genus_")]]
        reported_at_lower_than_genus.rename(
            columns={col: col.replace("most_supported_non_genus_", "") for col in reported_at_lower_than_genus.columns},
            inplace=True)
        resolved_by_majority = joined_data.loc[(joined_data.majority_support >= 0.25) & (
            ~joined_data.original_name.isin(reported_as_synonym.original_name)) & (~joined_data.original_name.isin(
            reported_at_lower_than_genus.original_name))][
            ["original_name", "matched_name", "genus", "family", "taxon_rank"]]
        resolved_by_majority = pd.concat([resolved_by_majority, reported_as_synonym, reported_at_lower_than_genus])
        return resolved_by_majority

    def get_best_names(self) -> pd.DataFrame:
        scores = []
        if self.reference_resolved_names is not None:
            self.reference_resolved_names.original_name = self.reference_resolved_names.original_name.str.lower()
            self.reference_resolved_names = self.reference_resolved_names[["original_name", "resolved_name"]].rename(
                columns={"resolved_name": "prev_resolved_name"})

        for i in range(len(self.resolved_names)):
            score_measures = []
            resolved_names_df = self.resolved_names[i]
            total_coverage = resolved_names_df.shape[0] / self.unresolved_names.shape[0]
            score_measures.append(total_coverage)
            one_to_one_unresolved_coverage = resolved_names_df.drop_duplicates("resolved_name").shape[0] / self.unresolved_names.shape[0]
            score_measures.append(one_to_one_unresolved_coverage)
            if self.reference_unresolved_names is not None:
                self.reference_unresolved_names.species_name = self.reference_unresolved_names.species_name.str.lower()
                one_to_one_unresolved_coverage = resolved_names_df.loc[resolved_names_df["original_name"].isin(self.reference_unresolved_names['species_name'])].drop_duplicates("resolved_name").shape[0] / \
                                                 self.reference_unresolved_names.shape[0]
                score_measures.append(one_to_one_unresolved_coverage)
                resolved_names_df["one_to_one_unresolved_coverage"] = one_to_one_unresolved_coverage
            if self.reference_resolved_names is not None:
                ref_resolved_to_df_resolved =  resolved_names_df.loc[resolved_names_df["original_name"].isin(self.reference_resolved_names.original_name)][["original_name", "resolved_name"]].merge(self.reference_resolved_names, left_on="original_name", right_on="original_name", how="left")
                ref_resolved_to_df_resolved = ref_resolved_to_df_resolved[["resolved_name", "prev_resolved_name"]].drop_duplicates()
                one_to_one_resolved_coverage = ref_resolved_to_df_resolved.drop_duplicates("resolved_name").shape[0] / ref_resolved_to_df_resolved.shape[0]
                score_measures.append(one_to_one_resolved_coverage)
                resolved_names_df["one_to_one_resolved_coverage"] = one_to_one_resolved_coverage
                self.resolved_names[i] = resolved_names_df
            resolved_names_df.to_csv(f"{self.work_dir}/resolved_names_df_{i}_with_scores.csv", index=False)
            score = np.mean(score_measures)
            scores.append(score)
        self.resolved_names = [df for _, df in sorted(zip(scores, self.resolved_names), reverse=True)]
        best_resolved_names = self.resolved_names[0]
        # complement missing data by the next best methods
        missing_data_fraction = best_resolved_names.shape[0] / self.unresolved_names.shape[0]
        index = 1
        while missing_data_fraction > 0 and index < len(self.resolved_names):
            complement_resolved_names = self.resolved_names[index].loc[~self.resolved_names[index]["original_name"].isin(best_resolved_names["original_name"])]
            logger.info(f"complementing {complement_resolved_names.shape[0]} resolved names from the {index}th best dataframe")
            complement_resolved_names = complement_resolved_names[[col for col in complement_resolved_names.columns if col in best_resolved_names.columns]]
            best_resolved_names = pd.concat([best_resolved_names, complement_resolved_names])
            index += 1
            missing_data_fraction = best_resolved_names.shape[0] / self.unresolved_names.shape[0]
        # for names mapped to genera, look for a species-wise mapping in one of the other resources
        names_mapped_to_genera = best_resolved_names.loc[~best_resolved_names.matched_name.str.contains(" "), "original_name"]
        index = 1
        while len(names_mapped_to_genera) > 0 and index < len(self.resolved_names):
            logger.info(f"current number of original names mapped to resolved names at genus level = {len(names_mapped_to_genera)}")
            complement_resolved_names = self.resolved_names[index]
            species_wise_complement_data = complement_resolved_names.loc[(complement_resolved_names["original_name"].isin(names_mapped_to_genera)) & (complement_resolved_names.matched_name.str.contains(" "))].set_index("original_name")["resolved_name"].to_dict()
            best_resolved_names.loc[~best_resolved_names.matched_name.str.contains(" ", na=False), "resolved_name"] = best_resolved_names.loc[~best_resolved_names.matched_name.str.contains(" "), "original_name"].apply(lambda q: species_wise_complement_data[q] if q in species_wise_complement_data else q)
            index += 1
            names_mapped_to_genera = best_resolved_names.loc[~best_resolved_names.matched_name.str.contains(" ", na=False), "original_name"]
        logger.info(f"# original names mapped to resolved names at genus level = { best_resolved_names.loc[~best_resolved_names.matched_name.str.contains(' ', na=False)].shape[0]}")
        logger.info(f"% coverage by the best resolved names df = {np.round(best_resolved_names.shape[0] / self.unresolved_names.shape[0] * 100, 2)}")
        return best_resolved_names

    @staticmethod
    def get_majority(record: pd.Series) -> Tuple[str, str, str, str, str, float]:
        resolved_name_fields = [col for col in record.index if col.endswith("_resolved_name")]
        resolved_record_names = record[resolved_name_fields].dropna().tolist()
        resolved_name_to_support = Counter(resolved_record_names)
        for name in resolved_name_to_support:
            resolved_name_to_support[name] = resolved_name_to_support[name] / len(resolved_record_names)
        if len(matched_names) > 0:
            majority_resolved_name = max(resolved_name_to_support, key=resolved_name_to_support.get)
            majority_support = np.round(resolved_name_to_support[majority_resolved_name], 2)
            winning_method = [col.replace("_resolved_name", "") for col in record.index if
                              col.endswith("_resolved_name") and record[col] == majority_resolved_name][0]
            majority_matched_name = record[f"{winning_method}_matched_name"]
            majority_genus = record[f"{winning_method}_genus"]
            majority_family = record[f"{winning_method}_family"]
            majority_taxon_rank = record[f"{winning_method}_taxon_rank"]
        else:
            majority_resolved_name = np.nan
            majority_matched_name = np.nan
            majority_genus = np.nan
            majority_family = np.nan
            majority_taxon_rank = np.nan
            majority_support = np.nan
        return majority_resolved_name, majority_matched_name, majority_genus, majority_family, majority_taxon_rank, majority_support

    @staticmethod
    def get_majority_accepted_name(record: pd.Series) -> Tuple[bool, str, str, str, str, str]:
        method_to_status = {col.replace("_orig_status", ""): record[col] for col in record.index if
                            col.endswith("orig_status")}
        is_synonym = np.any([status == "SYNONYM" for status in method_to_status.values()])
        most_frequent_synonym_resolved_name, most_frequent_synonym_matched_name, most_frequent_synonym_genus, most_frequent_synonym_family, most_frequent_synonym_taxon_rank = np.nan, np.nan, np.nan, np.nan, np.nan
        if is_synonym:
            methods_reporting_synonym = [method for method in method_to_status if method_to_status[method] == "SYNONYM"]
            synonyms = [record[f"{method}_resolved_name"] for method in methods_reporting_synonym]
            most_frequent_synonym_resolved_name = max(synonyms, key=synonyms.count)
            most_frequent_synonym_matched_name = \
            [record[f"{method}_matched_name"] for method in methods_reporting_synonym if
             f"{method}_matched_name" in record.index and record[
                 f"{method}_resolved_name"] == most_frequent_synonym_resolved_name][0]
            most_frequent_synonym_genus = [record[f"{method}_genus"] for method in methods_reporting_synonym if
                                           f"{method}_genus" in record.index and record[
                                               f"{method}_resolved_name"] == most_frequent_synonym_resolved_name][0]
            most_frequent_synonym_family = [record[f"{method}_family"] for method in methods_reporting_synonym if
                                            f"{method}_family" in record.index and record[
                                                f"{method}_resolved_name"] == most_frequent_synonym_resolved_name][0]
            most_frequent_synonym_taxon_rank = \
            [record[f"{method}_taxon_rank"] for method in methods_reporting_synonym if
             f"{method}_taxon_rank" in record.index and record[
                 f"{method}_resolved_name"] == most_frequent_synonym_resolved_name][0]
        return is_synonym, most_frequent_synonym_resolved_name, most_frequent_synonym_matched_name, most_frequent_synonym_genus, most_frequent_synonym_family, most_frequent_synonym_taxon_rank

    @staticmethod
    def get_lower_than_genus_rank(record: pd.Series) -> Tuple[bool, str, str, str, str, str]:
        relevent_record_data = record[[col for col in record.index if col.endswith('_taxon_rank')]].dropna().to_dict()
        method_to_rank = {col.replace("_taxon_rank", ""): relevent_record_data[col] for col in relevent_record_data}
        is_lower_than_genus = np.any([taxon_rank != "genus" for taxon_rank in method_to_rank.values()])
        methods_reporting_lower_than_genus = [method for method in method_to_rank if
                                              type(method_to_rank[method]) is str and method_to_rank[method] != "genus"]
        ranks = [method_to_rank[method] for method in methods_reporting_lower_than_genus]
        most_frequent_resolved_name, most_frequent_matched_name, most_frequent_genus, most_frequent_family, most_frequent_taxon_rank = np.nan, np.nan, np.nan, np.nan, np.nan
        if len(ranks) > 0:
            most_frequent_rank = max(ranks, key=ranks.count)
            resolved_names = record[[col for col in record.index if col.endswith('_resolved_name')]].dropna().to_dict()
            method_to_name = {col.replace("_resolved_name", ""): resolved_names[col] for col in resolved_names}
            most_frequent_resolved_name = \
            [method_to_name[method] for method in method_to_rank if method_to_rank[method] == most_frequent_rank][0]
            most_frequent_matched_name = [record[f"{method}_matched_name"] for method in method_to_rank if
                                          f"{method}_matched_name" in record.index and record[
                                              f"{method}_resolved_name"] == most_frequent_resolved_name][0]
            most_frequent_genus = [record[f"{method}_genus"] for method in method_to_rank if
                                   record[f"{method}_resolved_name"] == most_frequent_resolved_name][0]
            most_frequent_family = [record[f"{method}_family"] for method in method_to_rank if
                                    record[f"{method}_resolved_name"] == most_frequent_resolved_name][0]
            most_frequent_taxon_rank = [record[f"{method}_taxon_rank"] for method in method_to_rank if
                                        record[f"{method}_resolved_name"] == most_frequent_resolved_name][0]
        return is_lower_than_genus, most_frequent_resolved_name, most_frequent_matched_name, most_frequent_genus, most_frequent_family, most_frequent_taxon_rank

    def process_resolved_names(self, output_path: str):
        if self.selection_method == SelectionMethod.MAJORITY:
            resolved_names = self.get_majority_names()
        else:
            resolved_names = self.get_best_names()
        resolved_names = resolved_names[["original_name", "name_sources", "matched_name", "taxon_rank", "genus", "family"]].rename(columns={"original_name": "original_name", "matched_name": "resolved_name"})
        species_name_regex = re.compile("(\w*\s*\w*)")
        resolved_names["resolved_species_name"] = resolved_names.resolved_name.apply(lambda name: species_name_regex.search(name).group(1) if pd.notna(name) else name)
        resolved_names.to_csv(output_path, index=False)


if __name__ == '__main__':

    import sys
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s module: %(module)s function: %(funcName)s line %(lineno)d: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), ],
        force=True,  # run over root logger settings to enable simultaneous writing to both stdout and file handler
    )

    data_unresolved_names_path = "/groups/itay_mayrose/halabikeren/PloiDB/name_resolution/all_unresolved_names.csv"
    gnr_all_path = "/groups/itay_mayrose/halabikeren/PloiDB/name_resolution/resolved_names_different_methods/gnr_all_resolved_names_processed.csv"
    gnr_wfo_path = "/groups/itay_mayrose/halabikeren/PloiDB/name_resolution/resolved_names_different_methods/gnr_wfo_resolved_names_processed.csv"
    gnr_majority_path = "/groups/itay_mayrose/halabikeren/PloiDB/name_resolution/resolved_names_different_methods/resolved_names_by_majority_rule.csv"
    data_resolved_names_paths = [gnr_all_path, gnr_wfo_path, gnr_majority_path]
    data_selection_method = SelectionMethod.BEST
    reference_resolved_names_path = "/groups/itay_mayrose/halabikeren/PloiDB/ccdb/all_data.csv"
    reference_unresolved_names_path = "/groups/itay_mayrose/halabikeren/PloiDB/name_resolution/trees_unresolved_names.csv"
    data_output_path = "/groups/itay_mayrose/halabikeren/PloiDB/name_resolution/new_processed_resolved_names.csv"
    work_dir = "/groups/itay_mayrose/halabikeren/PloiDB/name_resolution/resolved_names_different_methods/"

    curator = ResolvedNamesCurator(unresolved_names_path=data_unresolved_names_path,
                                   resolved_names_paths=data_resolved_names_paths,
                                   selection_method=data_selection_method,
                                   reference_resolved_names=pd.read_csv(reference_resolved_names_path),
                                   reference_unresolved_names=pd.read_csv(reference_unresolved_names_path),
                                   work_dir=work_dir,
                                   is_processed=True)
    curator.process_resolved_names(output_path=data_output_path)