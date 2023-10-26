import logging

import click
import pandas as pd
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from services.pbs_service import PBSService

logger = logging.getLogger(__name__)

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())


@click.command()
@click.option(
    "--network_path",
    help="path to a specific networks file",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--null_networks_dir",
    help="directory to the null networks simulated based on the given network",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--script_path",
    help="path of script to execute per network",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
)
@click.option(
    "--work_dir",
    help="directory of intermediate io",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--log_path",
    help="path to log file",
    type=click.Path(exists=False),
    required=True,
)
@click.option(
    "--max_jobs_in_parallel",
    help="maximal number of jobs that should be executed in parallel",
    type=int,
    required=False,
    default=1900,  # for API request, this limitation will prevent too many simultaneous requests issue
)
@click.option(
    "--queue",
    help="queue to submit jobs to",
    type=str,
    required=False,
    default="itaym",
)
@click.option(
    "--mem_per_job",
    help="amount of memory in gb to assign to each single-network job",
    type=int,
    required=False,
    default=10,
)
def computer_features_across_queue(network_path: str,
                                    null_networks_dir: str,
                                    script_path: str,
                                    work_dir: str,
                                    log_path: str,
                                    max_jobs_in_parallel: int,
                                    queue: str,
                                    mem_per_job: int):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s module: %(module)s function: %(funcName)s line %(lineno)d: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(str(log_path)), ],
        force=True,  # run over root logger settings to enable simultaneous writing to both stdout and file handler
    )

    logger.info(f"working in {os.getcwd()}")

    alt_network_id = os.path.basename(network_path).split(".")[0]

    features_dir = f"{work_dir}/features_by_network/"
    output_dir = f"{features_dir}{alt_network_id}/"
    os.makedirs(output_dir, exist_ok=True)
    network_paths = [network_path] + [f"{null_networks_dir}/{p}" for p in os.listdir(null_networks_dir)]
    job_ids, networks_commands = [], []
    for path in network_paths:
        network_index = os.path.basename(path).replace(".csv", "")
        job_ids.append(network_index)
        if path == network_path:
            output_path = f"{output_dir}{network_index}_features.csv"
        else:
            output_path = f"{output_dir}{network_index}_null_features.csv"
        script_exec_cmd = f"Rscript --vanilla {script_path} {path} {null_networks_dir} {output_path} TRUE"
        if script_path.endswith(".py"):
            script_exec_cmd = f"python {script_path} --input_path={path} --output_dir={output_path}"
        if not os.path.exists(output_path):
            commands = [os.environ.get("CONDA_ACT_CMD", ""), f"cd {os.getcwd()}", script_exec_cmd]
            networks_commands.append(commands)

    logger.info(f"executing feature computation across {len(networks_commands)} networks corresponding to the "
                f"alternative and null networks of {alt_network_id}")
    res = PBSService.execute_job_array(work_dir=f"{work_dir}jobs/{alt_network_id}",
                                       output_dir=f"{work_dir}jobs_output/{alt_network_id}",
                                       job_ids=job_ids,
                                       jobs_commands=networks_commands,
                                       ram_gb_size=mem_per_job,
                                       max_parallel_jobs=max_jobs_in_parallel,
                                       queue=queue)

    res = os.system(f"mv {output_dir}{alt_network_id}_features.csv {features_dir}")
    alternative_network_features = pd.read_csv(f"{features_dir}{alt_network_id}_features.csv")
    null_network_features_paths = [f"{output_dir}/{p}" for p in os.listdir(output_dir) if "null" in p]
    null_network_features = pd.concat([pd.read_csv(p) for p in null_network_features_paths])
    null_network_features.to_csv(f"{features_dir}{alt_network_id}_features_across_null_networks.csv")
    # delta transform using null features
    features_to_transform = set(alternative_network_features.columns) - {"network", "degree", "d", "normalized.degree",
                                                                         "nestedness_contribution"}
    for feature in features_to_transform:
        alternative_network_features[f"standardized_{feature}"] = alternative_network_features[feature]-null_network_features[feature].mean()
    alternative_network_features.to_csv(f"{features_dir}{alt_network_id}_features.csv")


if __name__ == '__main__':
    computer_features_across_queue()



for p in to_make_weighted:
    new_index = nt_to_free_indices["weighted"].pop()
    new_p = f"/groups/itay_mayrose/halabikeren/plant_pollinator_networks/data/networks/fixed_all/weighted/{new_index}.csv"
    assert(not os.path.exists(new_p))
    assert(new_index not in nt_to_free_indices["weighted"])
    old_p = p.replace("../../","/groups/itay_mayrose/halabikeren/plant_pollinator_networks/")
    old_to_new_path[old_p] = new_p