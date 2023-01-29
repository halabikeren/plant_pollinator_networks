import logging
import sys
from typing import Optional
import click
import pandas as pd
import os
import sys
import glob
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from services.pbs_service import PBSService
logger = logging.getLogger(__name__)

from dotenv import load_dotenv, find_dotenv
load_dotenv(find_dotenv())


@click.command()
@click.option(
    "--networks_dir",
    help="directory to the networks files",
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
    "--features_path",
    help="path to hold the features across all networks",
    type=click.Path(exists=False),
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
    default=1900, # for API request, this limitation will prevent too many simultaneous requests issue
)
@click.option(
    "--queue",
    help="queue to submit jobs to",
    type=str,
    required=False,
    default="itaym",
)
def exec_by_network(networks_dir: str,
                    script_path: str,
                    work_dir: str,
                    features_path: str,
                    log_path: str,
                    max_jobs_in_parallel: int,
                    queue: str):

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s module: %(module)s function: %(funcName)s line %(lineno)d: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(str(log_path)), ],
        force=True,  # run over root logger settings to enable simultaneous writing to both stdout and file handler
    )

    output_dir = f"{work_dir}/features_by_network/"
    os.makedirs(output_dir, exist_ok=True)
    networks_commands = []
    logger.info(f"generating jobs for networks in {work_dir}jobs/")
    for network_path in os.listdir(networks_dir):
        if "reference" not in network_path:
            input_path = f"{networks_dir}{network_path}"
            output_path = f"{output_dir}{network_path.replace('.csv', '_features.csv')}"
            script_exec_cmd = f"Rscript --vanilla {script_path} {input_path} {output_path}"
            if script_path.endswith(".py"):
                output_path = f"{output_dir}{network_path.replace('.csv', '_features/')}"
                script_exec_cmd = f"python {script_path} --input_path={input_path} --output_dir={output_path}"
            if not os.path.exists(output_path):
                commands = [os.environ.get("CONDA_ACT_CMD", ""), script_exec_cmd]
                networks_commands.append(commands)

    logger.info(f"executing feature computation across {len(networks_commands)} networks")
    res = PBSService.execute_job_array(work_dir=f"{work_dir}jobs/",
                                 output_dir=f"{work_dir}jobs_output/",
                                 jobs_commands=networks_commands,
                                 ram_gb_size=1,
                                 max_parallel_jobs=max_jobs_in_parallel,
                                 queue=queue)
    
    if script_path.endswith(".py"):
        network_features = pd.concat([pd.read_csv(f"{output_dir}{path}/network_features.csv") for path in os.listdir(output_dir) if os.path.exists(f"{output_dir}{path}/network_features.csv")])
        network_features.to_csv(f"{features_path}/network_features.csv", index=False)
        plant_features = pd.concat([pd.read_csv(f"{output_dir}{path}/plant_features.csv") for path in os.listdir(output_dir) if os.path.exists(f"{output_dir}{path}/plant_features.csv")])
        plant_features.to_csv(f"{features_path}/plant_features.csv", index=False)
    else:
        features = pd.concat([pd.read_csv(f"{output_dir}{path}") for path in os.listdir(output_dir)])
        features.to_csv(features_path, index=False)
        
if __name__ == '__main__':
    exec_by_network()

