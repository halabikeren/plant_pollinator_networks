import logging
import os.path
import sys

import click
import pandas as pd

sys.path.append("..")
from services.pbs_service import PBSService

logger = logging.getLogger(__name__)

SCRIPT_PATH = "compute_shared_modules.R"


@click.command()
@click.option(
    "--pairs_data_path",
    help="directory to csv paths of the respective networks with their pairs data",
    type=str,
    required=True
)
@click.option(
    "--workdir",
    help="directory to create jobs in",
    type=str,
    required=False,
    default=f"{os.getcwd()}/compute_networks_shared_modules/"
)
@click.option(
    "--max_jobs_limit",
    help="maximal number of jobs to submit at once",
    type=int,
    required=False,
    default=1800
)
@click.option(
    "--queue",
    help="queue name to submit jobs to",
    type=str,
    required=False,
    default="itaym"
)
def compute_network_shared_modules(pairs_data_path: str,
                                   workdir: str,
                                   max_jobs_limit: int,
                                   queue: str):
    os.makedirs(workdir, exist_ok=True)
    log_path = f"{workdir}compute_network_shared_modules.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s module: %(module)s function: %(funcName)s line %(lineno)d: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(str(log_path)), ],
        force=True,  # run over root logger settings to enable simultaneous writing to both stdout and file handler
    )

    logger.info(f"writing jobs to {workdir}")
    pairs_data_dir = f"{workdir}/pairs_data/"
    os.makedirs(pairs_data_dir, exist_ok=True)
    jobs_dir = f"{workdir}/jobs/"
    os.makedirs(jobs_dir, exist_ok=True)
    jobs_output_dir = f"{workdir}/jobs_output/"
    os.makedirs(jobs_output_dir, exist_ok=True)

    pairs_commands = []
    job_names = []
    pairs_data = pd.read_csv(pairs_data_path)
    pairs_data_by_network = pairs_data.groupby("network_id")
    logger.info(
        f"creating jobs infrastructure feature computation across {len(pairs_data_by_network.groups.keys())} networks")
    for network_id in pairs_data_by_network.groups.keys():
        net_pairs_data_path = f"{pairs_data_dir}{network_id}.csv"
        pairs_data_by_network.get_group(network_id).to_csv(net_pairs_data_path)
        job_name = f"pairs_comp_{network_id}"
        cmds = [os.environ.get("CONDA_ACT_CMD", ""), f"cd {os.getcwd()}", f"Rscript --vanilla {SCRIPT_PATH} {net_pairs_data_path}"]
        pairs_commands.append(cmds)
        job_names.append(job_name)

    logger.info(f"executing modules sharing computation across {len(pairs_commands)} networks")
    res = PBSService.execute_job_array(work_dir=jobs_dir,
                                       output_dir=jobs_output_dir,
                                       job_ids=job_names,
                                       jobs_commands=pairs_commands,
                                       ram_gb_size=4,
                                       max_parallel_jobs=max_jobs_limit,
                                       queue=queue)

    all_pairs_data = []
    for network_id in pairs_data_by_network.groups.keys():
        net_pairs_data_path = f"{pairs_data_dir}{network_id}.csv"
        pairs_data = pd.read_csv(net_pairs_data_path)
        all_pairs_data.append(pairs_data)
    all_pairs_data = pd.concat(all_pairs_data)
    all_pairs_data.to_csv(pairs_data_path)



if __name__ == '__main__':
    compute_network_shared_modules()
