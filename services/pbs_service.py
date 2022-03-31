import logging
import os
import re
import typing as t
import getpass
import subprocess
import shutil
from time import sleep

import numpy as np

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

logger = logging.getLogger(__name__)


class PBSService:
    @staticmethod
    def create_job_file(
            job_path,
            job_name: str,
            job_output_dir: str,
            commands: t.List[str],
            queue: str = "itaym",
            priority: int = 0,
            cpus_num: int = 1,
            ram_gb_size: int = 4,
    ) -> int:
        os.makedirs(os.path.dirname(job_path), exist_ok=True)
        os.makedirs(os.path.dirname(job_output_dir), exist_ok=True)
        commands_str = "\n".join(commands)
        job_content = f"""# !/bin/bash
    #PBS -S /bin/bash
    #PBS -j oe
    #PBS -r y
    #PBS -q {queue}
    #PBS -p {priority}
    #PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
    #PBS -N {job_name}
    #PBS -e {job_output_dir}
    #PBS -o {job_output_dir}
    #PBS -r y
    #PBS -l select=ncpus={cpus_num}:mem={ram_gb_size}gb
    {commands_str}
    """
        with open(job_path, "w") as outfile:
            outfile.write(job_content)

        return 0

    @staticmethod
    def compute_curr_jobs_num() -> int:
        """
        :return: returns the current number of jobs under the shell username
        """
        username = getpass.getuser()
        proc = subprocess.run(f"qselect -u {username} | wc -l", shell=True, check=True, capture_output=True)
        curr_jobs_num = int(proc.stdout)
        return curr_jobs_num

    @staticmethod
    def execute_job_array(
            work_dir: str,
            output_dir: str,
            jobs_commands: t.List[t.List[str]],
            max_parallel_jobs: int = 1900,
    ):
        # set env: tree should be identical between input dir, work dir and output dir
        os.makedirs(work_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)

        # traverse input dir and for each create job file and job output dir with command that pushes output to
        # output file
        logger.info(f"# input paths to execute commands on = {len(jobs_commands)}")
        job_paths, job_output_paths = [], []
        for i in range(len(jobs_commands)):
            job_path = f"{work_dir}/{i}.sh"
            job_name = f"{i}.sh"
            job_output_path = f"{output_dir}/{i}.out"
            PBSService.create_job_file(
                job_path=job_path,
                job_name=job_name,
                job_output_dir=job_output_path,
                commands=[
                             os.environ.get(
                                 "conda_act_cmd",
                                 "source /groups/itay_mayrose/halabikeren/miniconda3/etc/profile.d/conda.sh; conda activate ppn")
                         ] + jobs_commands[i],
            )
            job_paths.append(job_path)
            job_output_paths.append(job_output_path)

        logger.info(f"# jobs to submit = {len(job_paths)}")

        # submit all the jobs, while maintaining the limit number on parallel jobs
        job_index = 0
        job_ids = []
        while job_index < len(job_paths):
            while PBSService.compute_curr_jobs_num() > max_parallel_jobs:
                sleep(2 * 60)
            try:
                res = subprocess.check_output(['qsub', f'{job_paths[job_index]}'])
                job_ids.append(re.search("(\d+)\.power\d", str(res)).group(1))
                job_index += 1
            except Exception as e:
                logger.error(f"failed to submit job at index {job_index} due to error {e}")
                exit(1)
            if job_index % 500 == 0:
                logger.info(f"submitted {job_index} jobs thus far")

        # wait for jobs to finish
        jobs_complete = np.all(os.system(f"stat -f {job_id}") != 0 for job_id in job_ids)
        while not jobs_complete:
            sleep(2 * 60)
            jobs_complete = np.all(os.system(f"stat -f {job_id}") != 0 for job_id in job_ids)

        # remove work dir
        shutil.rmtree(work_dir, ignore_errors=True)
