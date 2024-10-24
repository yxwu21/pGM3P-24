import dask

from typing import List
from dask_jobqueue import PBSCluster
from dask.distributed import Client
from dataclasses import dataclass, field


@dataclass
class ClusterConfigurations:
    queue: str = "all.q"

    name: str = "blackbox_water"

    cores: int = 4

    memory: str = "512GB"

    walltime: str = "00:10:00"

    job_script_prologue: List[str] = field(
        default_factory=lambda: ["cd /home/yxwu/blackbox_water_mirror"]
    )

    job_extra_directives: List[str] = field(default_factory=lambda: [])

    log_directory: str = "outputs/log"

    # total job number
    n_jobs: int = 4


class Cluster:
    def __init__(
        self,
        params: ClusterConfigurations,
    ) -> None:
        self.params = params

    def get_client(self):
        dask.config.set({"distributed.worker.daemon": False})
        cluster = PBSCluster(
            queue=self.params.queue,
            name=self.params.name,
            cores=self.params.cores,
            memory=self.params.memory,
            walltime=self.params.walltime,
            job_script_prologue=self.params.job_script_prologue,
            job_extra_directives=self.params.job_extra_directives,
            log_directory=self.params.log_directory,
            # local_directory="/home/yxwu/blackbox_water",
            death_timeout=3600 * 24 * 15,  # timeout after 15 days
        )
        print(cluster.job_script())

        cluster.scale(jobs=self.params.n_jobs)  # ask for n jobs
        client = Client(cluster)
        client.wait_for_workers(max(self.params.n_jobs // 8, 1))
        return client
