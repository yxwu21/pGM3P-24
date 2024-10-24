from smac import Callback
from smac.main.smbo import SMBO
from smac.runner.dask_runner import DaskParallelRunner
from dask.distributed import wait
from smac.utils.logging import get_logger

logger = get_logger(__name__)


class WaitRunningTrialsAfterSubmmision(Callback):
    def __init__(
        self, wait_submitted_num: int = 10, wait_pending_trails_num: int = None
    ) -> None:
        super().__init__()
        self.wait_submitted_num = wait_submitted_num
        self.wait_pending_trails_num = wait_pending_trails_num

    def on_ask_start(self, smbo: SMBO):
        submmited_trial_num = smbo.runhistory.submitted
        running_trial_num = smbo.runhistory.running
        if (
            isinstance(smbo._runner, DaskParallelRunner)
            and (submmited_trial_num + 1) % self.wait_submitted_num == 0
        ):
            logger.info(
                f"Waiting all {running_trial_num} running trials after submmited {submmited_trial_num} trials."
            )
            wait(
                smbo._runner._pending_trials[: self.wait_pending_trails_num],
                return_when="ALL_COMPLETED",
            )
            logger.info(
                f"Continuing trials after finished {running_trial_num if self.wait_pending_trails_num is None else self.wait_pending_trails_num} running trials."
            )
