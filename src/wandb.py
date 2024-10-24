import wandb
import smac
import os

from smac import Callback
from smac.runhistory import TrialInfo, TrialValue
from dataclasses import asdict, dataclass

os.environ["WANDB__SERVICE_WAIT"] = "300"


@dataclass
class WandBConfig:
    project: str = "blackbox_water"
    entity: str = "none"


class WandBCallback(Callback):
    def __init__(
        self,
        wandb_config: WandBConfig,
        config_to_save: dict = None,
        wandb_data_key="wandb",
        wandb_function=None,
        wandb_function_step: int = 50,
    ) -> None:
        self.wandb = wandb.init(
            project=wandb_config.project,
            entity=wandb_config.entity,
            config=config_to_save,
        )
        self.wandb_buffer = {}
        self.wandb_logging_key = wandb_data_key
        self.wandb_function = wandb_function
        self.wandb_function_step = wandb_function_step

        # wandb states
        self.best_config_wandb = {}
        self.step = 1

        super().__init__()

    def log(self, data: dict) -> None:
        self.wandb.log(data, step=self.step)

    def is_best_from_buffer(self, best_config_id: int) -> bool:
        # update best if necessary
        is_best_from_buffer = False
        if best_config_id in self.wandb_buffer:
            self.best_config_wandb[best_config_id] = self.wandb_buffer[best_config_id]
            self.wandb_buffer.clear()
            is_best_from_buffer = True

        return is_best_from_buffer

    def on_tell_start(
        self, smbo: smac.main.smbo.SMBO, info: TrialInfo, value: TrialValue
    ) -> bool | None:
        # pop out wandb specific keys into the buffer, because this information might #
        # not be able to be serialized e.g. image, etc.
        wandb_data = value.additional_info.pop(self.wandb_logging_key, {})
        self.wandb_buffer[info.config.config_id] = wandb_data

        info_dict = asdict(info)
        info_dict["config"] = dict(info_dict["config"])
        value_dict = asdict(value)

        # log state
        log_dict = {f"trial_state/{k}": v for k, v in (info_dict | value_dict).items()}

        # log config id
        log_dict["trial_state/config_id"] = info.config.config_id

        # log visualization
        if (
            self.wandb_function
            and wandb_data
            and self.step % self.wandb_function_step == 0
        ):
            for k, v in self.wandb_function(wandb_data).items():
                log_dict[f"trial_state/{k}"] = v

        self.log(log_dict)
        return super().on_tell_start(smbo, info, value)

    def on_tell_end(
        self, smbo: smac.main.smbo.SMBO, info: TrialInfo, value: TrialValue
    ) -> bool | None:
        # log the best result until now
        log_dict = {}
        trajectory = smbo.intensifier.trajectory
        if trajectory:
            best_traj = trajectory[-1]

            best_config_ids = best_traj.config_ids
            best_costs = best_traj.costs
            log_dict["best/config_num"] = len(best_config_ids)

            for i, (best_config_id, best_cost) in enumerate(
                zip(best_config_ids, best_costs)
            ):
                suffix = f"_{i}" if i > 0 else ""

                log_dict[f"best/config_id{suffix}"] = best_config_id
                log_dict[f"best/cost{suffix}"] = best_cost

                # Get the best configuration
                best_config = smbo.runhistory.get_config(best_config_id)
                for k, v in dict(best_config).items():
                    log_dict[f"best/{k}{suffix}"] = v

                # Log wandb specific information
                if self.wandb_function and self.is_best_from_buffer(best_config_id):
                    # if contains wandb specific information
                    wandb_data = self.best_config_wandb[best_config_id]
                    if wandb_data:
                        for k, v in self.wandb_function(wandb_data).items():
                            log_dict[f"best/{k}{suffix}"] = v

        self.log(log_dict)

        # update step when logging is done
        self.step += 1
        return super().on_tell_end(smbo, info, value)
