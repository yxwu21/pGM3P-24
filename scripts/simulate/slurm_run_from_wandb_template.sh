python simulate.py \
    --md-param-key {WANDB_LOG_NAME} \
    --md-param-step {LOG_STEP} \
    --directory-path "/home/yxwu/pGM_water_model/wandb_data/{WANDB_LOG_NAME}" \
    --md-config.output-dir /home/yxwu/pGM_water_model/MD_data/{THERMOSTAT}_{BAROSTAT}/{PIPELINE_NAME} \
    --md-config.case case_1 \
    --md-config.md-lengths {MD_LENGTHS} \
    --md-config.use-custom-case-folder-name \
    --md-config.custom-case-folder-name {WANDB_LOG_NAME}_{LOG_STEP} \
    --md-config.thermostat {THERMOSTAT} \
    --md-config.barostat {BAROSTAT} \
    --md-config.pipeline-name {PIPELINE_NAME} \
    --md-config.gamma-lns {GAMMA_LNS} \
    --slurm.mode slurm \
    --slurm.slurm-job-name {WANDB_LOG_NAME}_{LOG_STEP} \
    --slurm.slurm-partition "cpu" \
    --slurm.slurm-output-folder "/home/yxwu/pGM_water_model/MD_data/{THERMOSTAT}_{BAROSTAT}/{PIPELINE_NAME}/{WANDB_LOG_NAME}_{LOG_STEP}/slurm" \
    --slurm.tasks-per-node 1 \
    --slurm.cpus-per-task 2 \
    --slurm.mem 1G