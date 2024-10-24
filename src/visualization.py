import wandb

from matplotlib import pyplot as plt


def visualize_rdf(data: dict):
    fig, ax = plt.subplots(figsize=(12.5, 5))
    ax.plot(
        data["rdf"][:, 0],
        data["rdf"][:, 1],
        label="rdf",
        color="red",
    )
    ax.plot(
        data["exp_rdf"][:, 0],
        data["exp_rdf"][:, 1],
        label="expt_Skinner_2013",
        linestyle="--",
        color="black",
    )
    ax.legend()
    ax.set_xlabel("#Distance_(Ang)")
    ax.set_ylabel("goo")
    wandb_log_dict = {"rdf_plot": wandb.Image(fig)}
    plt.close(fig)
    return wandb_log_dict
