import matplotlib.pyplot as plt
import numpy as np
from fire import Fire
from matplotlib.colors import CenteredNorm, LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py


def heatmap(ax, data, **kwargs):
    imgplot = ax.imshow(data, origin="lower", **kwargs)

    ax.set_xlabel("$x$")
    ax.set_ylabel("$t$")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(imgplot, ax=ax, cax=cax)


@np.vectorize
def kink(t, x, V=0.0):
    gamma = 1 / np.sqrt(1 - V**2)
    x_prime = gamma * (x - V * t)

    if x_prime <= 0:
        return 0.0
    elif x_prime <= np.pi:
        return 1 - np.cos(x_prime)
    else:
        return 2.0


for l in np.linspace(0.5, 3, 6):
    for alpha in np.linspace(0, 0.75, 4):
        with h5py.File(f"data/{l=:0.2f},{alpha=:0.2f}.h5", "r") as f:
            field = np.array(f["field"])
            hamiltonian = np.array(f["hamiltonian"])
            t = np.array(f["t"])
            x = np.array(f["x"])

        extent = (x[0], x[-1], t[0], t[-1])
        figsize = (3.2, 3.1)

        fig, ax = plt.subplots(figsize=figsize, tight_layout=True, dpi=300)
        heatmap(ax, field, extent=extent)
        ax.set_xlim(-t[-1] / 2, t[-1] / 2)
        fig.savefig(f"fig/{l=:0.2f},{alpha=:0.2f}.pdf")
        fig.clear()

        fig, ax = plt.subplots(figsize=figsize, tight_layout=True, dpi=300)
        heatmap(
            ax,
            hamiltonian,
            extent=extent,
            norm=LogNorm(clip=True, vmin=1e-6),
            cmap="magma",
        )
        ax.set_xlim(-t[-1] / 2, t[-1] / 2)
        fig.savefig(f"fig/{l=:0.2f},{alpha=:0.2f}-energy.pdf")
        fig.clear()

        xx, tt = np.meshgrid(x, t)
        perturbation = field - kink(tt, xx + np.pi / 2)

        fig, ax = plt.subplots(figsize=figsize, tight_layout=True, dpi=300)
        heatmap(ax, perturbation, extent=extent, norm=CenteredNorm(), cmap="RdBu")
        ax.set_xlim(-t[-1] / 2, t[-1] / 2)
        fig.savefig(f"fig/{l=:0.2f},{alpha=:0.2f}-perturbation.pdf")
        fig.clear()

        plt.close("all")
