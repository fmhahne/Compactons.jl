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


for V in np.linspace(0, 0.99, 100):
    with h5py.File(f"data/{V=:0.2f}.h5", "r") as f:
        field = np.array(f["field"])
        hamiltonian = np.array(f["hamiltonian"])
        t = np.array(f["t"])
        x = np.array(f["x"])

    extent = (x[0], x[-1], t[0], t[-1])
    figsize = (3.2, 3.1 * t[-1] / 20.0)

    fig, ax = plt.subplots(figsize=figsize, tight_layout=True, dpi=300)
    heatmap(ax, field, extent=extent)
    ax.set_xlim(-10, 10)
    fig.savefig(f"fig/{V=:0.2f}.pdf")
    fig.clear()

    fig, ax = plt.subplots(figsize=figsize, tight_layout=True, dpi=300)
    heatmap(
        ax, hamiltonian, extent=extent, norm=LogNorm(clip=True, vmin=1e-6), cmap="magma"
    )
    ax.set_xlim(-10, 10)
    fig.savefig(f"fig/{V=:0.2f}-energy.pdf")
    fig.clear()

    plt.close("all")
