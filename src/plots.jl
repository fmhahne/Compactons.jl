using PyPlot
using PyCall

const mpl = PyPlot.matplotlib

rcParams = PyDict(PyPlot.matplotlib."rcParams")

rcParams["figure.autolayout"] = true
rcParams["figure.figsize"] = 4.5, 2.78
rcParams["figure.dpi"] = 300
rcParams["font.size"] = 8
rcParams["legend.frameon"] = false
rcParams["lines.linewidth"] = 1
rcParams["axes.linewidth"] = 0.5
rcParams["grid.linewidth"] = 0.5
rcParams["xtick.major.width"] = 0.5
rcParams["ytick.major.width"] = 0.5
rcParams["xtick.minor.width"] = 0.5
rcParams["ytick.minor.width"] = 0.5

mpl.use("agg")

axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

function heatmap!(ax, x, t, data; kwargs...)
    extent = (x[1], x[end], t[1], t[end])
    implot = ax.imshow(data'; origin="lower", extent=extent, kwargs...)

    ax.set_xlabel("\$x\$")
    ax.set_ylabel("\$t\$")

    divider = axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size=0.1, pad=0.05)
    plt.colorbar(implot, ax=ax, cax=cax)

    nothing
end
