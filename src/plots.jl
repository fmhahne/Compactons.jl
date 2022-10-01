ENV["QT_QPA_PLATFORM"] = "wayland"

using PyPlot
using PyCall

const mpl = PyPlot.matplotlib

rcParams = PyDict(PyPlot.matplotlib."rcParams")

rcParams["figure.autolayout"] = true
rcParams["figure.figsize"] = 4.5, 2.78
rcParams["figure.dpi"] = 250
rcParams["font.size"] = 8
rcParams["legend.frameon"] = false
rcParams["lines.linewidth"] = 1
rcParams["axes.linewidth"] = 0.5
rcParams["grid.linewidth"] = 0.5
rcParams["xtick.major.width"] = 0.5
rcParams["ytick.major.width"] = 0.5
rcParams["xtick.minor.width"] = 0.5
rcParams["ytick.minor.width"] = 0.5

axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

function heatmap!(ax, x, t, data; colorbar=false, kwargs...)
    extent = (x[begin], x[end], t[begin], t[end])
    implot = ax.imshow(data'; origin="lower", extent=extent, kwargs...)
    ax.set_xlabel(raw"$x$")
    ax.set_ylabel(raw"$t$")

    if colorbar
        divider = axes_grid1.make_axes_locatable(ax)
        cax = divider.append_axes("right"; size=0.1, pad=0.05)
        plt.colorbar(implot; ax=ax, cax=cax)
    end

    return nothing
end

function shared_colorbar!(fig, cb)
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    cax = fig.add_axes([0.87, 0.16, 0.015, 0.74])
    fig.colorbar(cb, cax=cax)

    return nothing
end

function show_kink_borders!(ax; color="black")
    ax.axvline(0; linewidth=0.5, color=color, linestyle="dashed")
    ax.axvline(float(π); linewidth=0.5, color=color, linestyle="dashed")

    return nothing
end
