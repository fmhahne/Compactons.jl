using PythonCall, PythonPlot
const plt = PythonPlot.pyplot
const mpl = PythonPlot.matplotlib

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
