using Compactons
using JLD2, CodecZlib
using DrWatson
include(srcdir("plots.jl"))

function fill_image!(image, energies, i, j; l, V, α, v₀)
    @unpack E₁, E₂, E₃ = energies[(l=l, V=V, α=α, v₀=v₀)]
    nc = findfirst(E₁ .< 0.1E₁[begin])
    Δt = 10.0 / length(E₁)

    if isnothing(nc)
        image[i, j] = NaN
    else
        n = nc + round(Int, (π + 1) / Δt)
        image[i, j] = n ≤ length(E₁) ? E₃[n] / E₁[begin] : NaN
    end
end

ls = 0.50:0.025:3.00
Vs = 0.0:1e-2:0.99
αs = 0.0:1e-2:0.99
v₀s = 0.0:1e-2:0.99

fig, axs = plt.subplots(2, 3; figsize=(6.2, 6.2 * 2 / (1 + √5)), tight_layout=false)

norm = mpl.colors.Normalize()
cmap = mpl.cm.get_cmap("plasma")
cb = mpl.cm.ScalarMappable(; norm=norm, cmap=cmap)

let l = 1.0, V = 0.75, ax = axs[1, 1]
    @unpack energies = load(datadir("kink_oscillon_scattering_energies", "l=$l,V=$V.jld2"))

    image = zeros(length(αs), length(v₀s))
    for (j, v₀) in enumerate(v₀s)
        for (i, α) in enumerate(αs)
            fill_image!(image, energies, i, j; l, V, α, v₀)
        end
    end

    ax.imshow(
        image;
        origin="lower",
        extent=[αs[begin], αs[end], v₀s[begin], v₀s[end]],
        cmap=cmap,
        aspect="auto",
        norm=norm,
    )
    ax.set_xlabel(raw"$v_0$")
    ax.set_ylabel(raw"$\alpha$")
    ax.set_title("\$ l=$l \$, \$ V=$V \$")

    cb.set_array(image)
    cb.autoscale()
end

let l = 1.0, α = 0.0, ax = axs[1, 2]
    @unpack energies = load(
        datadir("kink_oscillon_scattering_energies", "l=$l,alpha=$α.jld2")
    )

    image = zeros(length(v₀s), length(Vs))
    for (j, V) in enumerate(Vs)
        for (i, v₀) in enumerate(v₀s)
            fill_image!(image, energies, i, j; l, V, α, v₀)
        end
    end

    ax.imshow(
        image;
        origin="lower",
        extent=[v₀s[begin], v₀s[end], Vs[begin], Vs[end]],
        cmap=cmap,
        aspect="auto",
        norm=norm,
    )
    ax.set_xlabel(raw"$V$")
    ax.set_ylabel(raw"$v_0$")
    ax.set_title("\$ l=$l \$, \$ \\alpha=$α \$")

    cb.set_array(image)
    cb.autoscale()
end

let l = 1.0, v₀ = 0.0, ax = axs[1, 3]
    @unpack energies = load(
        datadir("kink_oscillon_scattering_energies", "l=$l,v0=$v₀.jld2")
    )

    image = zeros(length(αs), length(Vs))
    for (j, V) in enumerate(Vs)
        for (i, α) in enumerate(αs)
            fill_image!(image, energies, i, j; l, V, α, v₀)
        end
    end

    ax.imshow(
        image;
        origin="lower",
        extent=[αs[begin], αs[end], Vs[begin], Vs[end]],
        cmap=cmap,
        aspect="auto",
        norm=norm,
    )
    ax.set_xlabel(raw"$V$")
    ax.set_ylabel(raw"$\alpha$")
    ax.set_title("\$ l=$l \$, \$ v_0=$v₀ \$")

    cb.set_array(image)
    cb.autoscale()
end

let V = 0.75, α = 0.0, ax = axs[2, 1]
    @unpack energies = load(
        datadir("kink_oscillon_scattering_energies", "V=$V,alpha=$α.jld2")
    )

    image = zeros(length(v₀s), length(ls))
    for (j, l) in enumerate(ls)
        for (i, v₀) in enumerate(v₀s)
            fill_image!(image, energies, i, j; l, V, α, v₀)
        end
    end

    ax.imshow(
        image;
        origin="lower",
        extent=[ls[begin], ls[end], v₀s[begin], v₀s[end]],
        cmap=cmap,
        aspect="auto",
        norm=norm,
    )
    ax.set_xlabel(raw"$l$")
    ax.set_ylabel(raw"$v_0$")
    ax.set_title("\$ V=$V \$, \$ \\alpha=$α \$")

    cb.set_array(image)
    cb.autoscale()
end

let V = 0.75, v₀ = 0.0, ax = axs[2, 2]
    @unpack energies = load(
        datadir("kink_oscillon_scattering_energies", "V=$V,v0=$v₀.jld2")
    )

    image = zeros(length(αs), length(ls))
    for (j, l) in enumerate(ls)
        for (i, α) in enumerate(αs)
            fill_image!(image, energies, i, j; l, V, α, v₀)
        end
    end

    ax.imshow(
        image;
        origin="lower",
        extent=[ls[begin], ls[end], αs[begin], αs[end]],
        cmap=cmap,
        aspect="auto",
        norm=norm,
    )
    ax.set_xlabel(raw"$l$")
    ax.set_ylabel(raw"$\alpha$")
    ax.set_title("\$ V=$V \$, \$ v_0=$v₀ \$")

    cb.set_array(image)
    cb.autoscale()
end

let α = 0.00, v₀ = 0.0, ax = axs[2, 3]
    @unpack energies = load(
        datadir("kink_oscillon_scattering_energies", "alpha=$α,v0=$v₀.jld2")
    )

    image = zeros(length(Vs), length(ls))
    for (j, l) in enumerate(ls)
        for (i, V) in enumerate(Vs)
            fill_image!(image, energies, i, j; l, V, α, v₀)
        end
    end

    ax.imshow(
        image;
        origin="lower",
        extent=[ls[begin], ls[end], Vs[begin], Vs[end]],
        cmap=cmap,
        aspect="auto",
        norm=norm,
    )
    ax.set_xlabel(raw"$l$")
    ax.set_ylabel(raw"$V$")
    ax.set_title("\$ \\alpha=$α \$, \$ v_0=$v₀ \$")

    cb.set_array(image)
    cb.autoscale()
end

fig.tight_layout()
fig.subplots_adjust(; right=0.85)
cax = fig.add_axes([0.88, 0.15, 0.015, 0.7])
fig.colorbar(cb; cax=cax)
mkpath(plotsdir("kink_oscillon_scattering"))
fig.savefig(plotsdir("kink_oscillon_scattering", "energies.pdf"))
fig
