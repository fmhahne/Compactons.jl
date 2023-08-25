function producedata(model, ∂ₜφ₀, φ₀, tsave; dx, dt=dx / 10, sampling=10, callbacks=[])
    savedhamiltonian = SavedValues(Float64, Vector{Float64})
    cbhamiltonian = SavingCallback(gethamiltonian, savedhamiltonian; saveat=tsave)

    callback = CallbackSet(cbhamiltonian, callbacks...)

    N = length(φ₀)
    tspan = (tsave[begin], tsave[end])

    prob = SecondOrderODEProblem(fieldeq!, ∂ₜφ₀, φ₀, tspan, (model, dx))
    sol = solve(
        prob,
        RK4();
        adaptive=false,
        dt=dt,
        saveat=tsave,
        save_idxs=(N + 1):sampling:(2N),
        callback=callback,
    )

    φ = reduce(hcat, sol.u)
    H = reduce(hcat, savedhamiltonian.saveval)

    return φ, H
end

simulation(parameters) = error("Simulation $(typeof(parameters)) not implemented")

struct NonBPSKink{T<:Real}
    ϵ::T
end

function simulation(parameters::NonBPSKink; dx=1e-3, sampling=10)
    ϵ = parameters.ϵ
    tsave = 0.0:(dx * sampling):10.0

    x = (-tsave[end]):dx:tsave[end]
    xsave = x[begin:sampling:end]

    η₀ = kink.(x / (1.0 + ϵ) .+ π / 2)
    ∂ₜη₀ = zero(x)

    η, H = producedata(quadratic, ∂ₜη₀, η₀, tsave; dx, sampling)
    return Dict("x" => xsave, "t" => tsave, "η" => η, "H" => H)
end

@with_kw struct KinkOscillon{T<:Real}
    l::T
    V::T
    α::T
    v₀::T
    x₀::T = x_R(α, V; l, v₀)
    model::Model = quadratic
end

function getenergies(u, t, integrator)
    φ = @views u[(end ÷ 2 + 1):end]
    ∂ₜφ = @views u[begin:(end ÷ 2)]

    N = length(φ)
    model, dx = integrator.p

    N₁ = N ÷ 2
    N₂ = N ÷ 2 + round(Int, π / dx)

    E₁ = dx * sum(𝒯(∂ₜφ[i], (φ[i + 1] - φ[i - 1]) / (2dx)) + model.V(φ[i]) for i in 2:N₁)
    E₂ =
        dx *
        sum(𝒯(∂ₜφ[i], (φ[i + 1] - φ[i - 1]) / (2dx)) + model.V(φ[i]) for i in (N₁ + 1):N₂) -
        π / 2
    E₃ =
        dx * sum(
            𝒯(∂ₜφ[i], (φ[i + 1] - φ[i - 1]) / (2dx)) + model.V(φ[i]) for
            i in (N₂ + 1):(N - 1)
        )

    return [E₁; E₂; E₃]
end

function simulation(parameters::KinkOscillon; dx=1e-3, sampling=10)
    @unpack l, V, α, v₀, x₀, model = parameters
    tsave = 0.0:(dx * sampling):10.0

    x = (-tsave[end]):dx:tsave[end]
    xsave = x[begin:sampling:end]

    η₀ = oscillon.(l * α * γ(V), x .+ x₀, V; l, v₀)
    ∂ₜη₀ = ∂ₜoscillon.(l * α * γ(V), x .+ x₀, V; l, v₀)

    η₀ += if model == quadratic
        kink.(0, x)
    elseif model == toy
        toykink.(0, x)
    else
        error("Kink not implemented")
    end

    energies = SavedValues(Float64, Vector{Float64})
    cbenergies = SavingCallback(getenergies, energies; saveat=tsave)

    η, H = producedata(model, ∂ₜη₀, η₀, tsave; dx, sampling, callbacks=[cbenergies])
    E = reduce(hcat, energies.saveval)

    return Dict(
        "x" => xsave,
        "t" => tsave,
        "η" => η,
        "H" => H,
        "E₁" => E[1, :],
        "E₂" => E[2, :],
        "E₃" => E[3, :],
    )
end

@with_kw struct KinkAntikink{T<:Real}
    V::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 10.0
    sampling::Int = 10
end

function simulation(parameters::KinkAntikink)
    @unpack V, dx, dt, tmax, sampling = parameters
    tsave = 0.0:(dx * sampling):tmax

    x = (-tmax):dx:tmax
    xsave = x[begin:sampling:end]

    η₀ = kink.(0.0, -abs.(x) .+ π / γ(V), V)
    ∂ₜη₀ = ∂ₜkink.(0, -abs.(x) .+ π / γ(V), V)

    η, H = producedata(quadratic, ∂ₜη₀, η₀, tsave; dx, dt, sampling)
    return Dict("x" => xsave, "t" => tsave, "η" => η, "H" => H)
end

@with_kw struct KinkKink{T<:Real}
    V::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 10.0
    sampling::Int = 10
end

function simulation(parameters::KinkKink; dx=1e-3, sampling=10)
    @unpack V, dx, dt, tmax, sampling = parameters
    tsave = 0.0:(dx * sampling):tmax

    x = (-tmax):dx:tmax
    xsave = x[begin:sampling:end]

    η₀ = @. kink(0.0, x + π / γ(V), V) + kink(0.0, x, -V) - 2
    ∂ₜη₀ = @. ∂ₜkink(0, x + π / γ(V), V) + ∂ₜkink(0.0, x, -V)

    η, H = producedata(quadratic, ∂ₜη₀, η₀, tsave; dx, sampling=sampling, dt=0.1dx)
    return Dict("x" => xsave, "t" => tsave, "η" => η, "H" => H)
end

function get_right_idx(u, t, integrator)
    η = @views u[(end ÷ 2 + 1):end]
    return findfirst(abs.(η .- 2) .< 1e-7)
end

@with_kw struct KinkKinkBorder{T<:Real}
    V::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 20.0
    xmax::T = 20.0
end

function simulation(parameters::KinkKinkBorder)
    @unpack V, dx, dt, tmax, xmax = parameters
    tsave = 0.0:1e-3:tmax
    x = -5.0:dx:5.0

    xR = Float64[]
    η₀ = @. kink(0.0, x + π / γ(V), V) + kink(0.0, x, -V) - 2
    ∂ₜη₀ = @. ∂ₜkink(0, x + π / γ(V), V) + ∂ₜkink(0.0, x, -V)

    xR_idxs = SavedValues(Float64, Int64)
    cbidxs = SavingCallback(get_right_idx, xR_idxs; saveat=tsave)

    prob = SecondOrderODEProblem(fieldeq!, ∂ₜη₀, η₀, (0.0, tmax), (quadratic, dx))
    sol = solve(
        prob, RK4(); adaptive=false, dt=0.1dx, save_everystep=false, callback=cbidxs
    )

    return Dict("x" => x, "xR_idxs" => xR_idxs)
end

@with_kw struct TriangularDelta{T<:Real}
    k::T
    ϵ::T
    A::T = 1.0
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 10.0
    sampling::Int = 10
end

function δ(x; ϵ)
    if -ϵ < x < ϵ
        return 1 / ϵ - abs(x) / ϵ^2
    else
        return zero(x)
    end
end

function simulation(parameters::TriangularDelta)
    @unpack k, ϵ, A, dx, dt, tmax, sampling = parameters

    tsave = 0.0:(dx * sampling):tmax

    x = (-tmax):dx:tmax
    xsave = x[begin:sampling:end]

    ϕ₀ = zero(x)
    ∂ₜϕ₀ = A * δ.(x; ϵ)

    ϕ, H = producedata(tanh_gordon(k), ∂ₜϕ₀, ϕ₀, tsave; dx, dt, sampling)
    return Dict("x" => xsave, "t" => tsave, "ϕ" => ϕ, "H" => H)
end
