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

    prob = SecondOrderODEProblem(field_equation!, ∂ₜη₀, η₀, (0.0, tmax), (quadratic, dx))
    sol = solve(
        prob, RK4(); adaptive=false, dt=0.1dx, save_everystep=false, callback=cbidxs
    )

    return Dict("x" => x, "xR_idxs" => xR_idxs)
end
