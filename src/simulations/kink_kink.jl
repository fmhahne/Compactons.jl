@with_kw struct KinkKink{T<:Real}
    v::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 10.0
    sampling::Int = 10
end

function simulation(params::KinkKink)
    @unpack v, dx, dt, tmax, sampling = params
    tsave = 0.0:(dx * sampling):tmax

    x = (-tmax):dx:tmax
    xsave = x[begin:sampling:end]

    η₀ = @. kink(0.0, x + π / γ(v), v) + kink(0.0, x, -v) - 2
    ∂ₜη₀ = @. ∂ₜkink(0, x + π / γ(v), v) + ∂ₜkink(0.0, x, -v)

    η, H = produce_data(quadratic, ∂ₜη₀, η₀, tsave; dx, sampling, dt)
    return Dict("x" => xsave, "t" => tsave, "η" => η, "H" => H)
end

function get_right_idx(u, t, integrator)
    η = @views u[(end ÷ 2 + 1):end]
    return findfirst(abs.(η .- 2) .< 1e-7)
end

@with_kw struct KinkKinkBorder{T<:Real}
    v::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 20.0
    xmax::T = 20.0
end

function simulation(params::KinkKinkBorder)
    @unpack v, dx, dt, tmax, xmax = params
    tsave = 0.0:1e-3:tmax
    x = -5.0:dx:5.0

    η₀ = @. kink(0.0, x + π / γ(v), v) + kink(0.0, x, -v) - 2
    ∂ₜη₀ = @. ∂ₜkink(0, x + π / γ(v), v) + ∂ₜkink(0.0, x, -v)

    xR_idxs = SavedValues(Float64, Int64)
    cbidxs = SavingCallback(get_right_idx, xR_idxs; saveat=tsave)

    prob = SecondOrderODEProblem(field_equation!, ∂ₜη₀, η₀, (0.0, tmax), (quadratic, dx))
    solve(prob, RK4(); adaptive=false, dt, save_everystep=false, callback=cbidxs)

    return Dict("x" => x, "xR_idxs" => xR_idxs)
end
