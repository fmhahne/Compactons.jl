@with_kw struct KinkAntikinkKink{T<:Real}
    vL::T
    vR::T
    dx::T = 1e-3
    dt::T = 0.1dx
    tmax::T = 10.0
    sampling::Int = 10
end

function simulation(params::KinkAntikinkKink)
    @unpack vL, vR, dx, dt, tmax, sampling = params
    tsave = 0.0:(dx * sampling):tmax

    x = (-tmax - π / 2):dx:(tmax + π / 2)
    xsave = x[begin:sampling:end]

    η₀ =
        kink.(0.0, @.(x + π / γ(vL) + π / 2), vL) - kink.(0.0, @.(x + π / 2)) +
        kink.(0.0, @.(x - π / 2), vR)
    ∂ₜη₀ =
        ∂ₜkink.(0.0, @.(x + π / γ(vL) + π / 2), vL) - ∂ₜkink.(0.0, @.(x + π / 2)) +
        ∂ₜkink.(0.0, @.(x - π / 2), vR)

    η, H = produce_data(quadratic, ∂ₜη₀, η₀, tsave; dx, dt, sampling)
    return Dict("x" => xsave, "t" => tsave, "η" => η, "H" => H)
end
