@with_kw struct SymOscillonScattering{T<:Real}
    v::T
    α::T
    v₀::T
    dx = 1e-3
    dt = 0.1 * dx
    tmax = 10.0
    sampling = 10
end

function simulation(params::SymOscillonScattering)
    @unpack v, α, v₀, dt, dx, sampling = params
    tsave = 0.0:(dx * sampling):10.0

    x = (-tsave[end]):dx:tsave[end]
    xsave = x[begin:sampling:end]

    Δx = Compactons.x_R(α, v; v₀)

    ϕ₀ = oscillon.(α * γ(v), x .+ Δx, v; v₀) + oscillon.(α * γ(v), -x .+ Δx, v; v₀)
    ∂ₜϕ₀ = ∂ₜoscillon.(α * γ(v), x .+ Δx, v; v₀) + ∂ₜoscillon.(α * γ(v), -x .+ Δx, v; v₀)

    ϕ, H = produce_data(signum_gordon, ∂ₜϕ₀, ϕ₀, tsave; dx, sampling)

    return Dict("x" => xsave, "t" => tsave, "ϕ" => ϕ, "H" => H)
end
