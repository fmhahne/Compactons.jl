struct DeformedKink{T<:Real}
    ϵ::T
end

function simulation(params::DeformedKink; dx=1e-3, sampling=10)
    ϵ = params.ϵ
    tsave = 0.0:(dx * sampling):10.0

    x = (-tsave[end]):dx:tsave[end]
    xsave = x[begin:sampling:end]

    η₀ = kink.(x / (1.0 + ϵ) .+ π / 2)
    ∂ₜη₀ = zero(x)

    η, H = produce_data(quadratic, ∂ₜη₀, η₀, tsave; dx, sampling)
    return Dict("x" => xsave, "t" => tsave, "η" => η, "H" => H)
end
