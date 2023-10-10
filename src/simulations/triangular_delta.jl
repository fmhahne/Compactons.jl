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

    ϕ, H = produce_data(tanh_gordon(k), ∂ₜϕ₀, ϕ₀, tsave; dx, dt, sampling)
    return Dict("x" => xsave, "t" => tsave, "ϕ" => ϕ, "H" => H)
end
