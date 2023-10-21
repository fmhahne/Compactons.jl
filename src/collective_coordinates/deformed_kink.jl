@with_kw struct CCDeformedKink{T<:Real}
    ϵ::T
    tmax::T = 10.0
    saveat::T = 0.01
    quadgk_kwargs = (atol=1e-6, order=20)
end

function collective_coordinates(params::CCDeformedKink)
    @unpack ϵ, tmax, quadgk_kwargs, saveat = params

    function χ(x)
        if -π / 2 < x < π / 2
            return sin(2 * x) + 0.5 * sin(4 * x)
        else
            return 0.0
        end
    end

    function η(x, q)
        b, c = q
        return kink(b * x + π / 2) + c * χ(b * x)
    end

    function xspan(q)
        b, c = q
        return (-π / (2 * b), π / (2 * b))
    end

    q̇₀ = [0.0, 0.0]
    q₀ = [1 / (1 + ϵ), 0.0]

    prob = SecondOrderODEProblem(
        cceq!, q̇₀, q₀, (0.0, tmax), (quadratic, η, xspan, quadgk_kwargs)
    )
    sol = solve(prob; saveat, save_idxs=[3, 4])

    return Dict("t" => sol.t, "b" => sol[1, :], "c" => sol[2, :])
end
