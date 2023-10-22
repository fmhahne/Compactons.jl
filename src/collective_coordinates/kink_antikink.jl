@with_kw struct CCKinkAntikink{T<:Real}
    v::T
    tmax::T = 10.0
    saveat::T = 0.01
    method = VelocityVerlet()
    solver_kwargs = (; dt=0.005)
    quadgk_kwargs = (; rtol=1e-5, atol=1e-5)
end

function internal_mode(x)
    if -π / 2 < x < π / 2
        return sin(2 * x) + 0.5 * sin(4 * x)
    else
        return 0.0
    end
end

function rescale_internal_mode(a)
    if -π / 2 < a < π / 2
        return csc(a)
    else
        return sign(a)
    end
end

function xspanKAK(q)
    a = q[1]
    return (-π / 2 - abs(a), π / 2 + abs(a))
end

function ηKAK(x, a, c)
    if a == 0.0
        if -π / 2 < x < π / 2
            return 4 * c * (cos(2 * x) + cos(4 * x))
        else
            return 0.0
        end
    else
        return kink(x + π / 2 + a) - kink(x + π / 2 - a) +
               (internal_mode(x + a) - internal_mode(x - a)) * c * rescale_internal_mode(a)
    end
end

function ηKAK(x, q)
    a, c = q
    return ηKAK(x, a, c)
end

function collective_coordinates(params::CCKinkAntikink)
    @unpack v, tmax, saveat, method, quadgk_kwargs, solver_kwargs = params

    q̇0 = [-v, 0.0]
    q0 = [π / 2, 0.0]

    prob = SecondOrderODEProblem(
        cceq!, q̇0, q0, (0.0, tmax), (quadratic, ηKAK, xspanKAK, quadgk_kwargs)
    )
    sol = solve(prob, method; saveat, save_idxs=[3, 4], solver_kwargs...)

    return Dict("t" => sol.t, "a" => sol[1, :], "c" => sol[2, :])
end
