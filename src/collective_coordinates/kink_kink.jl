# Kink-kink, non-relativistic model

@with_kw struct CCKinkKinkNonRel{T<:Real}
    v::T
    tmax::T = 10.0
    saveat::T = 0.01
    method = Tsit5()
    solver_kwargs = (;)
    quadgk_kwargs = (; atol=1e-6, order=20)
end

function xspanKK_non_rel(q)
    a = q[1]
    return (-π / 2 - abs(a), π / 2 + abs(a))
end

function ηKK_non_rel(x, q)
    a = q[1]
    return -2.0 + kink(x + a + π / 2) + kink(x - a + π / 2)
end

function collective_coordinates(params::CCKinkKinkNonRel)
    @unpack v, tmax, saveat, method, quadgk_kwargs, solver_kwargs = params

    q̇0 = [-v]
    q0 = [π / 2]

    prob = SecondOrderODEProblem(
        cceq!,
        q̇0,
        q0,
        (0.0, tmax),
        (quadratic, ηKK_non_rel, xspanKK_non_rel, quadgk_kwargs),
    )
    sol = solve(prob, method; saveat, save_idxs=[2], solver_kwargs...)

    return Dict("t" => sol.t, "a" => sol[1, :])
end

# Kink-kink, non-relativistic model with internal mode

@with_kw struct CCKinkKinkNonRelMode{T<:Real}
    v::T
    tmax::T = 10.0
    saveat::T = 0.01
    method = Tsit5()
    solver_kwargs = (;)
    quadgk_kwargs = (; atol=1e-6, order=20)
end

function ηKK_non_rel_mode(x, a, c)
    return -2.0 +
           kink(x + a + π / 2) +
           kink(x - a + π / 2) +
           c * (internal_mode(x + a) + internal_mode(x - a))
end

function ηKK_non_rel_mode(x, q)
    return ηKK_non_rel_mode(x, q...)
end

function collective_coordinates(params::CCKinkKinkNonRelMode)
    @unpack v, tmax, saveat, method, quadgk_kwargs, solver_kwargs = params

    q̇0 = [-v, 0.0]
    q0 = [π / 2, 0.0]

    prob = SecondOrderODEProblem(
        cceq!,
        q̇0,
        q0,
        (0.0, tmax),
        (quadratic, ηKK_non_rel_mode, xspanKK_non_rel, quadgk_kwargs),
    )
    sol = solve(prob, method; saveat, save_idxs=[3, 4], solver_kwargs...)

    return Dict("t" => sol.t, "a" => sol[1, :], "c" => sol[2, :])
end

# Kink-kink, relativistic model

@with_kw struct CCKinkKinkRel{T<:Real}
    v::T
    tmax::T = 10.0
    saveat::T = 0.01
    method = Tsit5()
    solver_kwargs = (;)
    quadgk_kwargs = (; atol=1e-6, order=20)
end

function xspanKK_rel(q)
    a = q[1]
    b = q[2]
    return (-π / (2 * b) - abs(a), π / (2 * b) + abs(a))
end

function ηKK_rel(x, a, b)
    return -2.0 + kink(b * (x + a) + π / 2) + kink(b * (x - a) + π / 2)
end

function ηKK_rel(x, q)
    return ηKK_rel(x, q...)
end

function collective_coordinates(params::CCKinkKinkRel)
    @unpack v, tmax, saveat, method, quadgk_kwargs, solver_kwargs = params

    q̇0 = [-v, 0.0]
    q0 = [π / 2, γ(v)]

    prob = SecondOrderODEProblem(
        cceq!, q̇0, q0, (0.0, tmax), (quadratic, ηKK_rel, xspanKK_rel, quadgk_kwargs)
    )
    sol = solve(prob, method; saveat, save_idxs=[3, 4], solver_kwargs...)

    return Dict("t" => sol.t, "a" => sol[1, :], "b" => sol[2, :])
end

# Kink-kink, relativistic model with internal mode

@with_kw struct CCKinkKinkRelMode{T<:Real}
    v::T
    tmax::T = 10.0
    saveat::T = 0.01
    method = Tsit5()
    solver_kwargs = (;)
    quadgk_kwargs = (; atol=1e-6, order=20)
end

function ηKK_rel_mode(x, a, b, c)
    return -2.0 +
           kink(b * (x + a) + π / 2) +
           kink(b * (x - a) + π / 2) +
           c * (internal_mode(b * (x + a)) + internal_mode(b * (x - a)))
end

function ηKK_rel_mode(x, q)
    return ηKK_rel_mode(x, q...)
end

function collective_coordinates(params::CCKinkKinkRelMode)
    @unpack v, tmax, saveat, method, quadgk_kwargs, solver_kwargs = params

    q̇0 = [-v, 0.0, 0.0]
    q0 = [π / 2, γ(v), 0.0]

    prob = SecondOrderODEProblem(
        cceq!, q̇0, q0, (0.0, tmax), (quadratic, ηKK_rel_mode, xspanKK_rel, quadgk_kwargs)
    )
    sol = solve(prob, method; saveat, save_idxs=[4, 5, 6], solver_kwargs...)

    return Dict("t" => sol.t, "a" => sol[1, :], "b" => sol[2, :], "c" => sol[3, :])
end
