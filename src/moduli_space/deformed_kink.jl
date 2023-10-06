@with_kw struct DeformedKinkModuliSpace{T<:Real}
    ϵ::T
    tmax::T = 10.0
    dt::T = 1e-2
end

function moduli_space(params::DeformedKinkModuliSpace)
    @unpack ϵ, tmax, dt = params

    function g_bb(b, c)
        return (
            2250 * π * (π^2 - 6) +
            128 * (225 * π^2 - 3152) * c +
            375 * π * (48 * π^2 - 275) * c^2
        ) / (54000 * b^3)
    end

    function g_bc(b, c)
        return (1024 - 375 * π * c) / (1200 * b^2)
    end

    function g_cc(b, c)
        return 5 * π / (8 * b)
    end

    function metric(q)
        b, c = q
        return [g_bb(b, c) g_bc(b, c); g_bc(b, c) g_cc(b, c)]
    end

    function potential(q)
        b, c = q
        return (256 * (b^2 - 1) * c + 60 * π * (1 + b^2) + 15 * π * c^2 * (32 * b^2 - 5)) /
               (240 * b)
    end

    q0 = [1 / (1 + ϵ); 0.0]
    p0 = [0.0; 0.0]
    tspan = (0.0, tmax)
    tsave = 0.0:dt:tmax

    prob = HamiltonianProblem(moduli_space_hamiltonian, p0, q0, tspan, (metric, potential))
    sol = solve(prob; saveat=tsave)
    return Dict("t" => tsave, "b" => sol[3, :], "c" => sol[4, :])

    return nothing
end
