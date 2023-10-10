@with_kw struct KinkAntikinkRelModuliSpace{T<:Real}
    v::T
    tmax::T = 10.0
    saveat::T = 1e-2
end

function moduli_space(params::KinkAntikinkRelModuliSpace)
    @unpack v, tmax, saveat = params

    function g_aa(a, b)
        if abs(a) >= π / (2 * b)
            return π * b
        else
            return b * (sin(2 * abs(a) * b) + (π - 2 * abs(a) * b) * cos(2 * a * b) + π)
        end
    end

    function g_ab(a, b)
        if abs(a) >= π / (2 * b)
            return 0.0
        else
            return a * (sin(2 * abs(a) * b) + (π - 2 * abs(a) * b) * cos(2 * a * b))
        end
    end

    function g_bb(a, b)
        if abs(a) >= π / (2 * b)
            return π * (π^2 - 6) / (12 * b^3)
        else
            return (
                -(16 * b^3 * abs(a)^3 - 6(π^2 - 2) * b * abs(a) + π * (π^2 - 6)) *
                cos(2 * b * a) - 3 * (-4 * π * b * abs(a) + π^2 - 2) * sin(2 * b * abs(a)) +
                π * (π^2 - 6)
            ) / (12 * b^3)
        end
    end

    function metric(q)
        a, b = q
        return [g_aa(a, b) g_ab(a, b); g_ab(a, b) g_bb(a, b)]
    end

    function potential(q)
        a, b = q
        if abs(a) >= π / (2 * b)
            return π * (b^2 + 1) / (2 * b)
        else
            return -(
                (b^2 - 3) * sin(2 * abs(a) * b) +
                (b^2 - 1) * (π - 2 * b * abs(a)) * cos(2 * a * b) - 4 * b * abs(a) -
                π * b^2 + π
            ) / (2 * b)
        end
    end

    q0 = [π / (2 * γ(v)); γ(v)]
    p0 = metric(q0) * [-v; 0.0]
    tspan = (0.0, tmax)
    tsave = 0.0:saveat:tmax

    prob = HamiltonianProblem(moduli_space_hamiltonian, p0, q0, tspan, (metric, potential))
    sol = solve(prob; saveat)
    return Dict("t" => tsave, "a" => sol[3, :], "b" => sol[4, :])
end
