@with_kw struct KinkAntikinkNonRelModuliSpace{T<:Real}
    v::T
    tmax::T = 10.0
    dt::T = 1e-2
end

function moduli_space(params::KinkAntikinkNonRelModuliSpace)
    @unpack v, tmax, dt = params

    function metric(a)
        if abs(a) >= π / 2
            return π
        else
            return π + sin(2 * abs(a)) + (π - 2 * abs(a)) * cos(2 * a)
        end
    end

    function potential(a)
        if abs(a) >= π / 2
            return π
        else
            return 2 * abs(a) + sin(2 * abs(a))
        end
    end

    q0 = π / 2
    p0 = -metric(q0) * v
    tspan = (0.0, tmax)
    tsave = 0.0:dt:tmax

    prob = HamiltonianProblem(moduli_space_hamiltonian, p0, q0, tspan, (metric, potential))
    sol = solve(prob; saveat=tsave)
    return Dict("t" => tsave, "a" => sol[2, :])
end
