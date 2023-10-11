function moduli_space_hamiltonian(p, q, (metric, potential), t)
    return dot(p, inv(metric(q)), p) / 2 + potential(q)
end
