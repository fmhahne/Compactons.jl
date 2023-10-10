function produce_data(model, ∂ₜϕ₀, ϕ₀, tsave; dx, dt=dx / 10, sampling=10, callbacks=[])
    saved_hamiltonian = SavedValues(Float64, Vector{Float64})
    cb_hamiltonian = SavingCallback(get_hamiltonian, saved_hamiltonian; saveat=tsave)

    callback = CallbackSet(cb_hamiltonian, callbacks...)

    N = length(ϕ₀)
    tspan = (tsave[begin], tsave[end])

    prob = SecondOrderODEProblem(field_equation!, ∂ₜϕ₀, ϕ₀, tspan, (model, dx))
    sol = solve(
        prob,
        RK4();
        adaptive=false,
        callback,
        dt,
        saveat=tsave,
        save_idxs=(N + 1):sampling:(2N),
    )

    ϕ = reduce(hcat, sol.u)
    H = reduce(hcat, saved_hamiltonian.saveval)

    return ϕ, H
end

simulation(params) = error("Simulation $(typeof(params)) not implemented")
