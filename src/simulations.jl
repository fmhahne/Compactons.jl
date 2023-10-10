function produce_data(model, ∂ₜϕ₀, ϕ₀, tsave; dx, dt=dx / 10, sampling=10, callbacks=[])
    savedhamiltonian = SavedValues(Float64, Vector{Float64})
    cbhamiltonian = SavingCallback(gethamiltonian, savedhamiltonian; saveat=tsave)

    callback = CallbackSet(cbhamiltonian, callbacks...)

    N = length(ϕ₀)
    tspan = (tsave[begin], tsave[end])

    prob = SecondOrderODEProblem(field_equation!, ∂ₜϕ₀, ϕ₀, tspan, (model, dx))
    sol = solve(
        prob,
        RK4();
        adaptive=false,
        dt=dt,
        saveat=tsave,
        save_idxs=(N + 1):sampling:(2N),
        callback=callback,
    )

    ϕ = reduce(hcat, sol.u)
    H = reduce(hcat, savedhamiltonian.saveval)

    return ϕ, H
end

simulation(parameters) = error("Simulation $(typeof(parameters)) not implemented")
