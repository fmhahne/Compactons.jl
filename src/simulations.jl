function producedata(model, ∂ₜφ₀, φ₀, tsave; dx, dt=dx / 10, sampling=10, callbacks=[])
    savedhamiltonian = SavedValues(Float64, Vector{Float64})
    cbhamiltonian = SavingCallback(gethamiltonian, savedhamiltonian; saveat=tsave)

    callback = CallbackSet(cbhamiltonian, callbacks...)

    N = length(φ₀)
    tspan = (tsave[begin], tsave[end])

    prob = SecondOrderODEProblem(fieldeq!, ∂ₜφ₀, φ₀, tspan, (model, dx))
    sol = solve(
        prob,
        RK4();
        adaptive=false,
        dt=dt,
        saveat=tsave,
        save_idxs=(N + 1):sampling:(2N),
        callback=callback,
    )

    φ = reduce(hcat, sol.u)
    H = reduce(hcat, savedhamiltonian.saveval)

    return φ, H
end

simulation(parameters) = error("Simulation $(typeof(parameters)) not implemented")
