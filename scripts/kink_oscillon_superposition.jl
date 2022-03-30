using Printf
using Compactons

mkpath("data/kink_oscillon_superposition")

for l ∈ 0.5:0.5:3.0
    for α ∈ 0.00:0.25:0.75
        filename = @sprintf "data/kink_oscillon_superposition/l=%.2f,alpha=%.2f.h5" l α
        print("Producing $filename … ")

        xsave, tsave, field, hamiltonian = kink_oscillon_superposition(l, α)
        save_data(filename, xsave, tsave, field, hamiltonian)
        println("done")
    end
end
