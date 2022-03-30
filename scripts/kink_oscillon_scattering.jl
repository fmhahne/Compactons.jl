using Printf
using Compactons

mkpath("data/kink_oscillon_scattering")

for l ∈ 0.5:0.5:3.0
    for V ∈ 0.0:1e-1:0.9
        for α ∈ 0.00:0.25:0.75
            filename = @sprintf "data/kink_oscillon_scattering/l=%.2f,V=%.2f,alpha=%.2f.h5" l V α
            print("Producing $filename … ")

            xsave, tsave, field, hamiltonian = kink_oscillon_scattering(l, V, α)
            save_data(filename, xsave, tsave, field, hamiltonian)
            println("done")
        end
    end
end
