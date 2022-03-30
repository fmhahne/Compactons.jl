using Printf
using Compactons

mkpath("data/perturbed_kink")

for ϵ ∈ -0.50:0.01:0.50
    filename = @sprintf "data/perturbed_kink/eps=%.2f.h5" ϵ
    print("Producing $filename … ")

    xsave, tsave, field, hamiltonian = perturbed_kink(ϵ)
    save_data(filename, xsave, tsave, field, hamiltonian)
    println("done")
end
