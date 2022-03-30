using Printf
using Compactons

mkpath("data/kink_antikink_scattering")

for V ∈ 0.00:1e-2:0.99
    filename = @sprintf "data/kink_antikink_scattering/V=%.2f.h5" V
    print("Producing $filename … ")

    xsave, tsave, field, hamiltonian = kink_antikink_scattering(V)
    save_data(filename, xsave, tsave, field, hamiltonian)
    println("done")
end
