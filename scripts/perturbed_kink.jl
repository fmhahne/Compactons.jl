using Printf
using Compactons

mkpath("data/perturbed_kink")

for ϵ ∈ -0.50:0.01:0.50
    filename = @sprintf "data/perturbed_kink/eps=%.2f.h5" ϵ
    print("Producing $filename … ")

    x, t, η, H = perturbed_kink(ϵ)
    save_data(filename, x, t, η, H)
    println("done")
end
