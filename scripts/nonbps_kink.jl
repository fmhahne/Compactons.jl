using Printf
using JLD2, CodecZlib
using Compactons

mkpath("data/nonbps_kink")

for ϵ ∈ -0.30:0.05:0.30
    filename = @sprintf "data/nonbps_kink/eps=%.2f.jld2" ϵ
    print("Producing $filename … ")

    x, t, η, H = nonbps_kink(ϵ)
    jldsave(filename, true; x, t, η, H)
    println("done")
end
