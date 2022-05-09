# Kink
export kink, ∂ₜkink, ∂ₓkink

function kink(x)
    if x ≤ 0
        0.0
    elseif x ≤ π
        1 - cos(x)
    else
        2.0
    end
end

kink(t, x) = kink(x)

∂ₜkink(t, x) = 0.0

function ∂ₓkink(x)
    if 0 ≤ x ≤ π
        sin(x)
    else
        0.0
    end
end

∂ₓkink(t, x) = ∂ₓkink(x)

# Moving kink

function kink(t, x, V)
    _, x′ = boost(t, x, V)

    kink(x′)
end

function ∂ₜkink(t, x, V)
    _, x′ = boost(t, x, V)

    -V * γ(V) * ∂ₓkink(x′)
end


function ∂ₓkink(t, x, V)
    _, x′ = boost(t, x, V)

    γ(V) * ∂ₓkink(x′)
end
