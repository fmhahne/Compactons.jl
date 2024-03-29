x₀(k) = √(π / (2 - k)) * Γ(1 + 1 / (1 + k)) / Γ(1 / 2 + 1 / (1 + k))
function bps(η, x, k)
    return x₀(k) +
           (-1 + η) / √(2 - k) *
           _₂F₁(1 / 2, 1 / (1 + k), 1 + 1 / (1 + k), abs(1 - η)^(1 + k)) - x
end

"""
    generalized_kink(t, x; k)

Field ``η(t, x)`` for a kink at rest in the generalized model with parameter ``k``.
"""
function generalized_kink(t, x; k)
    if x <= 0
        0.0
    elseif x < 2 * x₀(k)
        find_zero(η -> bps(η, x, k), (0.0, 2.0))
    else
        2.0
    end
end
