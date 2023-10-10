"""
    toy_kink(t, x)

Field ``η(t, x)`` for a kink at rest in the toy model.
"""
function toy_kink(t, x)
    if x ≤ 0
        return 0.0
    elseif x ≤ √2
        return 0.5 * x^2
    elseif x ≤ 2 * √2
        return 2.0 - 0.5 * (x - 2 * √2)^2
    else
        return 2.0
    end
end
