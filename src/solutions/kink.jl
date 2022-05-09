# Kink
export kink, ∂ₜkink, ∂ₓkink

"""
    kink(x)

Field ``η(x)`` for a kink at rest in the quadratic periodic model.
"""
function kink(x)
    if x ≤ 0
        0.0
    elseif x ≤ π
        1 - cos(x)
    else
        2.0
    end
end

"""
    kink(t, x)

Field ``η(t, x)`` for a kink at rest in the quadratic periodic model.
"""
kink(t, x) = kink(x)

"""
    kink(t, x)

Partial derivative ``∂ₜη(t, x)`` for a kink at rest in the quadratic periodic model.
"""
∂ₜkink(t, x) = 0.0

"""
    ∂ₓkink(x)

Derivative ``η'(x)`` for a kink at rest in the quadratic periodic model.
"""
function ∂ₓkink(x)
    if 0 ≤ x ≤ π
        sin(x)
    else
        0.0
    end
end

"""
    ∂ₓkink(t, x)

Partial derivative ``∂ₓη(t, x)`` for a kink at rest in the quadratic periodic model.
"""
∂ₓkink(t, x) = ∂ₓkink(x)

# Moving kink

"""
    kink(t, x, V)

Field ``η(t, x)`` for a kink moving uniformly with velocity ``V`` in the quadratic periodic
model.
"""
function kink(t, x, V)
    _, x′ = boost(t, x, V)

    kink(x′)
end

"""
    ∂ₜkink(t, x, V)

Partial derivative ``∂ₜη(t, x)`` for a kink moving uniformly with velocity ``V`` in the
quadratic periodic model.
"""
function ∂ₜkink(t, x, V)
    _, x′ = boost(t, x, V)

    -V * γ(V) * ∂ₓkink(x′)
end

"""
    ∂ₓkink(t, x, V)

Partial derivative ``∂ₓη(t, x)`` for a kink moving uniformly with velocity ``V`` in the
quadratic periodic model.
"""
function ∂ₓkink(t, x, V)
    _, x′ = boost(t, x, V)

    γ(V) * ∂ₓkink(x′)
end
