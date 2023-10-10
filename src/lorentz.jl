"""
    γ(v)

Lorentz factor ``1 / √(1-v²)``.
"""
γ(v) = 1.0 / √(1 - v^2)

"""
    boost(t, x, v)

Return coordinates ``(t′, x′)`` in frame of reference moving uniformly with velocity ``v``.
"""
function boost(t, x, v)
    t′ = γ(v) * (t - v * x)
    x′ = γ(v) * (x - v * t)

    return (t′, x′)
end
