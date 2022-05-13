"""
    γ(V)

Lorentz factor ``1 / √(1-V²)``.
"""
γ(V) = 1.0 / √(1 - V^2)

"""
    boost(t, x, V)

Return coordinates ``(t′, x′)`` in frame of reference moving uniformly with velocity ``V``.
"""
function boost(t, x, V)
    t′ = γ(V) * (t - V * x)
    x′ = γ(V) * (x - V * t)

    (t′, x′)
end
