export γ, boost

γ(V) = 1.0 / √(1 - V^2)

function boost(t, x, V)
    t′ = γ(V) * (t - V * x)
    x′ = γ(V) * (x - V * t)

    (t′, x′)
end
