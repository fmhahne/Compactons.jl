# Oscillon at rest

"""
    x_L(t; l=1.0, v₀=0.0)

Position of the left border at time t ∈ [0, l/2] of a signum--Gordon oscillon at rest, with
size ``l`` and swaying velocity ``v₀``.
"""
x_L(t; l=1.0, v₀=0.0) = v₀ * t

"""
    x_L(t; l=1.0, v₀=0.0)

Position of the right border at time t ∈ [0, l/2] of a signum--Gordon oscillon at rest, with
size ``l`` and swaying velocity ``v₀``.
"""
x_R(t; l=1.0, v₀=0.0) = x_L(t; v₀, l) + l

function F(x; l=1.0, v₀=0.0)
    if 0 ≤ x ≤ (1 + v₀) * l / 2
        -0.25 / (1.0 + v₀) * x^2
    elseif (1 + v₀) * l / 2 < x ≤ l
        -l^2 / 8.0 + 0.25 / (1 - v₀) * (l - x)^2
    else
        0.0
    end
end

function f(x; l=1.0, v₀=0.0)
    if 0 ≤ x ≤ (1 + v₀) * l / 2
        -x / (1 + v₀)
    elseif (1 + v₀) * l / 2 < x ≤ l
        -(l - x) / (1 - v₀)
    else
        0.0
    end
end

τ(t, l) = abs(mod(t - l / 2, l) - l / 2)
σ(t, l) = sign(mod(t - l / 2, l) - l / 2)

"""
    oscillon(t, x; l=1.0, v₀=0.0)

Field ``ϕ(t,x)`` for a signum--Gordon oscillon at rest, with size ``l`` and swaying velocity
``v₀``.
"""
function oscillon(t, x; l=1.0, v₀=0.0)
    τₗ = τ(t, l)
    σₗ = σ(t, l)

    if (x_L(τₗ; v₀, l) ≤ x ≤ τₗ)
        σₗ * (F(x + τₗ; v₀, l) - F(x - τₗ + l; v₀, l) + τₗ^2 / 2.0 - l^2 / 8.0)
    elseif (τₗ < x ≤ l - τₗ)
        σₗ * (F(x + τₗ; v₀, l) - F(x - τₗ; v₀, l) + τₗ^2 / 2.0)
    elseif (l - τₗ < x ≤ x_R(τₗ; v₀, l))
        σₗ * (F(x + τₗ - l; v₀, l) - F(x - τₗ; v₀, l) + τₗ^2 / 2.0 - l^2 / 8.0)
    else
        0.0
    end
end

"""
    ∂ₜoscillon(t, x; l=1.0, v₀=0.0)

Partial derivative ``∂ₜϕ(t,x)`` for an signum--Gordon oscillon at rest, with size ``l`` and
swaying velocity ``v₀``.
"""
function ∂ₜoscillon(t, x; l=1.0, v₀=0.0)
    τₗ = τ(t, l)

    if x_L(τₗ; v₀, l) ≤ x ≤ τₗ
        0.5 * f(x + τₗ; v₀, l) + 0.5 * f(x - τₗ + l; v₀, l) + τₗ
    elseif τₗ < x ≤ l - τₗ
        0.5 * f(x + τₗ; v₀, l) + 0.5 * f(x - τₗ; v₀, l) + τₗ
    elseif l - τₗ < x ≤ x_R(τₗ; v₀, l)
        0.5 * f(x + τₗ - l; v₀, l) + 0.5 * f(x - τₗ; v₀, l) + τₗ
    else
        0.0
    end
end

"""
    ∂ₓoscillon(t, x; l=1.0, v₀=0.0)

Partial derivative ``∂ₜϕ(t,x)`` for an signum--Gordon oscillon at rest, with size ``l`` and
swaying velocity ``v₀``.
"""
function ∂ₓoscillon(t, x; l=1.0, v₀=0.0)
    τₗ = τ(t, l)
    σₗ = σ(t, l)

    if x_L(τₗ; v₀, l) ≤ x ≤ τₗ
        σₗ * 0.5 * (f(x + τₗ; v₀, l) - f(x - τₗ + l; v₀, l))
    elseif τₗ < x ≤ l - τₗ
        σₗ * 0.5 * (f(x + τₗ; v₀, l) - f(x - τₗ; v₀, l))
    elseif l - τₗ < x ≤ x_R(τₗ; v₀, l)
        σₗ * 0.5 * (f(x + τₗ - l; v₀, l) - f(x - τₗ; v₀, l))
    else
        0.0
    end
end

# Moving oscillon

"""
    oscillon(t, x, V; l=1.0, v₀=0.0)

Field ``ϕ(t,x)`` for a signum--Gordon oscillon moving uniformly with velocity ``V``. On its
frame of reference, the oscillon is of size ``l`` and swaying velocity ``v₀``.
"""
function oscillon(t, x, V; l=1.0, v₀=0.0)
    t′, x′ = boost(t, x, V)

    return oscillon(t′, x′; v₀, l)
end

"""
    ∂ₜoscillon(t, x, V; l=1.0, v₀=0.0)

Partial derivative ``∂ₜϕ(t,x)`` for a signum--Gordon oscillon moving uniformly with velocity
``V``. On its frame of reference, the oscillon is of size ``l`` and swaying velocity ``v₀``.
"""
function ∂ₜoscillon(t, x, V; l=1.0, v₀=0.0)
    t′, x′ = boost(t, x, V)

    return γ(V) * (∂ₜoscillon(t′, x′; v₀, l) - V * ∂ₓoscillon(t′, x′; v₀, l))
end

"""
    ∂ₓoscillon(t, x, V; l=1.0, v₀=0.0)

Partial derivative ``∂ₓϕ(t,x)`` for a signum--Gordon oscillon moving uniformly with velocity
``V``. On its frame of reference, the oscillon is of size ``l`` and swaying velocity ``v₀``.
"""
function ∂ₓoscillon(t, x, V; l=1.0, v₀=0.0)
    t′, x′ = boost(t, x, V)

    return γ(V) * (∂ₓoscillon(t′, x′; v₀, l) - V * ∂ₜoscillon(t′, x′; v₀, l))
end

function αₗ(V; v₀=0.0)
    if V ≥ 0
        (-1 + V * (2 + v₀)) / 2
    else
        (1 - V * (2 + v₀)) / 2
    end
end

function α₀(V; v₀=0.0)
    if V ≥ 0
        V
    else
        1 - V
    end
end

function αᵤ(V; v₀=0.0)
    if V ≥ 0
        (1 + V * (2 + v₀)) / 2
    else
        (3 - V * (2 + v₀)) / 2
    end
end

function αₛ(V; v₀=0.0)
    return (1 + v₀ * V) / 2
end

function ab(α, V; v₀=0.0)
    Vc = 1 / (2 + v₀)

    A = (1, 1)
    B = (-1, 0)
    C = (1, 0)
    D = (-1, 1)
    E = (1, -1)
    F = (-1, 2)

    if abs(V) == Vc
        if α ≤ α₀(V; v₀)
            V > 0 ? B : D
        else
            V > 0 ? C : E
        end
    elseif abs(V) < Vc
        if α ≤ α₀(V; v₀)
            V > 0 ? B : D
        elseif α ≤ αᵤ(V; v₀)
            V > 0 ? C : E
        else
            V > 0 ? D : F
        end
    else
        if α ≤ αₗ(V; v₀)
            V > 0 ? A : C
        elseif α ≤ α₀(V; v₀)
            V > 0 ? B : D
        else
            V > 0 ? C : E
        end
    end
end

function a′b′(α, V; v₀=0.0)
    if α ≤ αₛ(V; v₀)
        (1, 0)
    else
        (-1, 1)
    end
end

"""
    x_L(α, V; l=1.0, v₀=0.0)

Position of the left border of a signum--Gordon oscillon moving uniformly with velocity
``V`` when its phase is ``α``. On its frame of reference, the oscillon is of size ``l`` and
swaying velocity ``v₀``.
"""
function x_L(α, V; l=1.0, v₀=0.0)
    a′, b′ = a′b′(α, V; v₀)
    return l * γ(V) / (1 + a′ * v₀ * V) * ((1 - V^2) * (v₀ * b′) + α * (V + v₀ * a′))
end

"""
    x_R(α, V; l=1.0, v₀=0.0)

Position of the right border of a signum--Gordon oscillon moving uniformly with velocity
``V`` when its phase is ``α``. On its frame of reference, the oscillon is of size ``l`` and
swaying velocity ``v₀``.
"""
function x_R(α, V; l=1.0, v₀=0.0)
    a, b = ab(α, V; v₀)
    return l * γ(V) / (1 + a * v₀ * V) * ((1 - V^2) * (1 + v₀ * b) + α * (V + v₀ * a))
end

"""
    L(α, V; l=1.0, v₀=0.0)

Size of a signum--Gordon oscillon moving uniformly with velocity ``V`` when its phase is
``α``. On its frame of reference, the oscillon is of size ``l`` and swaying velocity ``v₀``.
"""
L(α, V; l=1.0, v₀=0.0) = x_R(α, V; l, v₀) - x_L(α, V; l, v₀)
