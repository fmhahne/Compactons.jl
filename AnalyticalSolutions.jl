# Lorentz factor

γ(V) = 1.0 / √(1 - V^2)

# Oscillon

x_L(t, v₀, l) = v₀ * t
x_R(t, v₀, l) = x_L(t, v₀, l) + l

function auxF(x, v₀, l)
    if 0 <= x && x <= (1 + v₀) * l / 2
        -0.25 / (1.0 + v₀) * x^2
    elseif (1 + v₀) * l / 2 < x && x <= l
        -l^2 / 8.0 + 0.25 / (1 - v₀) * (l - x)^2
    else
        0.0
    end
end

function auxf(x, v₀, l)
    if 0 <= x && x <= (1 + v₀) * l / 2
        -x / (1 + v₀)
    elseif (1 + v₀) * l / 2 < x && x <= l
        -(l - x) / (1 - v₀)
    else
        0.0
    end
end

τ(t, l) = abs(mod(t - l / 2, l) - l / 2)
σ(t, l) = sign(mod(t - l / 2, l) - l / 2)

function oscillon(t, x; v₀=0.0, l=1.0)
    τₗ = τ(t, l)
    σₗ = σ(t, l)

    if (x_L(τₗ, v₀, l) <= x && x <= τₗ)
        σₗ * (auxF(x + τₗ, v₀, l) - auxF(x - τₗ + l, v₀, l) + τₗ^2 / 2.0 - l^2 / 8.0)
    elseif (τₗ < x && x <= l - τₗ)
        σₗ * (auxF(x + τₗ, v₀, l) - auxF(x - τₗ, v₀, l) + τₗ^2 / 2.0)
    elseif (l - τₗ < x && x <= x_R(τₗ, v₀, l))
        σₗ * (auxF(x + τₗ - l, v₀, l) - auxF(x - τₗ, v₀, l) + τₗ^2 / 2.0 - l^2 / 8.0)
    else
        0.0
    end
end

function oscillon_t(t, x; v₀=0.0, l=1.0)
    τₗ = τ(t, l)

    if x_L(τₗ, v₀, l) <= x && x <= τₗ
        0.5 * auxf(x + τₗ, v₀, l) + 0.5 * auxf(x - τₗ + l, v₀, l) + τₗ
    elseif τₗ < x && x <= l - τₗ
        0.5 * auxf(x + τₗ, v₀, l) + 0.5 * auxf(x - τₗ, v₀, l) + τₗ
    elseif l - τₗ < x && x <= x_R(τₗ, v₀, l)
        0.5 * auxf(x + τₗ - l, v₀, l) + 0.5 * auxf(x - τₗ, v₀, l) + τₗ
    else
        0.0
    end
end

function oscillon_x(t, x; v₀=0.0, l=1.0)
    τₗ = τ(t, l)
    σₗ = σ(t, l)

    if x_L(τₗ, v₀, l) <= x && x <= τₗ
        σₗ * 0.5 * (auxf(x + τₗ, v₀, l) - auxf(x - τₗ + l, v₀, l))
    elseif τₗ < x && x <= l - τₗ
        σₗ * 0.5 * (auxf(x + τₗ, v₀, l) - auxf(x - τₗ, v₀, l))
    elseif l - τₗ < x && x <= x_R(τₗ, v₀, l)
        σₗ * 0.5 * (auxf(x + τₗ - l, v₀, l) - auxf(x - τₗ, v₀, l))
    else
        0.0
    end
end

function moving_oscillon(t, x, V; v₀=0.0, l=1.0)
    t_prime = γ(V) * (t - V * x)
    x_prime = γ(V) * (x - V * t)

    oscillon(t_prime, x_prime, v₀=v₀, l=l)
end

function moving_oscillon_t(t, x, V; v₀=0.0, l=1.0)
    t_prime = γ(V) * (t - V * x)
    x_prime = γ(V) * (x - V * t)

    γ(V) * (oscillon_t(t_prime, x_prime; v₀=v₀, l=l) - V * oscillon_x(t_prime, x_prime; v₀=v₀, l=l))
end

function moving_oscillon_x(t, x, V; v₀=0.0, l=1.0)
    t_prime = γ(V) * (t - V * x)
    x_prime = γ(V) * (x - V * t)

    γ(V) * (oscillon_x(t_prime, x_prime; v₀=v₀, l=l) - V * oscillon_t(t_prime, x_prime; v₀=v₀, l=l))
end

# Kink

function kink(x)
    if x <= 0
        return 0.0
    elseif x <= π
        return 1 - cos(x)
    else
        return 2.0
    end
end

function kink_x(x)
    if x >= 0 && x <= pi
        return sin(x)
    else
        return 0.0
    end
end


function moving_kink(t, x, V)
    γ(V) = 1.0 / sqrt(1 - V^2)
    x_prime = γ(V) * (x - V * t)

    return kink(x_prime)
end

function moving_kink_t(t, x, V)
    x_prime = γ(V) * (x - V * t)

    return -V * γ(V) * kink_x(x_prime)
end


function moving_kink_x(t, x, V)
    x_prime = γ(V) * (x - V * t)

    return γ(V) * kink_x(x_prime)
end
