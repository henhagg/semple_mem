function logpdf_prior_c(ci::Vector{T}, μ::Vector{T}, τ::Vector{T})::T where T<:AbstractFloat
    _logpdf = 0.0
    for i in eachindex(ci)
        _logpdf += logpdf(Normal(μ[i], sqrt(1.0/τ[i])), ci[i])
    end
    return _logpdf
end

function sample_μ(τ::Vector{T}, μ₀::Vector{T}, c::Vector{Vector{T}})::Vector{T} where T<:AbstractFloat
    n = length(c)
    out = similar(τ)
    for i in eachindex(out)
        ci = [_c[i] for _c in c]
        _mean = n*τ[i] / (n*τ[i] + τ[i]) * mean(ci) + τ[i]/(n*τ[i] + τ[i])*μ₀[i]
        _precision = n*τ[i] + τ[i]
        out[i] = rand(Normal(_mean, sqrt(1.0 / _precision)))
    end
    return out
end

function sample_τ(c::Vector{Vector{T}}, α::Vector{T}, β::Vector{T},
                  μ₀::Vector{T})::Vector{T} where T<:AbstractFloat
    n = length(c)
    out = similar(α)
    for i in eachindex(out)
        ci = [_c[i] for _c in c]
        _α = α[i] + 0.5*n
        _β = β[i] + 0.5*sum((ci .- mean(ci)).^2) + n/(2.0*(n+1))*((mean(ci) - μ₀[i])^2)
        # _β = β[i] + 0.5*sum((ci .- mean(ci)).^2) + n/(2.0*n)*(mean(ci) - μ₀[i])
        
        # To reduce debugging efforts, everyone should just agree on how to specify
        # Γ-distribution
        out[i] = rand(Gamma(_α, 1.0 / _β))
    end
    return out
end
