import AdaptiveMCMC

function take_mcmc_step!(adaptive_sampler, compute_ll::Function, compute_prior, iteration::Integer)::Nothing
    r, s = adaptive_sampler[:r], adaptive_sampler[:s]
    ll, chain = adaptive_sampler[:ll], adaptive_sampler[:chain]
    AdaptiveMCMC.draw!(r, s)
    ll_new = compute_ll(r.y) + compute_prior(r.y)
    ll[iteration-1] = compute_ll(r.x) + compute_prior(r.x)
    α = min(one(ll[iteration-1]), exp(ll_new - ll[iteration-1]))
    if rand() <= α
        ll[iteration] = ll_new
        AdaptiveMCMC.accept!(r)
    else
        ll[iteration] = ll[iteration-1]
    end
    AdaptiveMCMC.adapt!(s, r, α, iteration)
    chain[:, iteration] .= r.x
    return  nothing
end

function init_mcmc_stepper(x0, niteration::Integer, compute_ll::Function,
                           compute_prior)::Dict
    ll = Vector{Float64}(undef, niteration)
    chain = Matrix{Float64}(undef, (length(x0), niteration))
    ll[1] = compute_ll(x0) + compute_prior(x0)
    chain[:, 1] .= x0
    # x0 must be a vector for AdaptiveMCMC
    x0 = length(x0) == 1 ? [x0] : x0
    r = AdaptiveMCMC.RWMState(x0)
    s = AdaptiveMCMC.RobustAdaptiveMetropolis(x0)

    return Dict(:r => r, :s => s, :ll => ll, :chain => chain)
end
