using CSV, JSON, DataFrames, Distributions, LinearAlgebra, ProgressMeter, Plots, StatsPlots
import Random

include(joinpath(@__DIR__, "kalman.jl"))
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "adaptive_step.jl"))
include(joinpath(@__DIR__, "population_parameters.jl"))

function read_data(path_data)
    dataobs = CSV.read(path_data, DataFrame)
    tobs = dataobs[:, "time"]
    dt = tobs[2] - tobs[1]
    data_matrix = select!(dataobs, Not([:time])) |> Matrix
    return data_matrix, dt
end

function run_OU(μ_prior, α_prior, β_prior, ξ_prior, path_data, niterations)
    μ_start = [-0.7, 2.3, -0.9]
    τ_start = [4.0, 10.0, 4.0]
    c_start = [rand.(Normal.(μ_start, sqrt.(1.0 ./ τ_start))) for i in 1:40]
    ξ_start = log(0.3)
    compute_prior_ξ = (x) -> logpdf(ξ_prior, x)[1]

    μ_chain = Matrix{Float64}(undef, (3, niterations))
    τ_chain = Matrix{Float64}(undef, (3, niterations))
    μ_chain[:, 1] .= μ_start
    τ_chain[:, 1] .= τ_start

    dataobs, dt = read_data(path_data)
    nind = size(dataobs)[2]
    compute_ll_c_ξ = get_ll_func_c_ξ(dataobs, dt)

    # Initialise adaptive samplers for ξ and c
    adaptive_sampler_c = Vector{Dict}(undef, nind)
    _compute_prior_c = (c) -> logpdf_prior_c(c, μ_start, τ_start)
    for i in eachindex(adaptive_sampler_c)
        _compute_ll = get_ll_c_func(ξ_start, compute_ll_c_ξ[i])
        adaptive_sampler_c[i] = init_mcmc_stepper(c_start[i], niterations, _compute_ll,
                                                _compute_prior_c)
    end
    _compute_ll_vec, _compute_ll_sum = get_ll_ξ_func(c_start, compute_ll_c_ξ)
    adaptive_sampler_ξ = init_mcmc_stepper(ξ_start, niterations, _compute_ll_sum, compute_prior_ξ)

    # Main MCMC loop
    @showprogress "Running MCMC" for iteration in 2:niterations
        # Individial parameter c
        for j in 1:nind
            _compute_prior_c = (c) -> logpdf_prior_c(c, μ_chain[:, iteration-1], τ_chain[:, iteration-1])
            _compute_ll = get_ll_c_func(adaptive_sampler_ξ[:chain][iteration-1], compute_ll_c_ξ[j])
            take_mcmc_step!(adaptive_sampler_c[j], _compute_ll, _compute_prior_c, iteration)
        end

        # Population parmaeter ξ
        c_current = [as[:chain][:, iteration] for as in adaptive_sampler_c]
        _compute_ll_vec, _compute_ll_sum = get_ll_ξ_func(c_current, compute_ll_c_ξ)
        take_mcmc_step!(adaptive_sampler_ξ, _compute_ll_sum, compute_prior_ξ, iteration)
        # Need to update old ll for individual parameters if accepting

        # Conjugate priors for population parameters
        τ_chain[:, iteration] .= sample_τ(c_current, α_prior, β_prior, μ_prior)
        μ_chain[:, iteration] .= sample_μ(τ_chain[:, iteration], μ_prior, c_current)
    end

    return Dict(:μ => μ_chain, :τ => τ_chain, :ξ => adaptive_sampler_ξ, :c => adaptive_sampler_c)
end

random_seed = 123
Random.seed!(random_seed)
path_data = joinpath(@__DIR__, "..", "..", "ornstein_uhlenbeck_unperturbed_noise", "num_observation_40", "observations.csv")

ξ_prior_mean = -1.2
ξ_prior_sd = 1.0
ξ_prior = Normal(ξ_prior_mean, ξ_prior_sd)
α_prior = [6.0, 6.0, 6.0]
β_prior = [2.0, 1.0, 2.0]
μ_prior = [0.0, 1.5, 0.0]
niterations = 200_000

res = run_OU(μ_prior, α_prior, β_prior, ξ_prior, path_data, niterations)
a = 1


# PLOT RESULTS 
burnin = 1
density(res[:ξ][:chain][burnin:end])
vline!([log(0.3)])

plot(res[:ξ][:chain][burnin:end])
res[:μ]

density(res[:τ][1, burnin:end], title = "τ₁")
vline!([4.0])
density(res[:τ][2, burnin:end], title = "τ₂")
vline!([10.0])
density(res[:τ][3, burnin:end], title = "τ₃")
vline!([4.0])   

density(res[:μ][1, burnin:end], title = "μ₁")
vline!([-0.7])
density(res[:μ][2, burnin:end], title = "μ₂")
vline!([2.3])
density(res[:μ][3, burnin:end], title = "μ₃")
vline!([-0.9])

p = density(res[:c][1][:chain][3, :])
for i in 2:40
    p = density!(res[:c][i][:chain][3, :])
end
p
vline!([-0.9])

plot(res[:c][2][:chain][3, 1:burnin])

# WRITE RESULTS TO FILE
eta_matrix = vcat(res[:μ][:,:], res[:τ][:,:], res[:ξ][:chain][:,:])
CSV.write(joinpath(@__DIR__, "..", "kalman_post_samples_full.csv"),Tables.table(transpose(eta_matrix)), header=["mu1", "mu2", "mu3", "tau1", "tau2", "tau3", "xi1"])

# WRITE SETTINGS TO FILE
settings = Dict("alpha_prior" => α_prior,
    "beta_prior" => β_prior,
    "mu_prior" => μ_prior,
    "xi_prior_mean" => ξ_prior_mean,
    "xi_prior_sd" => ξ_prior_sd,
    "niterations" => niterations,
    "random_seed" => random_seed,
    "path_data" => path_data)
open("settings.json", "w") do f
    JSON.print(f, settings, 1)
end