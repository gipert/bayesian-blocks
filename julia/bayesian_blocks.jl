#= bayesian_blocks.jl
 =
 = Author: Luigi Pertoldi - pertoldi@pd.infn.it
 = Created: 29 Jun 2018
 =
 = The following software is distributed under the MIT licence
 = http://www.opensource.org/licenses/mit-license.php
 = Copyright (c) 2018 Luigi Pertoldi
 =
 = Use Julia >= v0.7!
 =#

using StatsBase, ProgressMeter

"""
    bayesian_blocks(data, logfitness=:cash, logprior=:p0, [gamma=0.01, p0=0.01])

Return an array of optimal change points for a set of one-dimensional data. This
is the implementation of the bayesian blocks algorithm, as outlined in [^1].

## Arguments
 * `data`: numeric array or a `StatsBase.Histogram`
 * `logfitness`: log of the block fitness function to be used, choose between
                 [:cash]
 * `logprior`: log of the prior distribution on the number of blocks to be used,
               choose between [:gamma, :p0]
 * `gamma`, `p0`...: set the parameter value for the specified prior distribution

## Example
```julia
using Distributions, StatsBase, Plots, LinearAlgebra

data = vcat(rand(Normal(0),1000),rand(Cauchy(5),1000))
data = data[(data .> -5) .& (data .< 10)]

h = fit(Histogram, gerda, -5:100:10, closed = :left)

# choose to use all data or an histogram of it!
hb = normalize(fit(Histogram, data, bayesian_blocks(data), closed = :left))
hb = normalize(fit(Histogram, data, bayesian_blocks(h),    closed = :left))

plot(data, st = :stephist, normalized = true, nbins=1000)
plot!(hb, st = :step, w = 3)
```

### Performance tips
You can convert your data container to a less precise representation to improve
the performance a bit, e.g.
```julia
x::Array{Float32} = [1.1, π, (√5-1)/2]
```

[^1]: Scargle, J et al. (2012) [https://arxiv.org/abs/1207.5578]
"""
function bayesian_blocks(x;
                         logfitness::Symbol=:cash, logprior::Symbol=:p0,
                         gamma=0.01, p0=0.01)

    if typeof(x) <: Array{<:Real,1}
        # take care of repeated data
        x_sorted = sort(x)
        x_unique = [x_sorted[1]]
        x_weight::Array{Int32} = [1]
        for i in 2:length(x_sorted)
            if x_sorted[i] == x_sorted[i-1]
                x_weight[end] += 1
            else
                push!(x_unique, x_sorted[i])
                push!(x_weight, 1)
            end
        end
    elseif typeof(x) <: Histogram{<:Integer,1}
        v = collect(x.edges[1])
        x_unique = [0.5 * (v[i] + v[i + 1]) for i = 1:length(v) - 1]
        x_weight = x.weights

        # delete empty bins
        deleteat!(x_unique, findall(iszero, x_weight))
        deleteat!(x_weight, findall(iszero, x_weight))
    else
        error("Unsupported input type: $(typeof(x))")
    end

    # final number of data points
    N = length(x_unique)

    # pre-defined (log)fitness functions
    logf_dict = Dict(
        :cash => (N_k, T_k) -> N_k * log(N_k/T_k)
    )

    # pre-defined (log)prior distributions on Nblocks
    logp_dict = Dict(
        # simple γ^Nblocks prior
        :gamma => (γ=gamma) -> log(γ),
        # Note that there was an error in this equation in the original Scargle
        # paper (the "log" was missing). The following corrected form is taken
        # from https://arxiv.org/abs/1304.2818
        :p0    => (Np=N, p=p0) -> log(73.53 * p * Np^(-0.478)) - 4
    )

    # check input
    !haskey(logf_dict, logfitness) && error("$logfitness function not defined!")
    !haskey(logp_dict, logprior)   && error("$logprior function not defined!")

    # save prior value for later computation
    ncp_prior = logp_dict[logprior]()

    # array of (all possible) block edges
    edges = vcat(x_unique[1],
                 0.5f0(x_unique[1:end-1]+x_unique[2:end]),
                 x_unique[end])

    # see Sec. 2.6 in [^1]
    best = Number[]
    last = Number[]

    # display progress bar for long computations
    # total number of steps: ∑n(n-1) = N(N^2-1)/3
    p = Progress(Integer(N*(N^2-1)/3), 2); m = 0

    for k in 1:N
        # define nice alias to mimic the notation used in [^1]
        F(r) = logf_dict[logfitness](cumsum(x_weight[r:k])[end], edges[k+1] - edges[r]) + ncp_prior

        # compute all possible configurations (Eq. (8) in [^1])
        A = [F(r) + (r == 1 ? 0 : best[r-1]) for r in 1:k]

        # save best configuration
        push!(last, argmax(A))
        push!(best, maximum(A))

        update!(p, m += k*(k-1))
    end

    # extract changepoints by iteratively peeling off the last block
    cp = Number[]
    i = N+1
    while i != 0
        push!(cp, i)
        i = (i == 1 ? 0 : last[i-1])
    end

    return [edges[j] for j in cp[end:-1:1]]
end
