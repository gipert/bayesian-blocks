using StatsBase#, Plots
#import GR

include("../bayesian_blocks.jl")

function test(range=(100, 5300), approx=false, binning=1)

    gerda = Float32[]
    gerda = readdlm("test.dat")[:,1]
    gerda = gerda[(gerda .> range[1]) .& (gerda .< range[2])]
    gerda = gerda[(gerda .< 2014) .| (gerda .> 2064)]

    if approx
        h = fit(Histogram, gerda, range[1]:binning:range[2], closed = :left)
    else
        h = normalize(fit(Histogram, gerda, range[1]:0.1:range[2], closed = :left))
    end

    edges = bayesian_blocks(approx ? h : gerda,
                            logfitness = :cash, logprior = :p0,
                            p0 = 0.05)

    hb = normalize(fit(Histogram, gerda, edges, closed = :left))

    approx && (h = normalize(h))

#    plot(
#         [h, hb],
#         st = :step,
#         w = [1 3],
#         linecolor = [:steelblue :orangered],
#         yscale = :log10,
#         lab = ["enriched detectors (1 keV)" "bayesian blocks (approx)"],
#         legend = true,
#         xlabel = "energy [keV]",
#         ylabel = "frequency"
#        )

#    savefig("test.pdf")
#    writedlm("edges.dat", hb.edges)

    return h, hb
end
