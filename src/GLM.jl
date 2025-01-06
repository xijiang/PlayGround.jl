#=
- Original data were in "~/R/x86_64-redhat-linux-gnu-library/4.0/GLMsData/data".
- I converted them to DataFrames and saved to this project.
=#
using RData, DataFrames, StatsPlots
"""
    function glmsdata()
---
The analysis of the `lungcap` data.
"""
function glmsdata()
    title("Lung capacity data")
    lungcap = load("/home/xijiang/R/x86_64-redhat-linux-gnu-library/4.0/GLMsData/data/lungcap.rda")["lungcap"]
    categorical!(lungcap, :Smoke, compress=true)
    categorical!(lungcap, :Gender, compress=true)
    println(first(lungcap, 5))
    describe(lungcap)
    @df lungcap scatter(:Age, :FEV, leg=false)
    @df lungcap boxplot(string.(:Gender), :FEV, leg=false)
    @df lungcap[lungcap.Smode.==1, :] scatter(:Age, :FEV, leg=false)

    title("Noisy miner data")
    data = load("/home/xijiang/R/x86_64-redhat-linux-gnu-library/4.0/GLMsData/data/nminer.rda")["nminer"]
end

    
