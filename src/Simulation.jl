using Random, DataFrames, Distributions

"""
    simu_u_ped(nid, ng)
---

# Description
Simulate a simple population of 0:`ng` generations:

- Each generation has `nid` breeders
- Half of them are males.
- Random mating
- Each pair has two offspring.

# Parameters
- `nid`: even (or, +1), positive (or, error)
- `ng`: positive integer

# Results
Return a two column matrix, which are `pa` and `ma`.
Row number = ID number.
Unknown parent = 0.

"""
function simu_u_ped(nid::Int, ng::Int)
    # subtitle("Simulate a pedigree of nid x ng = ($nid x $ng)")
    nid < 1 && @error "nid (= $nid) should > 0"
    ng  < 1 && @error "ng (= $ng) should > 0"
    nid % 2 == 1 && begin
        nid += 1
        @warn "nid increased by 1, = $nid now"
    end

    # generation 0
    ped = zeros(Int, 2nid)

    # generation 1:ng
    @inbounds for g in 1:ng
        fra = (g - 1) * nid + 1
        id = shuffle(fra : fra + nid - 1)
        for i in 1:2:nid
            append!(ped, [id[i], id[i+1], id[i], id[i+1]])
        end
    end
    #done()
    reshape(ped, 2, :)'
end


"""
    simu_p_bvy(ped, h2; μ = .0, Vp = 1.0)
---
Simulate breeding values (`bv`) and phenotypes (`y`) with polygenic method.
Given a pedigree `ped`, heritability `h2`, this function return two vectors, breeding values and phenotypes.

The population mean μ is 0.0, variance of phenotypes is 1.0 by default.
"""
function simu_p_bvy(ped, h2; μ = .0, Vp = 1.0)
    # subtitle("Simulate breeding and phenotype values with h2 = $h2, vp = $vp")
    ! (0.0 ≤ h2 ≤ 1.0) && throw(DomainError("h2 should be in [0, 1]"))
    Vp ≤ .0 && throw(DomainError("vp should be in (0, ∞)"))

    Va = Vp * h2
    RndAdd = Normal(.0, sqrt(Va/2.)) # when one parent is unknown
    Mendel = Normal(.0, sqrt(Va)/2.) # when one parent is known, add some Mendelian sampling error
    sde = sqrt(Vp * (1. - h2))
    
    nid = size(ped)[1]
    y   = rand(Normal(.0, sde), nid) # residuals first for phenotypes
    bv  = zeros(nid)
    @inbounds for i in 1:nid
        pa, ma = ped[i, :]
        bv[i] += pa > 0 ? bv[pa] / 2. + rand(Mendel) : rand(RndAdd)
        bv[i] += ma > 0 ? bv[ma] / 2. + rand(Mendel) : rand(RndAdd)
    end
    y = y + bv .+ μ
    # done()
    return bv, y 
end


"""
    simu_q_bvy(ped, h2; nqtl = 1000, keep = false, shape = 0.25, vp = 1.0, MAF = 0.01)
---
Simulate breeding values (`bv`) and phenotypes (`y`) with the QTL way.
This also a dirty simulation, which means there are some unrealistic hypotheses.
This simulation is however more realistic than `simu_p_bvy`.

1. Sample the biallelic QTL allele frequencies from `Beta(0.25)` in `[MAF, 1-MAF]`.
2. Sample QTL alleles from Uniform[0, 1] distribution on above frequencies.
3. Sample and adjust QTL effects.
4. Drop QTL alleles into the pedigree.
5. Calculate breeding values and add noises to create phenotypes.

By default,
1. This program uses 1000 QTL
2. QTL allele frequencies follow Beta distribution with a default shape (0.25, 0.25).
3. MAF = 0.01
4. tol = 1e-3, to adjust mean and variance

All above 3 can be modified.  This subroutine trys to simulate a population of Vp = 1.0.
"""
function simu_q_bvy(ped, h2; nQTL = 1000, Shape = 0.25, MAF = 0.01, vp = 1.0, tol=1e-3)
    # sample frequencies for QTL
    freq = Float64[]            # pop frequency of (MAF, 1-MAF), exclusive
    while length(freq) < nQTL   # sizeof is equivalent to capacity of C++
        f = rand(Beta(Shape))
        MAF < f < 1-MAF && push!(freq, f)
    end

    nID  = size(ped)[1]
    dfrq = Bernoulli.(freq)     # for an unknown parent, sample alleles from population
    allele = zeros(Bool, (2nID, nQTL)) # Only one byte each allele, saving memory
    
    for i in 1:2:2nID           # Simulate the alleles
        id = Int((i+1)/2)
        pa, ma = ped[id, :]
        if pa == 0
            allele[i, :] = rand.(dfrq)
        else
            ix = CartesianIndex.(rand(2pa-1:2pa, nQTL), 1:nQTL)
            allele[i, :] = allele[ix]
        end
        if ma == 0
            allele[i+1, :] = rand.(dfrq)
        else
            ix = CartesianIndex.(rand(2ma-1:2ma, nQTL), 1:nQTL)
            allele[i+1, :] = allele[ix]
        end
    end

    gt = begin                  # convert to 012 genotypes of Float64
        tmp = Float64[]
        for i in 1:2:2nID
            append!(tmp, convert.(Float64, allele[i, :] + allele[i+1, :]))
        end
        reshape(tmp, nID, :)
    end

    dQTL = Normal(0, sqrt(h2/nQTL))
    eQTL = rand(dQTL, nQTL)
    bv = gt * eQTL
    y = rand(Normal(0, sqrt(1-h2)), nID) + bv
    bv, y
    #=
    # sample QTL genotypes in generation 0
    qg = begin
        tmp = Bool[]            # only one byte each element
        for i in 1:nQTL
            append!(tmp, rand(Bernoulli(freq[i]), 2nid))
        end
        reshape(tmp, 2nid, :)
    end

    # sample QTL effects
    gt = begin              # convert allele types to 012 genotypes
        tmp = Float64[]
        for i in 1:2:2nid
            append!(tmp, convert.(Float64, qg[i, :] .+ qg[i+1, :]))
        end
        reshape(tmp, nid, :)
    end
    step = .1/nQTL
    eqtl = randn(nQTL)
    mu, vbv = begin
        m = Float64[]
        tmp = Float64[]
        for f in step:step:0.1
            e = eqtl .* f
            bv = gt * e
            push!(tmp, var(bv))
            push!(m, mean(bv))
        end
        m, tmp
    end
    =#
end
