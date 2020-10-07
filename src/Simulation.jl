using Random, DataFrames

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
    subtitle("Simulate a pedigree of nid x ng = ($nid x $ng)")
    nid < 1 && @error "nid (= $nid) should > 0"
    ng  < 1 && @error "ng (= $ng) should > 0"
    nid % 2 == 1 && begin
        nid += 1
        @warn "nid increased by 1, = $nid now"
    end

    # generation 0
    ped = [PM(0, 0) for i in 1:nid]

    # generation 1:ng
    @inbounds for g in 1:ng
        fra = (g - 1) * nid + 1
        id = shuffle(fra : fra + nid - 1)
        for i in 1:2:nid
            push!(ped, PM(id[i], id[i+1]))
            push!(ped, PM(id[i], id[i+1]))
        end
    end
    done()
    return ped
end


"""
    simu_g_p(h2, vp, ped; μ = 0.)
---
Given a pedigree `ped`, heritability `h2`, and phenotype variance `vp`, this function return two
vectors, breeding values and phenotypes.

The population mean μ is 0.0 by default.
"""
function simu_bv_p(h2, vp, ped; μ = 0.)
    subtitle("Simulate breeding and phenotype values with h2 = $h2, vp = $vp")
    ! (0.0 ≤ h2 ≤ 1.0) && throw(DomainError("h2 should be in [0, 1]"))
    vp ≤ 0. && throw(DomainError("vp should be in (0, ∞)"))
    
    sda = sqrt(vp * h2 / 2)     # note, this is only for one parent
    sde = sqrt(vp * (1 - h2))
    mse = sqrt(vp * h2) / 2.    # Mendel sampling error SD
    
    nid = length(ped)
    y   = randn(nid) .* sde     # residuals first for phenotypes
    bv  = zeros(nid)
    @inbounds for i in 1:nid
        pa, ma = ped[i].pa, ped[i].ma
        if pa > 0
            bv[i] = bv[i] + bv[pa] / 2. + randn() * mse
        else
            bv[i] += randn() * sda
        end
        if ma > 0
            bv[i] = bv[i] + bv[ma] / 2. + randn() * mse
        else
            bv[i] += randn() * sda
        end
    end
    y = y + bv .+ μ
    done()
    return bv, y 
end


"""
    simu_q_pop(nid, ng, h2; nqtl = 1000, keep = false, shape = 0.7, MAF = 0.01)
---

This also a dirty simulation, which means there are some unrealistic hypotheses.
This simulation is however more realistic than `simu_u_ped` + `simu_g_p`.

1. Create a pedigree with `simu_u_ped`.
2. Sample the biallelic QTL allele frequencies from `Beta(0.7)` in `[MAF, 1-MAF]`.
3. Sample QTL alleles from Uniform[0, 1] distribution on above frequencies.
4. Sample and adjust QTL effects.
5. Drop QTL alleles into the pedigree.
6. Calculate breeding values and add noises to create phenotypes.

By defaut,
1. This program uses 1000 QTL
2. Returned dataframe doesn't include QTL types.
3. Beta distribution with a default shape (0.7, 0.7).
4. MAF = 0.01

All above 4 can be modified.
"""
function simu_q_pop(nid, ng, h2; nqtl = 1000, keep = false, shape = 0.7, MAF = 0.01)
    println("hello")
end
