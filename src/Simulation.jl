using Random, Distributions, Roots

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
    reshape(ped, 2, :)'
end

"""
    simu_c_bvy(ped, h2)
---
Simulation with Cholesky decomposition of **A**.
"""
function simu_c_bvy(ped, h2)
    message("Hello, world")
end

"""
    simu_p_bvy(ped, h2; μ = .0, Vp = 1.0)
---
Simulate breeding values (`bv`) and phenotypes (`y`) with polygenic method.
Given a pedigree `ped`, heritability `h2`, this function return two vectors, breeding values and phenotypes.

The population mean μ is 0.0, variance of phenotypes is 1.0 by default.
"""
function simu_p_bvy(ped, h2; μ = .0, Vp = 1.0)
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
    return bv, y 
end


"""
    simu_q_founder(nID, h2; μ = .0, Vp = 1.0, nQTL = 1000, shape = 0.25, MAF = 0.01)
---
Simulate a founder population of size = `nid`.  This population will have
- Va = h2
= Vp = 1.0

unless specified otherwize.  This also assume the QTL
- are independent 
- frequencies follow a β(0.25) distribution
"""
function simu_q_founder(nID, h2; μ = .0, Vp = 1.0, nQTL = 1000, Shape = 0.25, MAF = 0.01)
    # sample frequencies for QTL
    freq = Float64[]            # pop frequency of (MAF, 1-MAF), inclusive
    while length(freq) < nQTL   # sizeof is equivalent to capacity of C++
        f = rand(Beta(Shape))
        MAF ≤ f ≤ 1-MAF && push!(freq, f)
    end

    allele = begin
        tmp = Bool[]
        for f in freq
            append!(tmp, rand(Bernoulli(f), 2nID))
        end
        reshape(tmp, 2nID, :)
    end

    eQTL = randn(nQTL)
    M = begin
        tmp = Float64[]
        for i in 1:2:2nID
            append!(tmp, convert.(Float64, allele[i, :] + allele[i+1, :]))
        end
        reshape(tmp, nQTL, :)
    end

    for _ in 1:3                # repeat 3 times to reach μ=0., var=h2*Vp
        fv(f) = begin
            tmp = eQTL .* f
            var(M'tmp) - h2*Vp
        end
        eQTL .*= fzero(fv, 1.0)
        
        fm(m) = begin
            tmp = eQTL .- m
            mean(M'tmp)
        end
        eQTL .-= fzero(fm, 0.)
    end
    allele, eQTL
end

"""
    simu_q_bvy(ped, sde, allele, eQTL)
---
1. Drop `allele` throw founders into `ped`.
2. Calculate `bv`
3. Add noise to `bv` → `y`
4. Return `bv`, `y`
"""
function simu_q_bvy(ped, sde, allele, eQTL)
    nFdr = size(allele)[1] ÷ 2
    nQTL = length(eQTL)
    nID  = size(ped)[1]
    y = randn(nID) .* sde

    @inbounds for id in nFdr+1:nID        # check pedigree
        if ped[id, 1] == 0 || ped[id, 2] == 0
            throw(DomainError("ID $id's parent(s) unknown"))
        elseif ped[id, 1] >= id || ped[id, 2] >= id
            throw(DomainError("ID $id's parents not before this ID"))
        end
    end
    
    mat = zeros(Bool, 2nID, nQTL)
    mat[1:2nFdr, :] = allele
    
    bv = begin
        # Drop alleles
        @inbounds for id in nFdr+1:nID
            pa, ma = ped[id, :]
            ix = CartesianIndex.(2pa .- rand(0:1, nQTL), 1:nQTL)
            mat[2id-1, :] = mat[ix]
            ix = CartesianIndex.(2ma .- rand(0:1, nQTL), 1:nQTL)
            mat[2id, :] = mat[ix]
        end
        M = begin
            tmp = Float64[]
            @inbounds for i in 1:2:2nID
                append!(tmp, convert.(Float64, mat[i, :] + mat[i+1, :]))
            end
            reshape(tmp, nQTL, :)
        end
        M'eQTL
    end
    y += bv
    bv, y
end
