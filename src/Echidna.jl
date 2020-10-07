"""
- Simulate a squence of datasets.
- Do the analysis with Echidna.
"""
function Echidna(h2, vp)
    subtitle("Test Echidna data analysis")
    item("Simulate a pedigree")
    ped = simu_u_ped(10, 2)
    done()
    item("Simulate breeding and phenotype values")
    bv, pv = simu_bv_p(h2, vp, ped)
    done()
    item("writing to Echidna formats")
    done()
end
