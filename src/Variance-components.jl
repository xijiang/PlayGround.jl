using Roots, Statistics
"""
    sv_sde(vec)
---
Given a vector which assumes to follow N(μ, σ^2), this function returns numerical ML sde of the data.
"""
function sv_sde(vec)
    m = mean(vec)
    v = var(vec)
end

