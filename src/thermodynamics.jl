free_energies(chem_potentials::AbstractVector, stoich_matrix::AbstractMatrix) =
    vec(chem_potentials'*stoich_matrix)

function free_energies(chem_potentials::AbstractVector, reactions::Vector{Reaction})
    ∇r, ∇p = stoichiometric_matrix(ReactionNetwork(reactions))

    return free_energies(chem_potentials, ∇r - ∇p)
end

# compute reaction rates from:
# free energy,
# overall rate (sqrt(kf*kr)),
# inverse temperature
function reaction_rates(free_e, ovr_rate, beta)
    kf = ovr_rate.*exp.(-0.5*beta*free_e)
    kr = ovr_rate.*exp.( 0.5*beta*free_e)
    return kf,kr
end