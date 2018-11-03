free_energies(chem_potentials::AbstractVector, stoich_matrix::AbstractMatrix) =
    vec(chem_potentials'*stoich_matrix)

function free_energies(chem_potentials::AbstractVector, reactions::Vector{Reaction})
    n_species = length(chem_potentials)

    ∇r, ∇p = stoichiometric_matrix(n_species, reactions)

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

# Note: we can replace a forcing rate k*exp(f(z))
# with an enzymatic complex, like so:
# A + B -> C + D  is transformed to:
#   A + B + E -> C1
#   C1 -> C + D + E,
#   E <-> ...  (very high rate)
# with the free energy of the complex and the enzyme production being
# chosen in a way such as to make the two reaction networks equivalent.