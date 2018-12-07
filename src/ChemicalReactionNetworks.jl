module ChemicalReactionNetworks

export
    Reaction,
    ReactionNetwork,
    SpeciesComplex,
    ==,
    hash,
    isredundant,
    stoichiometric_matrix,
    mass_action,
    mass_action!,
    reaction_current,
    jacobian,
    from_stoichiometric_matrix,
    normal_form,
    conservation_laws,
    cycles,
    n_conservation_laws,
    n_cycles,
    rand_concentrations,
    chemostatted!,
    relabel!,
    complexes,
    complex_decomposition,
    deficiency,
    complex_concentration,
    equilibrium_state,
    equilibrium_state_space,
    cycle_affinities,
    specificity,
    complex_basis_reactions,
    free_energies,
    reaction_rates,
    petri_net

include("crn.jl")
include("algebra.jl")
include("sampling.jl")
include("chemostat.jl")
include("complexes.jl")
include("thermodynamics.jl")
include("petri.jl")

"""
A Julia package for studying and simulating Chemical
Reaction Networks, of the kind that would be encountered in
the study of complex metabolic networks, self-replication,
and so on.

API Overview:
-`Reaction(...)` creates a reaction.
-`ReactionNetwork(...)` creates a network of reactions.
-`mass_action(reactions, concentration)` returns change in concentration over time.
-`chemostatted!(reactions, chemostatted_species, chemostatted_concentrations)` chemostats certain species and modifies the network.   
-`equilibrium_state(reaction_net)` returns detailed-balance equilibrium state.
-`conservation_laws(reaction_net)` returns conservation laws of net
-`cycles(reaction_net)` returns cycles of net.
-`deficiency(reaction_net)` returns network deficiency.
-`complex_decomposition(reaction_net)` decomposes a CRN in terms of a set of complexes and incidence matrix.
"""

end