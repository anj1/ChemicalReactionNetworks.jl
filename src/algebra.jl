# (Linear) algebraic functions for Chemical Reaction Nets.
# Compute conservation laws, cycles, emergent cycles.

using Nemo 

# generate stoichiometric matrices from list of reactions.
# matrices: (nspecies,nreactions)
# returns two matrices: reactants and products
function stoichiometric_matrix(n_species, reactions::Vector{Reaction})
    nr = length(reactions)

    ∇r = spzeros(Int, n_species, nr)
    ∇p = spzeros(Int, n_species, nr)

    for ri in 1:nr
        re = reactions[ri]
        for si in 1:length(re.reactants)
            s = re.reactants[si]
            c = re.stoichr[si]

            ∇r[s,ri] += c
        end

        for si in 1:length(re.products)
            s = re.products[si]
            c = re.stoichp[si]

            ∇p[s,ri] += c
        end
    end

    return ∇r,∇p
end

function stoichiometric_nullspace(n_species, reactions::Vector{Reaction}, tr)
    # We only want net stoichiometric matrix here.
    ∇r, ∇p = stoichiometric_matrix(n_species, reactions)
    ∇ = ∇r - ∇p

    ∇ = tr ? transpose(∇) : ∇

    S = MatrixSpace(ZZ, size(∇)...)

    ns,r = nullspace(S(full(∇)))
    m,n = size(ns)

    nsm = sparse(Int[ns[i,j] for i=1:m,j=1:n])

    return tr ? transpose(nsm) : nsm
end 

# returns dimensions of nullspace (dimension of kernel) without calculating
# the matrix explicitly. Used for 
function stoichiometric_dimker(n_species, reactions::Vector{Reaction}, tr)
    ∇r, ∇p = stoichiometric_matrix(n_species, reactions)
    ∇ = ∇r - ∇p

    ∇ = tr ? transpose(∇) : ∇

    return size(∇,2)-rank(full(∇))
end

conservation_laws(n_species, reactions::Vector{Reaction}) = 
    stoichiometric_nullspace(n_species, reactions, true)

n_conservation_laws(n_species, reactions::Vector{Reaction}) = 
    stoichiometric_dimker(n_species, reactions, true)

cycles(n_species, reactions::Vector{Reaction}) =
    stoichiometric_nullspace(n_species, reactions, false)

n_cycles(n_species, reactions::Vector{Reaction}) = 
    stoichiometric_dimker(n_species, reactions, false)

# Compute steady-states for both detailed-balanced 
# And complex-balanced states.

# The steady state concentration
# is one where the conserved moieties are distributed
# equally among all reactants.

# The equilibrium state is a steady-state where all
# concentration currents vanish.
# The solution is based on equation (18) in the paper:
#   kf/kr = prod_s z[s]^del_s,r for all r
# We take the logarithm, giving:
#   log(kf/kr) = sum del_s,r log z[s]
# And then we also append the conservation laws.
# This gives a linear system which can be solved.
function log_equilibrium_state(n_species, reactions::Vector{Reaction})
     # compute ratio of forward to backward reactions
    free_energy = log([r.kf/r.kr for r in reactions])

    dp,dr = stoichiometric_matrix(n_species,reactions)
    del = dr-dp

    logz = full(del')\free_energy

    return logz, del    
end 
function equilibrium_state(n_species, reactions::Vector{Reaction}) #, z0)
    logz, del = log_equilibrium_state(n_species, reactions)
    return exp(logz)
end

# The equilibrium state is often under-defined;
# this means we can add chemostats until it is fully defined,
# and get a unique solution.
# In addition, we can return the full solution space, which is:
# ss=nullspace(full(del'))
# Thus every solution can be written as:
# ss*x + logz0
# where x is some arbitrary vector and logz0=full(del')\log_k_ratio
function equilibrium_state_space(n_species, reactions::Vector{Reaction})
    logz, del = log_equilibrium_state(n_species, reactions)

    ss = nullspace(full(del'))
    return ss, logz
end

