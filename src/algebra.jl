# (Linear) algebraic functions for Chemical Reaction Nets.
# Compute conservation laws, cycles, emergent cycles.

using Nemo 
using SparseArrays 

full(a) = convert(Array,a)

# generate stoichiometric matrices from list of reactions.
# matrices: (nspecies,nreactions)
# returns two matrices: reactants and products
function stoichiometric_matrix(n_species::Integer, reactions::Vector{Reaction})
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

# generates a list of reactions from stoichiometric matrix
function from_stoichiometric_matrix(∇r::AbstractMatrix, ∇p::AbstractMatrix, kf::Vector{T}, kr::Vector{T}, names=[]) where T<:Real
    @assert size(∇r)==size(∇p)
    nreactions = size(∇r,2)

    reactions = Reaction[]
    for i = 1 : nreactions
        reactants,stoichr = findnz(∇r[:,i])
        products,stoichp  = findnz(∇p[:,i])
        nm = isempty(names) ? [] : names[i]
        push!(reactions, Reaction(reactants,stoichr,products,stoichp,kf[i],kr[i],nm))
    end

    return reactions
end

# return normal form of reaction net
# for example, S1 + S1 -> S2 becomes 2S1 -> S2
# and S3 + S1 -> S2 becomes S1 + S3 -> S2
function normal_form(n_species::Integer, reactions::Vector{Reaction})
    ∇r,∇p = stoichiometric_matrix(n_species, reactions)
    kf = [r.kf for r in reactions]
    kr = [r.kr for r in reactions]
    nm = [r.names for r in reactions]
    return from_stoichiometric_matrix(∇r,∇p,kf,kr,nm)
end 

function stoichiometric_nullspace(n_species::Integer, reactions::Vector{Reaction}, tr::Bool)
    # We only want net stoichiometric matrix here.
    ∇r, ∇p = stoichiometric_matrix(n_species, reactions)
    ∇ = ∇r - ∇p

    ∇ = tr ? transpose(∇) : ∇

    S = MatrixSpace(ZZ, size(∇)...)

    r,ns = nullspace(S(full(∇)))
    m,n = size(ns)

    nsm = sparse(Int[ns[i,j] for i=1:m,j=1:n])

    return tr ? transpose(nsm) : nsm
end 

# returns dimensions of nullspace (dimension of kernel) without calculating
# the matrix explicitly. Used for 
function stoichiometric_dimker(n_species::Integer, reactions::Vector{Reaction}, tr::Bool)
    ∇r, ∇p = stoichiometric_matrix(n_species, reactions)
    ∇ = ∇r - ∇p

    ∇ = tr ? transpose(∇) : ∇

    return size(∇,2)-rank(full(∇))
end

conservation_laws(n_species::Integer, reactions::Vector{Reaction}) = 
    stoichiometric_nullspace(n_species, reactions, true)

n_conservation_laws(n_species::Integer, reactions::Vector{Reaction}) = 
    stoichiometric_dimker(n_species, reactions, true)

cycles(n_species::Integer, reactions::Vector{Reaction}) =
    stoichiometric_nullspace(n_species, reactions, false)

n_cycles(n_species::Integer, reactions::Vector{Reaction}) = 
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
function log_equilibrium_state(n_species::Integer, reactions::Vector{Reaction})
     # compute ratio of forward to backward reactions
    free_energy = log.([r.kf/r.kr for r in reactions])

    dp,dn = stoichiometric_matrix(n_species,reactions)
    del = dn-dp

    logz = full(del')\free_energy

    return logz, del    
end 
function equilibrium_state(n_species::Integer, reactions::Vector{Reaction}) #, z0)
    logz, del = log_equilibrium_state(n_species, reactions)
    return exp.(logz)
end

# Note: if net doesn't have equilibrium state, then above functions
# will not return a steady state.
# This function however will always return a steady state
function steady_state(n_species::Integer, reactions::Vector{Reaction})
    
    z0 = equilibrium_state(n_species, reactions)

# todo: termination criteria
    for i=1:10000 
       z0 += 0.005*mass_action(reactions, z0)
    end

    return z0
end 

# The equilibrium state is often under-defined;
# this means we can add chemostats until it is fully defined,
# and get a unique solution.
# In addition, we can return the full solution space, which is:
# ss=nullspace(full(del'))
# Thus every solution can be written as:
# ss*x + logz0
# where x is some arbitrary vector and logz0=full(del')\log_k_ratio
function equilibrium_state_space(n_species::Integer, reactions::Vector{Reaction})
    logz, del = log_equilibrium_state(n_species, reactions)

    r,ss = nullspace(full(del'))
    return ss, logz
end

# -----------
function cycle_affinities(n_species::Integer, reactions::Vector{Reaction}, cycm::AbstractMatrix, z)
    cycaff = zeros(size(cycm,2))

    lnj = zeros(length(reactions))
    for rho = 1:length(reactions)
        ri = reactions[rho]
        jf,jr = reaction_currents(ri, z)
        lnj[rho] = log(jf/jr)
    end 

    return lnj'*cycm
end