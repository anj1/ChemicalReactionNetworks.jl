# (Linear) algebraic functions for Chemical Reaction Nets.
# Compute conservation laws, cycles, emergent cycles.

using Nemo 
using SparseArrays 

full(a) = convert(Array,a)

"""
    stoichiometric_matrix(rn::ReactionNetwork)

Generate stoichiometric matrices from list of reactions.
Output matrices are of dimension: (`n_species`,`n_reactions`).

Two matrices are returned, for reactants and products.
"""
function stoichiometric_matrix(rn::ReactionNetwork)
    n_species = rn.n_species 
    reactions = rn.reactions 

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

function findnz(l::AbstractVector)
    idx = findall(x->x!=0,l)
    return idx,l[idx]
end 

"""
    from_stoichiometric_matrix(∇r, ∇p, kf, kr, names=[])

Generates a list of reactions from stoichiometric matrices.
Inputs:
- ∇r, ∇p: stoichiometric matrices for reactants and products, respectively.
- kf, kr: Forward and backward rate constants.
- names (optional): names of the species, to be included in output network.
Returns:
- ReactionNetwork
"""
function from_stoichiometric_matrix(∇r::AbstractMatrix, ∇p::AbstractMatrix, kf::Vector{T}, kr::Vector{T}, names=[]) where T<:Real
    @assert size(∇r)==size(∇p)
    nspecies,nreactions = size(∇r)

    reactions = Reaction[]
    for i = 1 : nreactions
        reactants,stoichr = findnz(∇r[:,i])
        products,stoichp  = findnz(∇p[:,i])
        push!(reactions, Reaction(reactants,stoichr,products,stoichp,kf[i],kr[i]))
    end

    return ReactionNetwork(nspecies,reactions,names)
end

"""
    normal_form(rn::ReactionNetwork)

Returns normalized form of reaction net.

For example, S1 + S1 → S2 becomes 2S1 → S2
and S3 + S1 → S2 becomes S1 + S3 → S2
"""
function normal_form(rn::ReactionNetwork)
    ∇r,∇p = stoichiometric_matrix(rn)
    kf = [r.kf for r in rn.reactions]
    kr = [r.kr for r in rn.reactions]
    return from_stoichiometric_matrix(∇r,∇p,kf,kr,rn.names)
end 

# TODO: Bool arguments are a code smell
function stoichiometric_nullspace(rn::ReactionNetwork, tr::Bool)
    # We only want net stoichiometric matrix here.
    ∇r, ∇p = stoichiometric_matrix(rn)
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
function stoichiometric_dimker(rn::ReactionNetwork, tr::Bool)
    ∇r, ∇p = stoichiometric_matrix(rn)
    ∇ = ∇r - ∇p

    ∇ = tr ? transpose(∇) : ∇

    return size(∇,2)-rank(full(∇))
end

"""
    conservation_laws(rn::ReactionNetwork)

Returns the list of 'conservation laws' of the network, 
where each law is defined as a vector of length `n_species`.
"""
conservation_laws(rn::ReactionNetwork) = 
    stoichiometric_nullspace(rn, true)

n_conservation_laws(rn::ReactionNetwork) = 
    stoichiometric_dimker(rn, true)

"""
    cycles(rn::ReactionNetwork)

Returns the list of 'cycles' of the network, 
where each law is defined as a vector of length `n_reactions`.
"""
cycles(rn::ReactionNetwork) =
    stoichiometric_nullspace(rn, false)

n_cycles(rn::ReactionNetwork) = 
    stoichiometric_dimker(rn, false)

# The solution is based on equation (18) in the paper:
#   kf/kr = prod_s z[s]^del_s,r for all r
# We take the logarithm, giving:
#   log(kf/kr) = sum del_s,r log z[s]
# And then we also append the conservation laws.
# This gives a linear system which can be solved.
function log_equilibrium_state(rn::ReactionNetwork)
     # compute ratio of forward to backward reactions
    free_energy = log.([r.kf/r.kr for r in rn.reactions])

    dp,dn = stoichiometric_matrix(rn)
    del = dn-dp

    logz = full(del')\free_energy

    return logz, del    
end
"""
    equilibrium_state(rn::ReactionNetwork)

Compute steady-states for both detailed-balanced and
complex-balanced states.

The steady state concentration is one where the conserved
moieties are distributed equally among all reactants.

The equilibrium state is a steady-state where all
concentration currents vanish.
"""
function equilibrium_state(rn::ReactionNetwork)
    logz, del = log_equilibrium_state(rn)
    return exp.(logz)
end

# The equilibrium state is often under-defined;
# this means we can add chemostats until it is fully defined,
# and get a unique solution.
# In addition, we can return the full solution space, which is:
# ss=nullspace(full(del'))
# Thus every solution can be written as:
# ss*x + logz0
# where x is some arbitrary vector and logz0=full(del')\log_k_ratio
function equilibrium_state_space(rn::ReactionNetwork)
    logz, del = log_equilibrium_state(rn)

    r,ss = nullspace(full(del'))
    return ss, logz
end

# -----------
function cycle_affinities(rn::ReactionNetwork, cycm::AbstractMatrix, z)
    cycaff = zeros(size(cycm,2))

    lnj = zeros(length(rn.reactions))
    for rho = 1:length(rn.reactions)
        ri = rn.reactions[rho]
        jf,jr = reaction_currents(ri, z)
        lnj[rho] = log(jf/jr)
    end 

    return lnj'*cycm
end

function specificity(rn::ReactionNetwork, z)
    # remember: r is reactants, p is products
    ∇r, ∇p = stoichiometric_matrix(rn)
    ∇rp = cat(dims=2, ∇r, ∇p)

    cur = [reaction_current(r, z) for r in rn.reactions]

    # Concatenate forward/backward currents,
    # with order corresponding to ∇rp
    curfr = cat(dims=1,[c[1] for c in cur],[c[2] for c in cur])

    sumcur = ∇rp'*(∇rp*curfr)
    spec = curfr ./ sumcur 

    return reshape(spec, length(rn.reactions), 2)
end    
