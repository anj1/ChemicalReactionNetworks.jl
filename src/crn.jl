# This file contains the base definitions,
# Plus mass-action kinetics.

import Base.==,Base.hash,Base.reverse 

"""
    ChemicalReactionNetworks.Reaction 

Note: the default definition of a reaction is reversible.
This simplifies analyses for most reactions.
Note that no truly irreversible reactions exist in nature,
And in fact true irreversible reactions _can't_ exist because
this would violate microscopic reversibility.
However, conceptually it can sometimes be convenient to
think of a reaction as being purely unidirectional,
in which case one can set kf or kr to 0.0, which is treated
as a special case.
"""
const IntVec = AbstractArray{Unsigned,1}
const StrVec = AbstractArray{String,1}
mutable struct Reaction
    reactants::IntVec # vector of reactants
    stoichr::IntVec   # vector of stoichiometric coefficients of reactants
    products::IntVec  # vector of products
    stoichp::IntVec   # vector of stoichiometric coefficients of products
    kf::Number        # rate constant, forward. Zero means no forward reaction.
    kr::Number        # rate constant, backward. Zero means no backward reaction.
end

struct ReactionNetwork 
    n_species::Integer 
    reactions::Vector{Reaction}
    names::StrVec     # Species names
end 

ReactionNetwork(ns::Integer,r::Vector{Reaction},names=String[]) = 
    ReactionNetwork(ns,r,names)

ReactionNetwork(r::Vector{Reaction},names=String[]) = 
    ReactionNetwork(get_nspecies(r),r,names)

function get_nspecies(reactions::Vector{Reaction})
    ns = 0 
    for r in reactions
        ns = maximum(cat(dims=1,ns,r.reactants,r.products))
    end 
    return ns 
end

function show_reactants(io::IO, reactants::IntVec, stoich::IntVec, names::StrVec)
    strs = String[] 

    for si in 1:length(reactants)
        s = reactants[si]

        cstr = (stoich[si] == 1) ? "" : "$(stoich[si])"

        push!(strs, isempty(names) ? "$cstr(S$s) " : "$cstr$(names[s]) ")
    end 

    print(io, join(strs, "+ "))
end 

function Base.show(io::IO, r::Reaction, names=String[])
    showdict = Dict((true, true ) => ("↮ ",""),
                    (true, false) => ("← ","[k-:$(r.kr)]"),
                    (false,true ) => ("→ ","[k+:$(r.kf)]"),
                    (false,false) => ("⇋ ","[k+:$(r.kf) k-:$(r.kr)]"))
    eq_sym, k_str = showdict[(r.kf==0, r.kr==0)]

    # Show left-hand side
    show_reactants(io, r.reactants, r.stoichr, names)

    print(io, eq_sym)

    # Show right-hand side 
    show_reactants(io, r.products, r.stoichp, names)

    print(io, k_str)
end 

function Base.show(io::IO, rn::ReactionNetwork)
    println(io, "Reaction Network: ")
    for r in rn.reactions 
        show(io, r, rn.names)
        println(io)
    end
end 

function reaction_current_1side(reactants::IntVec, stoich::IntVec, z)
    j = 1

    for si = 1 : length(reactants)
        s = reactants[si]
        c = stoich[si]

        j *= z[s] ^ c 
    end

    return j 
end 

function reaction_current(ri::Reaction, z)
    jf = ri.kf*reaction_current_1side(ri.reactants, ri.stoichr, z)
    jr = ri.kr*reaction_current_1side(ri.products,  ri.stoichp, z)
    return jf,jr 
end 

function mass_action_1side!(dz, reactants::IntVec, stoich::IntVec, j, pos::Bool)
    for si = 1 : length(reactants)
        s = reactants[si]
        c = stoich[si]

        dz[s] += pos ? c*j : -c*j
    end
end 

function mass_action!(dz, ri::Reaction, z)
    jf,jr = reaction_current(ri, z)
    j=jf-jr

    # calculate change in concentration
    mass_action_1side!(dz, ri.reactants, ri.stoichr, j, false)

    mass_action_1side!(dz, ri.products,  ri.stoichp, j, true)
end 

function precompute_concpowers(z, n_powers)
    z_powers = zeros(length(z), n_powers)

    z_powers[:,1] = z
    for c = 2 : n_powers
        z_powers[:,c] = z_powers[:,c-1] .* z
    end

    return z_powers
end 

# evolve concentration of reactants based on mass-action kinetics.
# re: list of reactions
# z: concentration of reactants
# returns: dz, instantanous change in z over time
function mass_action(reactions::Vector{Reaction}, z)
    dz = zero(z)

    for ri in reactions
        mass_action!(dz, ri, z)
    end

    return dz 
end

function jacobian(rn::ReactionNetwork)
    jm = zeros(rn.n_species,rn.n_species)

    cm = mass_action(rn.reactions,zeros(rn.n_species))

    for i=1:rn.n_species
        z = zeros(rn.n_species)
        z[i]=1.0
        jm[:,i] = mass_action(rn.reactions,z)
    end

    return cm,jm
end

reverse(r::Reaction) = Reaction(r.products,r.stoichp,r.reactants,r.stoichr,r.kr,r.kf,r.names)

# Function to determine if a reaction is redundant;
# i.e. A + 2B -> 2B + A is redundant.
# Note that for this to work, the reaction must be in normal form,
# Which means the reactants and products must be a sorted list.
isredundant(r::Reaction) = (r.products == r.reactants) && (r.stoichp == r.stoichr )

function test_equality(r::Reaction,p::Reaction)
    return (r.stoichr == p.stoichr) && 
           (r.stoichp == p.stoichp) && 
           (r.reactants == p.reactants) && 
           (r.products  == p.products) &&
           (r.kf == p.kf) &&
           (r.kr == p.kr)
end

# Two reactions are equal if they are equal in either
# the forward or backward directions
function ==(r::Reaction,p::Reaction)
    return test_equality(r,p) || test_equality(reverse(r),p)
end

function hash1(r::Reaction)
    return foldr(hash, [r.reactants,r.stoichr,r.products,r.stoichp,r.kf,r.kr,r.names]; init=zero(UInt64))
end
function hash(r::Reaction)
    h1 = hash1(r)
    h2 = hash1(reverse(r))
    return h1+h2
end