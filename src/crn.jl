# This file contains the base definitions,
# Plus mass-action kinetics.

import Base.==,Base.hash,Base.reverse 

# Note: the default definition of a reaction is reversible.
# This simplifies analyses for most reactions.
# Note that no truly irreversible reactions exist in nature,
# And in fact true irreversible reactions _can't_ exist because
# this would violate microscopic reversibility.
# However, conceptually it can sometimes be convenient to
# think of a reaction as being purely unidirectional,
# in which case one can set kf or kr to 0.0, which is treated
# as a special case.
const IntVec = AbstractArray{Unsigned,1}
const StrVec = AbstractArray{String,1}
mutable struct Reaction
    reactants::IntVec # vector of reactants
    stoichr::IntVec   # vector of stoichiometric coefficients of reactants
    products::IntVec  # vector of products
    stoichp::IntVec   # vector of stoichiometric coefficients of products
    kf::Number        # rate constant, forward. Zero means no forward reaction.
    kr::Number        # rate constant, backward. Zero means no backward reaction.
    names::StrVec     # Array of names of reactants. Can be empty.
end

Reaction(a,b,c,d,kf,kr) = Reaction(a,b,c,d,kf,kr,String[])

function Base.show(io::IO, r::Reaction)
    if (r.kf == 0.0) && (r.kr == 0.0)
        print(io, "[Non-reaction]")
        return
    end 

    for si in 1:length(r.reactants)
        s = r.reactants[si]
        c = r.stoichr[si]

        cstr = c == 1 ? "" : "$c"

        if isempty(r.names)
            print(io, "$cstr(S$s) ")
        else
            nm = r.names[s]
            print(io, "$cstr$nm ")
        end 
        if si == length(r.reactants)
            break
        end 
        print(io, "+ ")
    end 

    if ((r.kf != 0.0) && (r.kr != 0.0))
        print(io,"⇋ ")
    end 
    if (r.kr == 0.0) && (r.kf != 0.0)
        print(io,"→ ")
    end
    if (r.kf == 0.0) && (r.kr != 0.0)
        print(io,"← ")
    end

    for si in 1:length(r.products)
        s = r.products[si]
        c = r.stoichp[si]

        cstr = c == 1 ? "" : "$c"

        if isempty(r.names)
            print(io, "$cstr(S$s) ")
        else
            nm = r.names[s]
            print(io, "$cstr$nm ")
        end 
        if si == length(r.products)
            break
        end 
        print(io, "+ ")
    end 

    kf = r.kf
    kr = r.kr
    print(io, "[")
    if kf != 0.0
        print(io, "k+:$kf")
    end
    if kr != 0.0
        if kf != 0.0
            print(io, " ")
        end
        print(io, "k-:$kr")
    end 
    print(io, "]")
end 

function reaction_current(ri::Reaction, z)
    jf = ri.kf
    for si = 1 : length(ri.reactants)
        s = ri.reactants[si]
        c = ri.stoichr[si]

        jf *= z[s] ^ c 
    end

    jr = ri.kr
    for si = 1 : length(ri.products)
        s = ri.products[si]
        c = ri.stoichp[si]

        jr *= z[s] ^ c
    end

    return jf,jr 
end 

function mass_action!(dz, ri::Reaction, z)
    jf,jr = reaction_current(ri, z)
    j=jf-jr

    # calculate change in concentration
    for si = 1 : length(ri.reactants)
        s = ri.reactants[si]
        c = ri.stoichr[si]

        dz[s] -= c*j
    end

    for si = 1 : length(ri.products)
        s = ri.products[si]
        c = ri.stoichp[si]

        dz[s] += c*j
    end
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

function jacobian(n_species::Integer,reactions::Vector{Reaction})
    jm = zeros(n_species,n_species)

    cm = mass_action(reactions,zeros(n_species))

    for i=1:n_species
        z = zeros(n_species)
        z[i]=1.0
        jm[:,i] = mass_action(reactions,z)
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