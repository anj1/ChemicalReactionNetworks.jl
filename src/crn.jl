# This file contains the base definitions,
# Plus mass-action kinetics.

# Note: the default definition of a reaction is bidirectional.
# This simplifies analyses for most reactions.
# Note that no truly unidirectional reactions exist in nature,
# And in fact unidirectional reactions _can't_ exist because
# this would violate microscopic reversibility.
# However, unidirectional reactions can be simulated by setting
# kf or kr to 0.0, which most functions can deal with correctly.
type Reaction
    reactants    # vector of reactants
    stoichr      # vector of stoichiometric coefficients of reactants
    products     # vector of products
    stoichp      # vector of stoichiometric coefficients of products
    kf::Real     # rate constant, forward. Zero means no forward reaction.
    kr::Real     # rate constant, backward. Zero means no backward reaction.
end

type DrivenReaction
    r::Reaction  # base reaction
    j_coupling   # coupling matrix
    c_offset     # offset matrix
end

function Base.show(io::IO, r::Reaction)
    if (r.kf == 0.0) && (r.kr == 0.0)
        print(io, "[Non-reaction]")
        return
    end 

    for si in 1:length(r.reactants)
        s = r.reactants[si]
        c = r.stoichr[si]

        print(io, "$c(S$s) ")
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

        print(io, "$c(S$s) ")
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
    println(io, "]")
end 

function Base.show(io::IO, re::DrivenReaction)
    print(io, "driven reaction: ")
    show(io, re.r)
end 

function concentration_currents(ri::Reaction, z)
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
    jf,jr = concentration_currents(ri, z)
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

function mass_action(reactions::Vector{DrivenReaction}, z)
    dz = zero(z)

    for ri in reactions 
        mass_action!(dz, ri.r, z)
    end

    return dz
end
