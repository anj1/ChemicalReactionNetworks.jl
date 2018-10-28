# Code for working with complexes,
# and complex-balanced systems

#=
immutable SpeciesComplex 
    spec    # vector of species
    stoich  # vector of stoichiometric coefficients of species
end 
=#

#Base.isequal(a::SpeciesComplex,b::SpeciesComplex) = isequal(a.spec,b.spec) && isequal(a.stoich,b.stoich)
#Base.hash(a::SpeciesComplex) = hash([hash(a.spec),hash(a.stoich)])

# (species,stoichiometric coefficients)
typealias SpeciesComplex Tuple{Array{Int64,1},Array{Int64,1}}

# composition matrix has dimensions n_species, n_complexes
function composition_matrix(n_species, cmplx::Vector{SpeciesComplex})
    Γ = spzeros(Int64,n_species, length(cmplx))
    for i in 1:length(cmplx) 
        spec   = cmplx[i][1]
        stoich = cmplx[i][2]
        for j = 1 : length(spec)
            Γ[spec[j],i] = stoich[j]
        end
    end 
    return Γ
end

# return composition and incidence matrices of reaction net.
# composition matrix has dimensions n_species, n_complexes.
# incidence matrix has dimensionality n_complex, n_reactions,
# the entries indicate if a complex is the reactant or product
function complex_decomposition(n_species,reactions::Vector{Reaction})
    c = SpeciesComplex[]
    inc_mat = spzeros(Int64,2*length(reactions),length(reactions))  # we can have at most 2*rho complexes
    n_complex = 0

    idx = 0 

    for i = 1:length(reactions)
        r = reactions[i]

        # Reactant complex
        rc = (r.reactants, r.stoichr)

        # Find complex in list of complexes
        idx = findfirst(c,rc)
        if idx==0  # If not already in, put it in
            append!(c, [rc])
            n_complex += 1
            idx = n_complex
        end

        inc_mat[idx, i] += -1

        # Product complex
        pc = (r.products, r.stoichp)

        # Find complex in list of complexes
        idx = findfirst(c,pc)
        if idx==0  # If not already in, put it in
            append!(c, [pc])
            n_complex += 1
            idx = n_complex
        end 

        inc_mat[idx, i] +=  1
    end

    return composition_matrix(n_species, c), inc_mat[1:n_complex,:]
end 

#=

# return the set of complexes of a reaction net.
function complexes(reactions::Vector{Reaction})
    # reactant complexes
    cr = SpeciesComplex[(r.reactants,r.stoichr) for r in reactions]
    cp = SpeciesComplex[(r.products, r.stoichp) for r in reactions]
    # combine and sort out unique ones
    c = unique(cat(1,cr,cp))
    # remove empty complexes
    filter!(x -> x!=(Int64[],Int64[]), c)
    return c
end




function incidence_matrix(n_species, reactions::Vector{Reaction})
    ∇r, ∇p = stoichiometric_matrix(n_species, reactions)

    # Composition matrix
    comp_mat = unique(sparse(cat(2, ∇r, ∇p)), 2)

    # Incidence matrix

end


# returns composition matrix and incidence matrix,
    # get list of complexes
    cmplx = complexes(reactions)
=#

# Compute deficiency of reaction network
function deficiency(n_species, reactions::Vector{Reaction})
    comp_mat, inc_mat = complex_decomposition(n_species,reactions)

    return n_cycles(n_species,reactions) - (size(inc_mat,2)-rank(full(inc_mat)))
end
