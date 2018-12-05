# Code for working with complexes,
# and complex-balanced systems

# (species,stoichiometric coefficients)
const SpeciesComplex = Tuple{IntVec,IntVec}

# composition matrix has dimensions n_species, n_complexes
function composition_matrix(n_species::Integer, cmplx::Vector{SpeciesComplex})
    Γ = spzeros(Int64,n_species, length(cmplx))
    for i in 1:length(cmplx) 
        spec   = cmplx[i][1]
        stoich = cmplx[i][2]
        for j = 1 : length(spec)
            Γ[spec[j],i] += stoich[j]
        end
    end 
    return Γ
end

# return composition and incidence matrices of reaction net.
# composition matrix has dimensions n_species, n_complexes.
# incidence matrix has dimensionality n_complex, n_reactions,
# the entries indicate if a complex is the reactant or product
function complex_decomposition(rn::ReactionNetwork)
    n_species = rn.n_species 
    reactions = rn.reactions 

    c = SpeciesComplex[]
    inc_mat = spzeros(Int64,2*length(reactions),length(reactions))  # we can have at most 2*rho complexes
    n_complex = 0

    idx = 0 

    for i = 1:length(reactions)
        r = reactions[i]

        # Reactant complex
        rc = (r.reactants, r.stoichr)

        # Find complex in list of complexes
        idx = findfirst(x->x==rc,c)
        if idx==nothing  # If not already in, put it in
            append!(c, [rc])
            n_complex += 1
            idx = n_complex
        end

        inc_mat[idx, i] += -1

        # Product complex
        pc = (r.products, r.stoichp)

        # Find complex in list of complexes
        idx = findfirst(x->x==pc,c)
        if idx==nothing  # If not already in, put it in
            append!(c, [pc])
            n_complex += 1
            idx = n_complex
        end 

        inc_mat[idx, i] +=  1
    end

    return composition_matrix(n_species, c), inc_mat[1:n_complex,:]
end 

# Compute deficiency of reaction network
function deficiency(rn::ReactionNetwork)
    comp_mat, inc_mat = complex_decomposition(rn)

    return n_cycles(rn) - (size(inc_mat,2)-rank(full(inc_mat)))
end

# Take the species concentration vector,
# And return a vector of concentrations of complexes
# z: species concentration
complex_concentration(z, composition_matrix) = vec(exp((log(z)'*composition_matrix)))

# Create a reaction system for the complexes
function complex_basis_reactions(reactions::Unsigned, incidence_matrix)
    r = Reaction[]
    n_complexes = size(incidence_matrix,1)
    names = ["(C$x)" for x=1:n_complexes]
    for i = 1 : length(reactions)
        nzind,nzval = findnz(incidence_matrix[:,i])

        reactants = nzval[1] < 0 ? [nzind[1]] : [nzind[2]]
        products  = nzval[2] > 0 ? [nzind[2]] : [nzind[1]]

        append!(r, [Reaction(reactants,[1],products,[1], reactions[i].kf, reactions[i].kr, names)])
    end 
    return r
end 