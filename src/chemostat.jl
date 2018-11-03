# Functions for chemostatting reactants,
# and producing simplified reaction nets.

# Make a chemostatted version of the reaction system,
# setting species spec to concentration z
function chemostatted!(reactions::Vector{Reaction}, spec::Vector{Int}, z)
    for i = 1:length(reactions)
        r = reactions[i]

        nr = length(r.reactants)
        j=1
        while (j <= nr)
            idx = findfirst(x->x==r.reactants[j], spec)
            if idx>0
                nr -= 1
                r.kf *= z[idx]^r.stoichr[j]

                deleteat!(r.reactants, j)
                deleteat!(r.stoichr,   j)
            else
                j += 1
            end
        end 

        np = length(r.products)
        j=1
        while (j <= np)
            idx = findfirst(x->x==r.products[j], spec)
            if idx>0
                np -= 1
                r.kr *= z[idx]^r.stoichp[j]

                deleteat!(r.products, j)
                deleteat!(r.stoichp,  j)
            else
                j += 1
            end
        end 
    end
end

# re-label species according to permutation prm,
# i.e. [1,3,2] maps [1,2,3]->[1,3,2]
function relabel!(reactions::Vector{Reaction}, prm)
    for i = 1:length(reactions)
        r = reactions[i]

        nr = length(r.reactants)
        for j=1:nr
            r.reactants[j] = prm[r.reactants[j]]
        end

        np = length(r.products)
        for j=1:np
            r.products[j] = prm[r.products[j]]
        end     
    end
end