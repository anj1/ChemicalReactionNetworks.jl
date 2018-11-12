# Compute petri net adjacency matrix

using SparseArrays

# The adjacency matrix that is returned has dimensions (n_species + n_reactions)^2
function petri_net(n_species::Integer,reactions::Vector{Reaction})
	n_reactions = length(reactions)

	n = n_species + n_reactions

	# start with empty adjacency matrix
	adjm = spzeros(Int,n,n)

	for j = 1:length(reactions)
		r = reactions[j]
		nr = length(r.reactants)
		for i=1:nr
			adjm[r.reactants[i],n_species+j] = r.stoichr[i]
		end
		nr = length(r.products)
		for i=1:nr
			adjm[n_species+j,r.products[i]] = r.stoichp[i]
		end
	end

	return adjm
end