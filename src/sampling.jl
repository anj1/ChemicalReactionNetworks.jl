# From Horowitz paper:
# choosing initial conditions uniformly
# from the simplex of species concentrations with total
# sum of concentrations fixed to 25
# (on average 1 mol/vol of each species).
function rand_concentrations(n_species, total_conc)
    # see 'Uniform sampling from a simplex' from StackExchange
    # Generate a list of sorted numbers from 0 to 1
    lst = cat(dims=1, [0.0], sort(rand(n_species-1)), [1.0])
    # Take differences between consecutive elements
    return total_conc*(lst[2:end]-lst[1:end-1])
end 