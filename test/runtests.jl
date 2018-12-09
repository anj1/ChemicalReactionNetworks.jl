using ChemicalReactionNetworks

###
chems = ["CH₄", "O₂", "CO₂", "H₂O"]
rn = ReactionNetwork([Reaction([1,2],[1,2],[3,4],[1,2],1.0,0.0)], chems)
print(rn)

###

reactions = [Reaction([1,2],[1,1],[3],[1],1.0,1.0),
            Reaction([1,1],[1,1],[2],[2],1.0,1.0)]

cl = conservation_laws(ReactionNetwork(reactions))

@assert convert(Array,cl) == [1 1 2]

reactions = [Reaction([1,2],[1,1],[3],[1],1.0,1.0),
            Reaction([1,1],[1,1],[2],[1],1.0,1.0)]

cl = conservation_laws(ReactionNetwork(reactions))

@assert convert(Array,cl) == [1 2 3]

# Figure 1 in Rao paper 

reactions = [Reaction([1],[1],[2,3],[2,1],1.0,1.0),
            Reaction([3,4],[1,1],[5],[1],1.0,1.0)]

cmm=conservation_laws(ReactionNetwork(reactions))

@assert convert(Array,cmm) ==
[ 2  1  0  0  0
 0  0  0  1  1
 1  0  1  0  1]

# Example 5 in Rao paper, without chemostatting
# species: [Y1, Xa, Xb, Xc, Y2]
reactions = [Reaction([1,2],[1,1],[3],[1],1.0,1.0),
            Reaction([3],[1],[2,5],[1,1],1.0,1.0),
            Reaction([2,5],[1,1],[4],[1],1.0,1.0),
            Reaction([4],[1],[1,2],[1,1],1.0,1.0)]

clm = cycles(ReactionNetwork(reactions))

@assert vec(convert(Array,clm)) == [1,1,1,1]


# Example 5 in Rao paper, with chemostatting.
# A chemostatted species can be represented as a 'reaction',
# with no reactants, and the species as product.
# The equilibrium concentration is k+/k-,
# And the speed at which it reaches equilibrium is proportional to k+
reactions = [Reaction([1,2],[1,1],[3],[1],1.0,1.0),
            Reaction([3],[1],[2,5],[1,1],1.0,1.0),
            Reaction([2,5],[1,1],[4],[1],1.0,1.0),
            Reaction([4],[1],[1,2],[1,1],1.0,1.0),
            Reaction([],[],[1],[1],100.0,100.0),
            Reaction([],[],[5],[1],100.0,100.0)]

clm = cycles(ReactionNetwork(reactions))

# Note that Rao's conclusion of 'one linearly independent emergent cycle (20) arises'
# is actually wrong; there are two emergent cycles, but they can be combined into one.
@assert convert(Array,clm[1:4,:]) ==
[ 0   1
  0   1
  1   0
  1   0]

# There is also only one conservation law in the chemostatted system 
clm = conservation_laws(ReactionNetwork(reactions))
@assert vec(convert(Array,clm)) == [0,1,1,1,0]

#-----
# Test edge cases

# Case where one species does not react
reactions = [Reaction([1],[1],[2],[1],1.0,1.0),
             Reaction([2],[1],[4],[1],1.0,1.0)]
clm = conservation_laws(ReactionNetwork(reactions))
@assert convert(Array,clm) ==
[ 1  1  0  1
  0  0  1  0]

#-----
# Test complex decomposition
cm,imm=complex_decomposition(ReactionNetwork(reactions))
dp,dr=stoichiometric_matrix(ReactionNetwork(reactions))
@assert cm*imm == dr-dp


# Test (29) in paper
# [Xa,Xb,Xc]
reactions = [Reaction([1],[1],[2],[1],1.0,1.0),
            Reaction([1,2],[1,1],[2],[2],1.0,1.0),
            Reaction([2],[2],[3],[1],1.0,1.0)]
cm,imm=complex_decomposition(ReactionNetwork(reactions))
dp,dr=stoichiometric_matrix(ReactionNetwork(reactions))
@assert cm ==
[ 1  0  1  0  0
  0  1  1  2  0
  0  0  0  0  1
]
@assert imm ==
[ -1   0   0
   1   0   0
   0  -1   0
   0   1  -1
   0   0   1
]
@assert cm*imm == dr-dp

# Example 9 in paper.
# Here we 'chemostat' species by simply removing them from the equations;
# incorporating the chemostatted concentrations into the rate equations.
# This procedure is described:
# "This regrouping corresponds to the equivalent CRN only made of
#  internal species with effective rate constant ruling each reaction"

# Here the chemostatted species are 1 and 3.
# The only species left is 2, which we re-label as 1.
reactions = [Reaction([ ],[ ],[1],[1],1.0,1.0),
            Reaction([1],[1],[1],[2],1.0,1.0),
            Reaction([1],[2],[ ],[ ], 1.0,1.0)]
# we can also get the above by doing:
# chemostatted!(reactions, [1,3], [1.0, 1.0])
cm,imm=complex_decomposition(ReactionNetwork(reactions))
@assert cm ==  [0  1  2]
@assert imm ==
[ -1   0   1
   1  -1   0
   0   1  -1]


# ---------------
# Test chemostatting, with previous example.
####
reactions = [Reaction([1],[1],[2],[1],1.0,1.0),
            Reaction([1,2],[1,1],[2],[2],1.0,1.0),
            Reaction([2],[2],[3],[1],1.0,1.0)]
chemostatted!(reactions, [1,3],[1.0,1.0])
relabel!(reactions, [0,1,0])
cm,imm=complex_decomposition(ReactionNetwork(reactions))
@assert cm ==  [0  1  2]
@assert imm ==
[ -1   0   1
   1  -1   0
   0   1  -1]


#----------------
# [Ya,Xb,Xc,Xd,Ye]
reactions = [Reaction([1],[1],[2],[1],1.0,1.0),
            Reaction([2],[1],[3,4],[1,1],1.0,1.0),
            Reaction([3,4],[1,1],[5],[1],1.0,1.0)]

@assert deficiency(ReactionNetwork(reactions)) == 0

chemostatted!(reactions, [1,5], [1.0, 1.0])

relabel!(reactions, [0,1,2,3,0])

@assert deficiency(ReactionNetwork(reactions)) == 0

#------------------
# Futile cycle enzyme example.
# The deficiency is 1.
n_species=4
reactions = [Reaction([1,2],[1,1],[3,2],[1,1],1.0,1.0),
            Reaction([3,4],[1,1],[1,4],[1,1],1.0,1.0)]

@assert deficiency(ReactionNetwork(reactions)) == 1

#-------------------
# Examples given on:
# https://reaction-networks.net/wiki/Reaction_graph#Deficiency

reactions = [Reaction([],[],[1],[1],1.0,0.0),
            Reaction([1],[1],[2],[1],1.0,0.0),
            Reaction([2],[1],[],[],1.0,0.0)]

@assert deficiency(ReactionNetwork(reactions)) == 0

reactions = [Reaction([1],[2],[2],[2],1.0,0.0),
            Reaction([2],[2],[1,2],[1,1],1.0,0.0),
            Reaction([1,2],[1,1],[1],[2],1.0,0.0)]
@assert deficiency(ReactionNetwork(reactions)) == 1

reactions = [Reaction([1,2],[2,1],[1],[3],1.0,0.0),
            Reaction([1],[3],[1,2],[1,2],1.0,0.0),
            Reaction([1,2],[1,2],[2],[3],1.0,0.0),
            Reaction([2],[3],[1,2],[2,1],1.0,0.0)]
@assert deficiency(ReactionNetwork(reactions)) == 2

# ----------

# example taken from Sarkar paper
reactions = [Reaction([1,1],[1,1],[2],[1],1.0,1.0),
            Reaction([1,2],[1,1],[3],[1],1.0,1.0),
            Reaction([1,3],[1,1],[4],[1],1.0,1.0),
            Reaction([1,3],[1,1],[2],[2],2.0,1.0),
            Reaction([1,4],[1,1],[2,3],[1,1],1.0,1.0),
            Reaction([2,4],[1,1],[3],[2],1.0,1.0)]
chem_potentials=ones(4)
fe=free_energies(chem_potentials, reactions)
kf,kr = reaction_rates(fe, ones(6), 1.0)
for i=1:6
    reactions[i].kf = kf[i]
    reactions[i].kr = kr[i]
end 
z=equilibrium_state(ReactionNetwork(reactions))
dz=mass_action(reactions, z)
@assert all(dz .< 1e-9)   # the network should admit an equilibrium state


# ------------------------------------
# Testing mass action: the Oregonator

chems = ["HBrO₂", "Br⁻", "Ce(IV)", "BrO₃⁻", "CH₂(COOH)₂", "HOBr"]
#chems = ["X", "Y", "Z", "A", "B", "P"]
k1 = 1.0
k2 = 1.0
k3 = 1.0
k4 = 1.0
k5 = 1.0

reactions = [Reaction([4,2],[1,1],[1,6],[1,1],k1,0),
            Reaction([1,2],[1,1],[6],[2],k2,0),
            Reaction([4,1],[1,1],[1,3],[2,2],k3,0),
            Reaction([1],[2],[4,6],[1,1],k4,0),
            Reaction([5,3],[2,2],[2],[1],k5,0)]

@show reactions


# -----------------------------------
# Autocatalysis reaction
chems = ["A₂","B","AB","Z","AZ","BZ"]
n_species=length(chems)
reactions = [Reaction([1,4],[1,2],[5],[2],1.0,1.0),
            Reaction([2,4],[1,1],[6],[1],1.0,1.0),
            Reaction([5,6],[1,1],[3,4],[1,2],1.0,1e-4)]
# Chemostat A₂, B, and AB 
chemostatted!(reactions, [1,2,3], [1.0,1.0,1.0])
f(z,p,t) = mass_action(reactions,z)
#=function eom!(znew,z,p,t)
    znew[:] = mass_action(reactions,z)
end
ds=ContinuousDynamicalSystem(eom!, ones(n_species), nothing)
=#
function eom(z,p,t)
    dz = mass_action(reactions,z)
    for i=1:length(z)
        if z[i] < 0
            dz[i] = exp(-z[i])
        end
    end
    return dz
end
# ssz=solve(SteadyStateProblem(f, equilibrium_state(6, reactions)))


# ---------------------------------
# Test specificity

reactions = [Reaction([1],[1],[2],[1],1.0,1.0),
             Reaction([1],[1],[4],[1],1.0,1.0)]
spmat = [0.5 1.0
         0.5 1.0]
isapprox(specificity(ReactionNetwork(reactions), rand(n_species)), spmat)

reactions = [Reaction([1],[1],[2],[1],1.0,1.0),
             Reaction([1],[1],[4],[1],0.0,1.0)]
spmat = [1.0 1.0
         0.0 1.0]
isapprox(specificity(ReactionNetwork(reactions), rand(n_species)), spmat)

# ----------------------------------
# Test normal form
reactions = [Reaction([1,1],[1,2],[3],[1],1.0,1.0),
             Reaction([1],[2],[4],[1],1.0,1.0)]

rnf = normal_form(ReactionNetwork(reactions))

@assert Reaction([1],[3],[3],[1],1.0,1.0) == rnf.reactions[1]
