using ChemicalReactionNetworks

reactions = [Reaction([1,2],[1,1],[3],[1],1.0,1.0),
             Reaction([1,1],[1,1],[2],[2],1.0,1.0)]

cl = conservation_laws(3, reactions)

assert(full(cl) == [1 1 2])

reactions = [Reaction([1,2],[1,1],[3],[1],1.0,1.0),
             Reaction([1,1],[1,1],[2],[1],1.0,1.0)]

cl = conservation_laws(3, reactions)

assert(full(cl) == [1 2 3])

# Figure 1 in Rao paper 

reactions = [Reaction([1],[1],[2,3],[2,1],1.0,1.0),
             Reaction([3,4],[1,1],[5],[1],1.0,1.0)]

cmm=conservation_laws(5, reactions)
cmm[2,:]=cmm[2,:]+2*cmm[1,:]
cmm[3,:]=2*cmm[3,:]+cmm[2,:]-2*cmm[1,:]

assert(full(cmm) ==
[0  0   0  1  1
 0  1  -2  2  0
 2  1   0  0  0])

# Example 5 in Rao paper, without chemostatting
# species: [Y1, Xa, Xb, Xc, Y2]
reactions = [Reaction([1,2],[1,1],[3],[1],1.0,1.0),
             Reaction([3],[1],[2,5],[1,1],1.0,1.0),
             Reaction([2,5],[1,1],[4],[1],1.0,1.0),
             Reaction([4],[1],[1,2],[1,1],1.0,1.0)]

clm = cycles(5, reactions)

assert(vec(full(clm)) == [1,1,1,1])


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

clm = cycles(5, reactions)

# Note that Rao's conclusion of 'one linearly independent emergent cycle (20) arises'
# is actually wrong; there are two emergent cycles, but they can be combined into one.
assert(full(clm[1:4,:]) ==
[ 0   1
  0   1
  1   0
  1   0])

# There is also only one conservation law in the chemostatted system 
clm = conservation_laws(5, reactions)
assert(vec(full(clm)) == [0,1,1,1,0])

#-----
# Test edge cases

# Case where one species does not react
reactions = [Reaction([1],[1],[2],[1],1.0,1.0),
             Reaction([2],[1],[4],[1],1.0,1.0)]
clm = conservation_laws(4, reactions)
assert(full(clm) ==
[0  0  1  0
 1  1  0  1])

#-----
# Test complex decomposition
cm,imm=complex_decomposition(5, reactions)
dp,dr=stoichiometric_matrix(5, reactions)
assert(cm*imm == dr-dp)


# Test (29) in paper
# [Xa,Xb,Xc]
reactions = [Reaction([1],[1],[2],[1],1.0,1.0),
             Reaction([1,2],[1,1],[2],[2],1.0,1.0),
             Reaction([2],[2],[3],[1],1.0,1.0)]
cm,imm=complex_decomposition(3, reactions)
dp,dr=stoichiometric_matrix(3, reactions)
assert(cm ==
[ 1  0  1  0  0
  0  1  1  2  0
  0  0  0  0  1
])
assert(imm == 
[ -1   0   0
   1   0   0
   0  -1   0
   0   1  -1
   0   0   1
])
assert(cm*imm == dr-dp)

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
cm,imm=complex_decomposition(1, reactions)
assert(cm ==  [0  1  2])
assert(imm == 
[ -1   0   1
   1  -1   0
   0   1  -1])


# ---------------
# Test chemostatting, with previous example.
####
reactions = [Reaction([1],[1],[2],[1],1.0,1.0),
             Reaction([1,2],[1,1],[2],[2],1.0,1.0),
             Reaction([2],[2],[3],[1],1.0,1.0)]
chemostatted!(reactions, [1,3],[1.0,1.0])
relabel!(reactions, [0,1,0])
cm,imm=complex_decomposition(1, reactions)
assert(cm ==  [0  1  2])
assert(imm == 
[ -1   0   1
   1  -1   0
   0   1  -1])


#----------------
# [Ya,Xb,Xc,Xd,Ye]
n_species=5
reactions = [Reaction([1],[1],[2],[1],1.0,1.0),
             Reaction([2],[1],[3,4],[1,1],1.0,1.0),
             Reaction([3,4],[1,1],[5],[1],1.0,1.0)]


assert(deficiency(n_species,reactions) == 0)

chemostatted!(reactions, [1,5], [1.0, 1.0])

relabel!(reactions, [0,1,2,3,0])

assert(deficiency(n_species,reactions) == 0)

#------------------
# Futile cycle enzyme example.
# The deficiency is 1.
n_species=4
reactions = [Reaction([1,2],[1,1],[3,2],[1,1],1.0,1.0),
             Reaction([3,4],[1,1],[1,4],[1,1],1.0,1.0)]

assert(deficiency(n_species,reactions) == 1)

#-------------------
# Examples given on:
# https://reaction-networks.net/wiki/Reaction_graph#Deficiency

n_species=2
reactions = [Reaction([],[],[1],[1],1.0,0.0),
             Reaction([1],[1],[2],[1],1.0,0.0),
             Reaction([2],[1],[],[],1.0,0.0)]

assert(deficiency(n_species,reactions) == 0)

n_species=2
reactions = [Reaction([1],[2],[2],[2],1.0,0.0),
             Reaction([2],[2],[1,2],[1,1],1.0,0.0),
             Reaction([1,2],[1,1],[1],[2],1.0,0.0)]
assert(deficiency(n_species,reactions) == 1)

n_species=2
reactions = [Reaction([1,2],[2,1],[1],[3],1.0,0.0),
             Reaction([1],[3],[1,2],[1,2],1.0,0.0),
             Reaction([1,2],[1,2],[2],[3],1.0,0.0),
             Reaction([2],[3],[1,2],[2,1],1.0,0.0)]
assert(deficiency(n_species,reactions) == 2)