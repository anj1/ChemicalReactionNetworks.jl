ChemicalReactionNetworks.jl
------

[![Build Status](https://travis-ci.org/anj1/ChemicalReactionNetworks.jl.svg?branch=master)](https://travis-ci.org/anj1/ChemicalReactionNetworks.jl)
[![Coverage Status](https://coveralls.io/repos/github/anj1/ChemicalReactionNetworks.jl/badge.svg?branch=master)](https://coveralls.io/github/anj1/ChemicalReactionNetworks.jl?branch=master)

#### Introduction

This is a package for simulating Chemical Reaction Networks, of the kind that would be encountered in the study of complex metabolic networks, self-replication, and so on. A high-altitude overview of the capabilities of this package are:

- Basic functions: Specifying reaction networks and solving their time-evolution via mass-action kinetics.
- Reaction network manipulation: Adding external sources or sinks for chemicals, grouping species in terms of chemical complexes.
- Numerical functions: Finding equilibrium/steady states, calculating the contribution of various metabolic cycles.
- Investigation: Identifying conserved species, identifying closed and open cycles.

#### Detailed Reference

A chemical reaction network (CRN from now on) is a list of *reactions*, each reaction having a set of *reactants* and a set of *products*, with each reactant or product being a distinct chemical species. Reactions can be either *reversible* or *irreversible*, with reversible reactions capable of happening in both directions, with a characteristic rate for either direction. An irreversible reaction can be thought of as a reversible reaction with a very small rate in the opposite direction, and in fact in this package the default way of representing reactions is that they can happen in both directions. From this viewpoint, the distinction between 'reactant' and 'product' is arbitrary and depends which direction one is considering. A detailed overview of CRN theory is beyond the scope of this document but excellent references can be found in [[1](#gun2003)] and [[2](#feinberg1979)].

In CRNs, reactions are rarely considered in isolation; one is almost always interested in how different reactions interact and compete with one another in a closed or open system. Because of this, it is assumed that each reaction net has a 'global' list of reactants, and so in each reaction one does not need to explicitly specify the reactants by composition, instead one just refers to them. Thus, for instance, if our list of reactants is: 1. Methane (CH₄), 2. Dioxygen (O₂), 3. Carbon dioxide (CO₂), and 4. H₂O, one can specify the combustion of methane;

<p style="text-align: center;">
 CH₄ + 2O₂ → CO₂ + 2H₂O
</p>

As the following set of lists:

- Reactants: `[1, 2]`
- Stoichiometric coefficients of reactants: `[1, 2]`
- Products: `[3, 4]`
- Stoichiometric coefficents of products: `[1, 2]`

And in this package, this is represented as follows:

```julia
Reaction([1, 2], [1, 2], [3, 4], [1, 2], k, 0)
```

With `k` being the forward reaction rate, and 0 being the effective reverse reaction rate. Note that `k` is dependent on temperature, and at high enough temperature there is a non-zero reverse reaction rate (at very high temperature the forward and reverse rates are effectively the same).

For various reasons, it is advantageous to index species by number rather than by a pointer, symbolic reference, or other kind of arbitrary reference. For example, it makes explicit the ordering used when doing various algebraic manipulations on the net, as shall be seen below.

#### Time-Evolution of CRNs

Perhaps the most basic thing one might want to do with a CRN is to see how the net evolves in time -- how species are produced and consumed, and what the final concentration of species is. The function to do this is `mass_action`, and it works as follows:

```julia
dz = mass_action(reactions, z)
```

Where reactions is a list of `Reaction`s, and `z` is some instantaneous concentration of all the chemical species. This function returns the *rate of change* of the instantaneous concentration. Thus by feeding this function into an ODE solver, such as a stiff 2nd-order solver, one can solve for the full time-dynamics of the network. An example is given below.


#### Steady States

One can find the steady states 

A special class of steady states are *equilibrium* states; these are states where each reaction happens equally in the forward and backward directions. Obviously equilibrium states can only occur when all reactions are reversible. If a CRN has an equilibrium state, it can be found by solving a linear problem; the function `equilibrium_state` does this:

```julia
equilibrium_state(n_species, reactions)
```

If the CRN does not have an equilibrium state, then the above function will not return a steady state. In such a scenario, `SteadyStateProblem` from DifferentialEquations.jl can be used to find a steady state by iteratively solving until converging, starting from an initial guess. To do this we wrap `mass_action` in an equations-of-motion function, which also ensures that only positive solutions are found:

```julia
function eom(z,p,t)
    dz = mass_action(reactions,z)
    for i=1:length(z)
        if z[i] < 0
            dz[i] = exp(-z[i])
        end
    end
    return dz
end
steady_z = solve(SteadyStateProblem(f, z0))
```



#### Catalytic System

Consider the simple catalytic system:

```julia
chems = ["A₂","B","AB","Z","AZ","BZ"]
reactions = [Reaction([1,4],[1,2],[5],[2],1.0,1.0,chems),
             Reaction([2,4],[1,1],[6],[1],1.0,1.0,chems),
             Reaction([5,6],[1,1],[3,4],[1,2],1.0,0.0,chems)]
@show reactions
```
```
3-element Array{Reaction,1}:
 A₂ + 2Z ⇋ 2AZ [k+:1.0 k-:1.0]
 B + Z ⇋ BZ [k+:1.0 k-:1.0]   
 AZ + BZ → AB + 2Z [k+:1.0]   
```

This is a simple chemical system that has multiple steady states. 


#### Pseudo-Reactions and Chemostatting

Before diving into the more advanced abilities of this package, we must first describe an important concept that is used often: *pseudo-reactions*. A pseudo-reaction is a reaction that involves *hidden* or unseen reactants. For example, we can take the combustion of methane:

<p style="text-align: center;">
 CH₄ + 2O₂ → CO₂ + 2H₂O
</p>

 and make O₂ an unseen reactant:

<p style="text-align: center;">
 CH₄ → CO₂ + 2H₂O
</p>

Note that this is no longer a **physical** reaction - methane does not spontaneously produce CO₂ in the absence of oxidant. If you were to write down this reaction in a high school chemistry exam, you would probably fail. However, there is nothing stopping us from implementing this as a reaction in code, and indeed, from the code's point of view, this is just as good as any other reaction!

Why would one do this? One situation where pseudo-reactions are useful is when describing situations where some of the reactants are so *abundant* and available (either by being present as well-mixed solvent, or being deliberately replenished or removed externally to maintain a constant concentration) that we simply model their concentration by a fixed value. In these kinds of cases, the reaction rate is thus solely determined by the other, non-hidden reactants. Note that a pseudo-reaction does not need to have any reactants at all! Thus the following:

<p style="text-align: center;">
⇋ H₂O
</p>

Given by `Reaction([],[],[1],[1],1,1)` (With forward and backward rate of 1, and assuming water is reactant no. 1) is a perfectly acceptable reaction. It can be interpreted as 'water is always present in a constant amount' and is useful in CRNs that take place in an aqueous environment and in which some species may react with water. We can imagine that the empty side of the reaction is referring to something happening 'off-screen' that is constantly replenishing our supply of water. Another example might be:

<p style="text-align: center;">
CO₂ →
</p>

Which may, for example, represent CO₂ bubbling out of solution and thus its concentration always being small in solution.

More generally, the idea of fixing the concentration of a species externally is called *chemostatting*, and pseudo-reactions are a way of performing chemostatting without having to write special cases in code. A CRN where some reactants are chemostatted is called an *open* CRN. Otherwise, it is called a *closed* CRN.

In this package, one can fix the concentration of a species by calling the function `chemostatted!`. This function takes a CRN as input, along with a set of species to chemostat and what concentrations to chemostat them to. It then modifies the CRN by removing all the chemostatted species, and adjusting the rate constants of the reactions they participate in so that the network dynamics are the same as before.

For example:


#### Quasi-Steady Simulation

Some CRNs can be very complicated and include many reactions happening at different timescales. A way of simplifying the time-evolution of these CRNs is to take the reactions that happen quickly as happening *instantaneously*, so that both sides of those reactions are set to be equal, and thus those reactions (along with any species that only show up in those reactions) are eliminated. This is the *quasi-steady* simplification. This can be very useful, for instance, if the intermediate reactions involve short-lived, high-energy radicals.

In ChemicalReactionNetworks, a CRN can be simplified to a quasi-steady CRN by using the function `quasi_steady`. This function takes the CRN and some reaction rate k, and eliminates all reactions that happen at rate faster than k.

Note that in general this process isn't as simple as just equating the products and reactants, as the law of mass action i

#### Inquiring Reaction Nets

Given a CRN, one can find its *conservation laws*, which describe the underlying species that are conserved by all the reactions in the CRN. For instance, going back to the example of combustion of methane, 


Similarly, one can find the *cycles* of a CRN, which are the set of reactions that, when performed in a certain order and a certain number of times, result in a *zero* change of concentration of the reactants. For example:



There is far more capability in this package than merely calculating time-dynamics.

#### Demos

This package has been used to implement the models presented in the following research papers:

- "Design of conditions for emergence of self-replicators"
- "Spontaneous fine-tuning to environment in many-species chemical reaction networks"

#### References

1. <a name="gun2003">Jeremy Gunawardena (2003).</a> ["Chemical reaction network theory for in-silico biologists"](http://vcp.med.harvard.edu/papers/crnt.pdf).
2. <a name="feinberg1979">Martin Feinberg (1979).</a> ["Lecture notes for chemical reaction networks"](https://crnt.osu.edu/LecturesOnReactionNetworks)
