Modification to CPMD-4.1 for high-through-put alchemy
=====================================================
For first order alchemical prediction, it is necessary to load the
previously converged wavefunction file to do a single integral
evaluation. It could be more than trillions of target molecule
for which we want to get the alchemical derivative. As a result,
we will be reloading the same wavefunction file trillions of times.
To improve it, we have to modify the source code such that
it can do a single reload of the wavefunction, than quickly modify
the atomic coordinate/pseudopotential parameters and gives a 
fast prediction. Effectively, we will be able to do high through put
calculations useing alchemy in CPMD

This folder contains the modified source file from unpatched CPMD-4.1
