Alchemcy hack in Quantum Espresso
=================================
By default, 
restart files are written only for 'properly stopped' jobs.
To do first order alchemy, 
it is necessary to load the wavefunction of the reference system. 
Therefore it is necessary to change the source code 
such that the restart files will be written.

Moreover, ks energy will not be printed 
if the calculation is not converged.
Since we need to hijack the routine
to evaluate &lt;psi|H|psi>,
this behaviour also need to be changed.

The above two modification can be done simply modify the 
`/path/to/espresso/PW/src/electrons.f90` 
file with dirty hacks, which also change the overall behaviour. 
USE WITH CAUTION...
