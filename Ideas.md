# Collection of Ideas for Concretization of New Programs.


* QM-Systems:
    * pure HEG - CPIMC
        * too large Kink-Weights for certain Parameters ← allow to rescale the Hamiltonian (i.e. Matrix elements) (automically)


* Procedure:
    * Option: Store full configuration instead of Samples: N "initial" occupation numbers and K kink-times and K 4-tuples (Kink-indices) for each sample vs. Samples for all observables. Could be less memory-consuming for low particle numbers N. Easier statistics over multiple program runs.

* Parallelization: Realize program with multiple parellel MCRuns. Thus all Samples (or Configurations) can be averaged together.

* Structure:
    * Results: use ABC to realize Canonical/GrandCanonical/... Ensembles.