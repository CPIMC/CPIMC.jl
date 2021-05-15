# Overview

This package provides a collection of types and routines for realizing path-integral Monte Carlo procedures in occupation number representation. In particular, one samples many-particle trajectories in imaginary time which are represented by the type `Configuration{T}`.

It is given the set of orbitals occupied at ``\tau = 0``, and a collection of one- and two-particle excitations and their respective imaginary times. The type used for representing single-particle basis states is arbitrary and corresponds to the type variable `T`.

One generally needs two kinds of functions:
- Updates, which propose changes to a `Configuration`, thus providing the building blocks of the Markov chain.
- Estimators, which give the contribution of a particular `Configuration` to an expectation value.

The Monte Carlo process is carried out by the function `sweep!`, which generates the Markov chain from which the expectation values are calculated.

A many-particle model is specified by the one- and two-particle matrix elements ``\epsilon_{ij}`` and ``w_{ijkl}`` with respect to the chosen basis. These can be used to calculate the weight changes for the acceptance probability.




