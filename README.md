# CPIMC2020

New implementation of Configuration path-integral Monte Carlo (CPIMC) written in Julia.

Compared to other implementations we put particular emphasis on flexibility. The goal of this project is to provide an environment for trying out new ideas that is both easy to use and extend.

## Features
- full CPIMC simulation for the uniform electron gas using a plane waves basis
- no artificial restriction on the number of basis functions
- no explicit ordering of orbitals required
- multi-threading

## Getting started

A full description of the formalism can be found in Tim Schoof's PhD thesis[^1]. For a more modern perspective which is closer to this implementation one should also consider Kai Hunger's Master thesis [WIP].

As a gentle introduction we prepared a guided implementation of the ideal Fermi gas in `examples/ideal.jl` which should get you accustomed to the core functionality of this package. This should prepare you to explore on your own how the model of the uniform electron gas is realized in `src/UniformElectronGas.jl`.

## Documentation

Documentation and tutorials can be found [here](url).


## References

[^1]: Schoof, T. (2016). Configuration Path Integral Monte Carlo: Ab initio simulations of fermions in the warm dense matter regime [PhD thesis, CAU Kiel] [https://d-nb.info/1133492177/34](https://d-nb.info/1133492177/34)





