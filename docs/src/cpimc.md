# CPIMC.jl

```@meta
CurrentModule = CPIMC
```

```@docs
CPIMC
```

## sweep!, measure!, update!

```@docs
sweep!
measure!
update!
print_results
```

## Model

```@docs
Model
Orbital
```

## Configuration

```@docs
Configuration
ImgTime
Δ
Kink
T2
T4
basis
orbs
creators
annihilators
excite
excite!
occupations_at
next
prev
```

## weight calculation

```@docs
ΔWdiag_element
ΔW_diag
ΔWoffdiag_element
ΔT_element
Woffdiag_element
sign_offdiagonal_product
ladder_operator_order_factor
wminus
```

## constructing updates

```@docs
Step
apply_step
apply_step!
drop
drop!
add
add!
Update
```

## utility functions

```@docs
time_ordered_orbs
kinks_affecting_orbs
adjacent_kinks
adjacent_kinks_affecting_orbs
τ_prev_affecting
τ_next_affecting
τ_borders
isunaffected
isunaffected_in_interval
```
