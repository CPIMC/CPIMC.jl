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
Kinks
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
add_orbs
add_orbs!
drop_orbs
drop_orbs!
add_kinks
add_kinks!
remove_kinks
remove_kinks!
UpdateCounter
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
prev
next
prev_affecting
next_affecting
isunaffected
isunaffected_in_interval
```
