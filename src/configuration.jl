export Configuration, Orbital, Kink, T2, T4, ImgTime, excite!, excite, Kinks, haskey, Kinks, keys, values, getindex

"""
Abstract type for single-particle basis states, implementation is required for each model.
"""
abstract type Orbital end

"""
Type alias for imaginary time.
`FixedPointNumbers` are used since these are stable for `==` and are thus stable as keys in `Dict`.
"""
const ImgTime = Fixed{Int64,60}


@doc raw"""
    T2{T}

Parametric type representing a 1-particle excitation by specifying
a transition from one state (annihilator) to another (creator).
In the occupation number representation, this reads

$a^{\dagger}_i a_j$

with a creator orbital `i` and an annihilator orbital `j`.
The single-particle basis is represented by the type parameter `T`.
"""
struct T2{T}
  " creator "
  i :: T
  " annihilator "
  j :: T
end

@doc raw"""
    T4{T}

Parametric type representing a 2-particle excitation by specifying
a transition from two state (annihilators) to two other (creator) states.
In the occupation number representation, this reads

$a^{\dagger}_i a^{\dagger}_j a_k a_l$

with creator orbitals `i`, `j` and an annihilator orbitals `k`, `l`.
The single-particle basis is represented by the type parameter `T`.
"""
struct T4{T}
  " creator "
  i :: T
  j :: T
  " annihilator "
  k :: T
  l :: T
end

"""
    const Kink{T} = Union{T2{T}, T4{T}}

Parametric type representing either a one- or a two-particle scattering event
given by types 'T2{T}' and 'T4{T}', respectively.
This term is motivated by the geometrical interpretation of a path
in Configuration Path Integral Monte Carlo formulation of the partition function
in occupation number representation.
Discrete occupations are given by horizontal lines of (quasi-)particles
occupying the corresponding single-particle states and where transitions between
these states are thus indicated by vertical 'kinks' connecting the orbital lines
that are part of the transition.
"""
const Kink{T} = Union{T2{T}, T4{T}}

# outer constructor method to construct a T2 kink, inferring the type parameter from the arguments
Kink(i,j) = T2(i,j)
# outer constructor method to construct a T4 kink, inferring the type parameter from the arguments
Kink(i,j,k,l) = T4(i,j,k,l)

# outer constructor method to extract a kink from a pair where the second element is a kink.
Kink(p::Pair{S,T} where {T<:Kink} where {S}) = p[2]# first substitute S, then T

"""
    const Kinks{T} = Vector{Pair{ImgTime, Kink{T}}}

Structure for storing excitations and their imaginary times (kinks).

Can also be constructed by passing such pairs directly.

    Kinks(pairs::Pair{ImgTime,<:Kink{T}}...)
"""
const Kinks{T} = Vector{Pair{ImgTime, Kink{T}}}

Kinks(pairs::Pair{ImgTime,<:Kink{T}}...) where {T}  = reduce(push!, pairs, init=Kinks{T}())

"""
    values(ck::Kinks)

Return a list of the excitations of a `Kinks`-object, used to allow the use of dictionary syntax.
"""
values(ck::Kinks) = last.(ck)

"""
    keys(ck::Kinks)

Return a list of the imaginary-times of a `Kinks`-object, used to allow the use of dictionary syntax.
"""
keys(ck::Kinks) = first.(ck)

"""
    Base.haskey(ck::Kinks, key::ImgTime)

Check if a `Kinks`-object contains a kink at a specific time.
"""
Base.haskey(ck::Kinks, key::ImgTime) = in(key, keys(ck))

"""
    Base.getindex(kinks::Kinks{T}, τ::ImgTime) where {T}

Allow the `kinks[τ]`-syntax to retrieve a `Kink` at a specific `ImgTime`.
"""
function Base.getindex(kinks::Kinks{T}, τ::ImgTime) where {T}
    searchsortedfirst(kinks, by=first, τ)
    τ_next, k = kinks[searchsortedfirst(kinks, by=first, τ)]
    if τ_next != τ
        Throw(KeyError(τ))
    else
        return k
    end
end


"""
Multi-particle trajectory in imaginary time.

`Configuration{T}` is a parametric type depending on the single-particle basis `T <: Orbital`.

**Fields**
- `occupations :: Set{T}`  -- orbitals which are occupied initially (`τ=0`)
- `kinks :: Kinks{T}`      -- collection of excitations and their imaginary times (kinks)

"""
mutable struct Configuration{T}
  occupations :: Set{T}
  kinks :: Kinks{T}
end

" outer constructor method for a configuration with occupations given by o and kinks given by k. "
Configuration(o::Set{T}, k) where {T} = Configuration(o, Kinks(k) )
" outer constructor method for a configuration with occupations given by o and no kinks. "
Configuration(o::Set{T}) where {T} = Configuration(o, Kinks{T}())
" outer constructor method for a configuriation with no occupations and kinks given by k. "
Configuration(k::Kinks{T}) where {T} = Configuration(Set{T}(), k)
" outer constructor method for a configuration with no occupations and kinks as given by varargs p... of Pair{ImgTime,<:Kink}, which are passed to the Kinks constructor "
Configuration(p::Pair{ImgTime,<:Kink{T}}...) where {T} = Configuration(Kinks(p...))
" outer constructor method for a configuration with occupations given by o and kinks as given by varargs p... of Pair{ImgTime,<:Kink}, which are passed to the Kinks constructor "
Configuration(o::Set{T}, p::Pair{ImgTime,<:Kink{T}}...) where {T} = Configuration(o, Kinks(p...))

" outer constructor method for empty Configurations{T} "
Configuration{T}() where T = Configuration(Set{T}())

"""
    basis(c::Configuration{T})

Return type parameter `T` of the configuration.
"""
basis(c::Configuration{T}) where T = T

"""
    orbs(::T2)

Return a Tuple of all orbitals that are affected by a `T2`-kink in the conventional ordering i, j.
"""
orbs(x::T2) = x.i, x.j

@doc raw"""
    orbs(::T4)

Return a Tuple of all orbitals that are affected by a `T4`-kink in the conventional ordering i, j, k, l.
This corresponds to the matrix element $w_{ijkl}$ in the same order.
"""
orbs(x::T4) = x.i, x.j, x.k, x.l

"""
    creators(::T4)

return a set of the two creators which are affected by a T4 kink
"""
creators(x::T4) = x.i, x.j

"""
    annihilators(x::T4)

return a set of the two annihilators which are affected by a T4 kink
"""
annihilators(x::T4) = x.k, x.l

"""
    ImgTime(t::Tuple{Nothing,Nothing})

Outer constructor method for a `Tuple{ImgTime,ImgTime}` from `Tuple{Nothing,Nothing}` which returns the bounds of the `ImgTime` interval as `Tuple (ImgTime(0), ImgTime(1))`
"""
ImgTime(t::Tuple{Nothing,Nothing}) = (ImgTime(0), ImgTime(1))

"""
    ImgTime(::Pair{ImgTime,<:Kink})

outer constructor method for an ImgTime from Pair{ImgTime,<:Kink}
this returns the first argument of the Pair
useful for automatic conversion
"""
ImgTime(p::Pair{ImgTime,<:Kink}) = first(p)


"""
    ImgTime(t::Tuple{Pair{ImgTime,<:Kink},Pair{ImgTime,<:Kink}})

Outer constructor method for a `Tuple{ImgTime,ImgTime}` from `Tuple{Pair{ImgTime,<:Kink},Pair{ImgTime,<:Kink}}` which returns the `ImgTime` of the Pairs as `Tuple{ImgTime,ImgTime}`
"""
ImgTime(t::Tuple{Pair{ImgTime,<:Kink},Pair{ImgTime,<:Kink}}) = (first(first(t)), first(last(t)))

@doc raw"""
    Δ(::ImgTime,::ImgTime)

return the length of the periodic interval between the two ::ImgTimes

the periodic distance of two imaginary times:

    `Δ(τ1,τ2) = τ2 - τ1` if `τ2 >= τ1` and
    `Δ(τ1,τ2) = 1 - (τ1 - τ2) = 1 + τ2 - τ1` else
"""
Δ(τ1::ImgTime,τ2::ImgTime) = τ1 > τ2 ? ImgTime(1) + τ2 - τ1 : τ2 - τ1

@doc raw"""
    excite(::Set{T}, ::T4{T})

Apply a kink to a set of basis states, i.e. return a set of basis states where
the states specified by the creators of the kink are added and the states specified
by the annihilators of the kink are dropped from the given set of basis states.
This has the physical meaning of a two-particle scattering event where two (quasi-)particles
in a many-body state change the single-particle states they occupy.
In occupation number representation this reads

$a^{\dagger}_i a^{\dagger}_j a_k a_l |\{n\}\rangle$

for creator orbitals `i` and `j` and annihilator orbitals `k` and `l`.
This function assumes fermionic particle statistics (Pauli principle),

$a^{\dagger}_i a^{\dagger}_i = a_i a_i = 0$

i.e. the target (creator) states must not be occupied and
the initial (annihilator) states must be occupied in the given set of states.
"""
function excite(o::Set{T}, κ::T4{T}) where T
  @assert ( in(κ.k, o) & in(κ.l, o) ) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the annihilators $(κ.k), $(κ.l) is not occupied."
  @assert ( !in(κ.i, o) & !in(κ.j, o) ) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the creators $(κ.i), $(κ.j) is already occupied. (Pauli principle)"
  union(setdiff(o, Set([κ.k,κ.l])), Set([κ.i, κ.j]))
end

"""
    excite(o::Set{T}, κ::Pair{ImgTime,T4{T}})

Apply a `T4`-kink to a set of basis states for a pair of a time and a `T4`-kink. This is useful for iterating over a `Kinks`-object when calculating the occupations at a specific time.
"""
excite(o::Set{T}, κ::Pair{Fixed{Int64,60},Union{T2{T}, T4{T}}}) where T = excite(o, last(κ))
excite(o::Set{T}, κ::Pair{ImgTime,T4{T}}) where T = excite(o, last(κ))

"""
    excite!(::Set{T}, i::T, j::T, k::T, l::T) where {T}
    excite!(o::Set{T}, κ::T4{T})

Apply an excitation in-place to a set of basis states that is given by creating orbitals i, j
and annihilating the orbitals k, l.
"""
function excite!(o::Set{T}, i::T, j::T, k::T, l::T) where {T}
    @assert (in(k, o) & in(l, o)) "Kink ($(i),$(j),$(k),$(l)) cannot be applied: one or two of the annihilators $(k), $(l) is not occupied. (Pauli-Principle)"
    @assert (!in(i, o) & !in(j, o)) "Kink ($(i),$(j),$(k),$(l)) cannot be applied: one or two of the creators $(i), $(j) is already occupied. (Pauli-Principle)"
    delete!(o, k)
    delete!(o, l)
    push!(o, i)
    push!(o, j)
end

excite!(o::Set{T}, κ::T4{T}) where T = excite!(o, κ.i, κ.j, κ.k, κ.l)


"""
    occupations_at(o::Set{T}, kinks::Kinks{T})

Return the occupied orbitals after applying all kinks to initial occupation.
"""
function occupations_at(o::Set{T}, kinks::Kinks{T}) :: Set{T} where {T}
  foldl(excite, kinks; init=o)
end

"""
    occupations_at(c::Configuration, τ::ImgTime)

Return the occupied orbitals to the right of τ.
"""
function occupations_at(c::Configuration, τ::ImgTime)
  occupations_at(c.occupations, filter(x -> first(x) <= τ, c.kinks))
end

"""
    next(x, τ)

Return first kink after to the given `τ`. If there is no kink with `::ImgTime` larger than `τ`, return first kink.
"""
function next(x, τ)
    if isempty(x)
        throw(DomainError(" prev(::Excitation_arr, τ) is not supported for empty x. It makes no sense to get the element previous to some element for an empty collection. "))
    end
    index = findfirst(p -> first(p) > τ, x)
    if isnothing(index)# if no entry is found with ImgTime < τ
        x[begin]# return last entry (assuming x is sorted with respect to ImgTime)
    else
        x[index]
    end
end


"""
    prev(x, τ)

Return first kink earlier than `τ::ImgTime`. If there is no such kink the last kink is returned.
If `ck` is empty this will throw an error.

This function assumes that the kinks in `x` are ordered by their imaginary time!
"""
function prev(x, τ)
    if isempty(x)
        throw(DomainError(" prev(::Excitation_arr, τ) is not supported for empty x. It makes no sense to get the element previous to some element for an empty collection. "))
    end
    index = findlast(p -> first(p) < τ, x)
    if isnothing(index)# if no entry is found with ImgTime < τ
        x[end]# return last entry (assuming x is sorted with respect to ImgTime)
    else
        x[index]
    end
end


"""
     kinks_affecting_orbs(itr::Any, itr::Any)

Return all kinks in the first argument that affect one ore more of the orbitals the second argument.
The first itr is expected to contain tuples or at least types on which `last()` is defined
and `last(i)` is intended to contain a ::T4 for all elements `i ∈ itr` or at least a type
which has fields with the names `i`, `j`, `k`, `l`.
"""
kinks_affecting_orbs(ck, os) = filter( p -> any( Set([last(p).i,last(p).j,last(p).k,last(p).l]) .∈ (os,) ), ck)
kinks_affecting_orbs(c::Configuration, os) = kinks_affecting_orbs(c.kinks, os)


"""
    adjacent_kinks(itr::Any, ::Any)

Return a tuple of the kink in `itr`
that is closest right of the time in the second argument.
This expects that methods for the functions `next()` and `prev()` are defined
for the given argument types.
Used for getting a tuple of the neighbouring kinks from some imaginary time.
"""
adjacent_kinks(ck::Any, τ::Any) = prev(ck, τ), next(ck, τ)


"""
    adjacent_kinks_affecting_orbs(::Any, ::Any, ::Any)

Return a tuple of the closest kink to the right and the closest kink to the left of some imaginary time
that affect one of the orbitals in some collection. If there are no such kinks, return (nothing,nothing).
The first argument is expected to be an iterable which contains pairs as elements, with
the first element of each pair containing an imaginary time and the second element of each pair containing a kink.
The second argument is expected to be a collection of orbitals, that are to be considered regarding affection with one of the
kinks in the first argument.
The third argument is expected to be an imaginary time from which the neighbouring kinks affected by any orbitals from the
second argument are to be determined."""
function adjacent_kinks_affecting_orbs(ck, os, τ)
    k = kinks_affecting_orbs(ck, os)
    if isempty(k)
        return (nothing,nothing)
    else
        return adjacent_kinks(k, τ)
    end
end
adjacent_kinks_affecting_orbs(c::Configuration, os, τ) = adjacent_kinks_affecting_orbs(c.kinks, os, τ)


"""
     τ_prev_affecting(::Any, ::Any, ::Any)

Return the ImgTime of the closest kink to the left of τ
that affects one of the orbitals in os.
If no orbital in os is affected by and kink,
return the lower interval bound ImgTime(0).
"""
function τ_prev_affecting(ck, os, τ)
    κs = kinks_affecting_orbs(ck, os)
    ImgTime(prev(κs, τ))
end


"""
     τ_next_affecting(::Any, ::Any, ::Any)

Return the ImgTime of the closest kink to the right of τ
that affects one of the orbitals in os.
If no orbital in os is affected by and kink,
return the lower interval bound ImgTime(0).
"""
function τ_next_affecting(ck, os, τ)
    κs = kinks_affecting_orbs(ck, os)
    ImgTime(next(κs, τ))
end

"""
     τ_borders(::Kinks{T}, orbs, ::ImgTime) where {T <: Orbital}

Return a tuple of
the ImgTime of the closest kink to the right and
the ImgTime of the closest kink to the left of τ
that affect one of the orbitals in orbs.
If no orbital in os is affected by and kink from the collection in the first argument,
return a tuple of the interval bounds (ImgTime(0), ImgTime(1))."""
τ_borders(ck::Kinks{T}, orbs, τ::ImgTime) where {T <: Orbital} = ImgTime(adjacent_kinks_affecting_orbs(ck, orbs, τ))

τ_borders(c::Configuration{T}, orbs, τ::ImgTime) where {T <: Orbital} = τ_borders(c.kinks, orbs, τ)


"""
    isunaffected(Kinks, orb)

Return if an orbital is not affected by any kink.
"""
isunaffected(kinks, orb) = all(kink -> orb ∉ last(kink), kinks)

Base.in(i, k::T2) = i == k.i || i == k.j
Base.in(i, k::T4) = i == k.i || i == k.j || i == k.k || i == k.l

"""
    in_open_interval(τ, τ_first, τ_last)

Check if `τ` lies in the open `ImgTime`-interval `(τ_first,`τ_last)`, assuming that `τ_first` is the left border and `τ_last` the right one.
"""
in_open_interval(τ, τ_first, τ_last) = (τ_first != τ != τ_last != τ_first) & ((τ_first < τ_last) == ((τ < τ_first) != (τ < τ_last)))

"""
    isunaffected_in_interval(kinks, orb, τ_first::ImgTime, τ_last::ImgTime)

Return if an orbital is not affected by any of the kinks from `ck` in the open interval `(τ_first,τ_last)`.
"""
isunaffected_in_interval(kinks, orb, τ_first::ImgTime, τ_last::ImgTime) = all(kink -> ((orb ∉ last(kink)) && in_open_interval(first(kink), τ_first, τ_last) ), kinks)

## Method definitions for function drop
# This function is mostly used for calculating the changes proposed by an update
# in order to calculate proposal and/or acceptance probabilities.

"""
    drop(c::Configuration, n::Nothing)

Return the configuration c without dropping anything.
"""
drop(c::Configuration, n::Nothing) = c

"""
    drop(oc::Set{T}, o::T) where {T <: Orbital}

Return a set with the orbital o dropped from oc.
"""
drop(oc::Set{T}, o::T) where {T <: Orbital} = setdiff(oc, Set([o]))

"""
    drop(c::Configuration{T}, o::T) where {T <: Orbital}

Return a configuration with the `Orbital` `o` dropped from `c.occupations`.
"""
drop(c::Configuration{T}, o::T) where {T <: Orbital} = Configuration(drop(c.occupations, Set([o])),c.kinks)

"""
    drop(oc1::Set{T}, oc2::Set{T}) where {T <: Orbital}

Return a `Set` with the orbitals in `oc2` dropped from `oc1`.
"""
drop(oc1::Set{T}, oc2::Set{T}) where {T <: Orbital} = setdiff(oc1, oc2)

"""
    drop(oc1::Set, oc2::Tuple)

Return a `Set` with the orbitals in `oc2` dropped from `oc1`.
"""
drop(oc1::Set, oc2::Tuple) = setdiff(oc1, oc2)

"""
    drop(c::Configuration{T}, oc::Set{T}) where {T <: Orbital}

Return a `Configuration` with the orbitals in `oc` dropped from `c.occupations`.
"""
drop(c::Configuration{T}, oc::Set{T}) where {T <: Orbital} = Configuration(drop(c.occupations, oc), c.kinks)

"""
    drop(ck1, ck2)

Return a container with the elements in `ck2` dropped from `ck1`.
"""
drop(ck1, ck2) = setdiff(ck1, ck2)

"""
    drop(c::Configuration{T}, ck::Kinks)

Return a `Configuration` with the pairs in `ck` dropped from `c.kinks`.
"""
drop(c::Configuration{T}, ck::Kinks) where {T <: Orbital} = Configuration(c.occupations, drop(c.kinks, ck))

"""
    drop(ck::Kinks, ps::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital}

Return a `Kinks`-object with the pairs `ps...` dropped from `ck`.
"""
drop(ck::Kinks, ps::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital} = setdiff(ck, Kinks(ps))

"""
    drop(c::Configuration{T}, ps::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital}

Return a `Configuration` with the pairs `ps...` dropped from `c.kinks`.
"""
drop(c::Configuration{T}, ps::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital} = Configuration(c.occupations, drop(c.kinks, ps...))

"""
    drop(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}

Return a `Configuration` with occupations and kinks in `c2` dropped from `c1`.
"""
drop(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital} = Configuration(drop(c1.occupations,c2.occupations), drop(c1.kinks,c2.kinks))


## Method definitions for function drop!
# This function is mostly for applying the changes determined in an MC Step Δ
# to the current configuration c in function apply_step!(c, Δ). cf. CPIMC.jl
# These methods change the first argument in-place.

"""
    drop!(c::Configuration, n::Nothing)

Drop nothing.
"""
function drop!(c::Configuration, n::Nothing)
    nothing
end

"""
    drop!(c::Configuration{T}, o::T) where {T <: Orbital}

Drop a single orbital.
"""
function drop!(c::Configuration{T}, o::T) where {T <: Orbital}
    delete!(c.occupations, o)
end

"""
    drop!(c::Configuration{T}, oc::Set{T}) where {T <: Orbital}

Drop all orbitals given in a set `oc` from `c.occupations`.
"""
function drop!(c::Configuration{T}, oc::Set{T}) where {T <: Orbital}
    for o in oc
        drop!(c, o)
    end
end


"""
    drop!(c::Configuration, τ::ImgTime)

Drop a single kink at time τ from `c.kinks`.
"""
function drop!(c::Configuration, τ::ImgTime)
    setdiff!(c.kinks, [τ => c.kinks[τ]])
end

"""
    drop!(c::Configuration{T}, ck::Kinks) where {T <: Orbital}

Drop all kinks given in `ck` from `c.kinks`.
"""
function drop!(c::Configuration{T}, ck::Kinks) where {T <: Orbital}
    for (τ, k) in ck
        drop!(c, τ)
    end
end

"""
    drop!(c::Configuration{T}, pk::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital}

Drop all kinks given by `Varargs{Pair{ImgTime,<:Kink{T}}}` from `c.kinks`.
"""
function drop!(c::Configuration{T}, pk::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital}
    for (τ, k) in pk
        drop!(c, τ)
    end
end

"""
    drop!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}

Drop occupations and kinks given by second argument `c2` from first argument `c1`.
"""
function drop!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}
    drop!(c1, c2.occupations)
    drop!(c1, c2.kinks)
end


## Method definitions for function add.
# This function is mostly used for calculating the changes proposed by an update
# in order to calculate proposal and/or acceptance probabilities.

"""
    add(c::Configuration, n::Nothing)

Return the `Configuration` `c` without adding anything.
"""
add(c::Configuration, n::Nothing) = c

"""
    add(oc::Set{T}, o::T) where {T <: Orbital}

Return a set with the orbital o added to oc.
"""
add(oc::Set{T}, o::T) where {T <: Orbital} = union(oc, Set([o]))

"""
    add(c::Configuration{T}, o::T) where {T <: Orbital}

Return a configuration with the orbital o added to c.occupations.
"""
add(c::Configuration{T}, o::T) where {T <: Orbital} = Configuration(add(c.occupations, o), c.kinks)

"""
    add(oc1::Set{T}, oc2::Set{T}) where {T <: Orbital}

Return a set with the orbitals in oc2 added to oc1.
"""
add(oc1::Set{T}, oc2::Set{T}) where {T <: Orbital} = union(oc1, oc2)

"""
    add(c::Configuration{T}, oc::Set{T}) where {T <: Orbital}

Return a configuration with the orbitals in oc added to c.occupations.
"""
add(c::Configuration{T}, oc::Set{T}) where {T <: Orbital} = Configuration(union(c.occupations, oc), c.kinks)


"""
    add(ck::Kinks, p::Pair)

Return a `Kinks`-object with the pair `p` added to `ck` with respect to the sorting.
"""
add(ck::Kinks, p::Pair) = insert(ck, searchsortedfirst(ck, by=first, first(p)), p)


"""
    add(ck::Kinks, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital}

Return a `Kinks`-object with the pairs `ps...` added to `ck`.
"""
function add(ck::Kinks, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital}
    ck_copy = copy(ck)
    for k in ps
        add!(ck_copy, k)
    end
    return ck_copy
end

"""
    add(ck1::Kinks, ck2::Kinks)

Return a `Kinks`-object with the pairs from `ck2` added to `ck1`.
"""
function add(ck1::Kinks, ck2::Kinks) where {T <: Orbital}
    ck1_copy = copy(ck1)
    for k in ck2
        add!(ck1_copy, k)
    end
    return ck1_copy
end

"""
    add(c::Configuration{T}, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital}

Return a configuration with the pairs ps... added to `c.kinks`.
"""
add(c::Configuration{T}, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital} = Configuration(c.occupations, add(c.kinks, ps...))

"""
    add(c::Configuration{T}, ck::Kinks) where {T <: Orbital}

Return a configuration with the pairs in ck added to `c.kinks`.
"""
add(c::Configuration{T}, ck::Kinks) where {T <: Orbital} = add(c, ck...)

"""
    add(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}

Return a configuration with occupations and kinks in c2 dropped from c1.
"""
add(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital} = Configuration(add(c1.occupations, c2.occupations), add(c1.kinks, c2.kinks))


## Method definitions for function add!.
# This function is mostly for applying the changes determined in an MC Step Δ
# to the current configuration c in function apply_step!(c, Δ). cf. CPIMC.jl
# These methods change the first argument in-place.


"""
    add!(c::Configuration, n::Nothing)

Add nothing.
"""
function add!(c::Configuration, n::Nothing)
    nothing
end

"""
    add!(c::Configuration{T}, o::T) where {T <: Orbital}

Add a single orbital.
"""
function add!(c::Configuration{T}, o::T) where {T <: Orbital}
    push!(c.occupations, o)
end

"""
    add!(c::Configuration{T}, oc::Set{T}) where {T <: Orbital}

Add all orbitals given in a set oc to c.occupations.
"""
function add!(c::Configuration{T}, oc::Set{T}) where {T <: Orbital}
    union!(c.occupations, oc)
end

"""
    add!(ck::Kinks, τ::ImgTime, k::Kink)

Add the pair `(τ, k)` to the `ck` with respect to the sorting.
"""
add!(ck::Kinks, τ::ImgTime, k::Kink) = insert!(ck, searchsortedfirst(ck, by=first, τ), τ => k)

function Base.setindex!(ck::Kinks, kink::Kink, τ::ImgTime)
    add!(ck::Kinks, τ::ImgTime, kink::Kink)
end

"""
  add!(c::Configuration{T}, p::Pair ) where {T <: Orbital}

Add a single kink at time first(p) to `c.kinks`.
"""
function add!(c::Configuration{T}, p::Pair ) where {T <: Orbital}
    add!(c.kinks, first(p), last(p))
end



"""
    add!(ck::Kinks, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital}

Add all kinks given by the pairs `ps` to `ck`.
"""
function add!(ck::Kinks, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital}
    for k in ps
        add!(ck, first(k), last(k))
    end
    return ck
end

"""
    add!(c::Configuration{T}, ck::Kinks) where {T <: Orbital}

Add all kinks given in `ck` to `c.kinks`.
"""
function add!(c::Configuration{T}, ck::Kinks) where {T <: Orbital}
    add!(c.kinks, ck...)
end

"""
    add!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}

Add occupations and kinks given by second argument `c2` from first argument `c1`.
"""
function add!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}
    add!(c1, c2.occupations)
    add!(c1, c2.kinks)
end

"""
    time_ordered_orbs(ck::Kinks{T}) where {T <: Orbital}

Get a list of orbitals that affect each kink in the time-ordering of the kinks and in the conventional ordering i, j, k, l.
"""
function time_ordered_orbs(ck::Kinks{T}) where {T <: Orbital}
    if isempty(ck)
        return Array{T,1}()
    else
        return collect( Iterators.flatten( orbs(k) for k in values(ck) ) )
    end
end

@doc raw"""
    ladder_operator_order_factor(sortedlist::Array{<:Orbital,1})

Returns $1$ or $-1$ depending on the order of all ladder operators as given by a list of orbitals.
"""
function ladder_operator_order_factor(sortedlist::Array{<:Orbital,1})
    phase_factor = 1
    while !isempty(sortedlist)
        index = 1
        op_orb = sortedlist[index]
        while sortedlist[index + 1] != op_orb
            sortedlist[index:index+1] = [sortedlist[index+1], sortedlist[index]]
            phase_factor *= -1
            index += 1
            @assert(index <= length(sortedlist))
        end
        deleteat!(sortedlist, [index, index+1])
    end
    return phase_factor
end


@doc raw"""
    ladder_operator_order_factor(ck::Kinks)

Returns $1$ or $-1$ depending on the order of all ladder operators as given by the orbitals that affect each kink in the time-ordering of the kinks and in the conventional ordering i, j, k, l,
used in the sign estimator.
"""
ladder_operator_order_factor(ck::Kinks) = ladder_operator_order_factor(time_ordered_orbs(ck))


"""
    is_type_1(left_kink::T4, right_kink::T4)

Returns True if left_kink and right_kink are entangled in a Type-1 way.
This does not check wether the two kinks are neighbouring."""
function is_type_1(left_kink::T4, right_kink::T4)
  if ((length(intersect(Set([left_kink.i, left_kink.j]),Set([right_kink.k, right_kink.l]))) == 1) &
        (length(intersect(Set([left_kink.k, left_kink.l]), Set([right_kink.i, right_kink.j]))) == 0)) ||
    ((length(intersect(Set([left_kink.i, left_kink.j]),Set([right_kink.k, right_kink.l]))) == 0) &
          (length(intersect(Set([left_kink.k, left_kink.l]), Set([right_kink.i, right_kink.j]))) == 1))
    return true
  else
    return false
  end
end

"""
    right_type_1_chain_length(ck, τ, count = 0)

Returns the length of the chain of type-1-entaglements starting with the Kink at τ counting to the right.
"""
function right_type_1_chain_length(ck::Kinks, τ, counted_τs = [])
    kink = ck[τ]
    next_kink = next(kinks_affecting_orbs(ck, Set([kink.i, kink.j, kink.k, kink.l])), τ)
    if is_type_1(ck[τ], last(next_kink)) & !in(τ, counted_τs)
        push!(counted_τs,τ)
        return right_type_1_chain_length(ck, first(next_kink), counted_τs)
    else
        return length(counted_τs)
    end
end

"""
    left_type_1_chain_length(ck::Kinks, τ, counted_τs = [])

Returns the length of the chain of type-1-entaglements starting with the Kink at τ counting to the left.
"""
function left_type_1_chain_length(ck::Kinks, τ, counted_τs = [])
    kink = ck[τ]
    prev_kink = prev(kinks_affecting_orbs(ck, Set([kink.i, kink.j, kink.k, kink.l])), τ)
    if is_type_1(last(prev_kink), ck[τ]) & !in(τ, counted_τs)
        push!(counted_τs,τ)
        return left_type_1_chain_length(ck, first(prev_kink), counted_τs)
    else
        return length(counted_τs)
    end
end

"""
    longest_type_1_chain_length(ck::Kinks)

Returns the longest chain of type-1-entaglements in ck.
"""
function longest_type_1_chain_length(ck::Kinks)
    longest_length = 0
    for (τ, kink) in ck
        longest_length = max(right_type_1_chain_length(ck, τ),
                                left_type_1_chain_length(ck, τ), longest_length)
    end
    return longest_length
end

"""
    right_type_1_count(ck::Kinks)

Returns the number of kinks that are type 1 entagled with their right neighbour.
"""
function right_type_1_count(ck::Kinks)
    count = 0
    for (τ, kink) in ck
        if is_type_1(kink, last(next(kinks_affecting_orbs(ck, Set([kink.i, kink.j, kink.k, kink.l])),τ)))
            count += 1
        end
    end
    return count
end

"""
    random_shuffle(::Any, ::Any)

shuffle or not shuffle the given tuple with equal propability:
Return either the same tuple, or the tuple in reverse order,
with equal probability 0.5.
"""
function random_shuffle(i, j)
    if rand() < 0.5
        return i, j
    else
        return j, i
    end
end


"""
    shuffle_annihilators(::Any, ::Any, ::Any, ::Any)

Return a tuple where the last two arguments are randomly shuffled with equal probability 0.5
in comparison to the ordering of the input.
"""
shuffle_annihilators(i, j, k, l) = i, j, random_shuffle(k, l)...
shuffle_annihilators(κ::T4) = T4(shuffle_annihilators(κ.i, κ.j, κ.k, κ.l)...)

"""
    shuffle_creators(::Any, ::Any, ::Any, ::Any)

Return a tuple where the first two arguments are randomly shuffled with equal probability 0.5
in comparison to the ordering of the input.
"""
shuffle_creators(i, j, k, l) = random_shuffle(i, j)..., k, l
shuffle_creators(κ::T4) = T4(shuffle_creators(κ.i, κ.j, κ.k, κ.l)...)

"""
    shuffle_indices(::Any, ::Any, ::Any, ::Any)

Return a tuple where both the first two arguments and the last two arguments
are randomly shuffled with equal probability 0.5 respectively
in comparison to the ordering of the input.
"""
shuffle_indices(i, j, k, l) = random_shuffle(i, j)..., random_shuffle(k, l)...
shuffle_indices(κ::T4) = T4(shuffle_indices(κ.i, κ.j, κ.k, κ.l)...)


"""
    kinks_from_periodic_interval(::Kinks, τ1, τ2)

return kinks with τ ∈ (τ1,τ2) if τ1 < τ2 and τ ∈ (τ2,1) ∪ (0,τ1) if τ1 > τ2
"""
function kinks_from_periodic_interval(ck::Kinks, τ1, τ2)
    filter(x-> in_open_interval(first(x), τ1, τ2), ck)
end

"""
    times_from_periodic_interval(::Kinks, ::ImgTime, ::ImgTime)

return a list of all times of kinks with τ ∈ (τ1,τ2) if τ1 < τ2 or τ ∈ (τ2,1) ∪ (0,τ1) if τ1 > τ2
in the periodic ordering suggested by the relation of the first time-argument τ1 to the second time-argument τ2
"""
function times_from_periodic_interval(ck::Kinks, τ1::ImgTime, τ2::ImgTime)
    keys(filter(x-> in_open_interval(first(x), τ1, τ2), ck))
end
