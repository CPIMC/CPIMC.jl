export Configuration, Orbital, Kink, T2, T4, ImgTime, excite!, excite, Kinks, haskey, Kinks, times, excitations, keys, values, getindex, (==)

"""
Abstract type for single-particle basis states, implementation is required for each model.
"""
abstract type Orbital end

"""
Type alias for imaginary time.
`FixedPointNumbers` are used since these are stable for `==` and are thus stable as times in `Dict`.
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
    excitations(ck::Kinks)

Return a list of the excitations of a `Kinks`-object, used to allow the use of dictionary syntax.
"""
excitations(ck::Kinks) = last.(ck)


"""
    times(ck::Kinks)

Return a list of the imaginary-times of a `Kinks`-object, used to allow the use of dictionary syntax.
"""
times(ck::Kinks) = first.(ck)

"""
    Base.haskey(ck::Kinks, key::ImgTime)

Check if a `Kinks`-object contains a kink at a specific time.
"""
Base.haskey(ck::Kinks, key::ImgTime) = in(key, times(ck))

"""
    Base.getindex(kinks::Kinks{T}, ??::ImgTime) where {T}

Allow the `kinks[??]`-syntax to retrieve a `Kink` at a specific `ImgTime`.
"""
function Base.getindex(kinks::Kinks{T}, ??::ImgTime) where {T}
    ??_next, k = kinks[searchsortedfirst(kinks, by=first, ??)]
    if ??_next != ??
        throw(KeyError(??))
    else
        return k
    end
end


function Base.setindex!(ck::Kinks, kink::Kink, ??::ImgTime)
    add!(ck::Kinks, ??::ImgTime, kink::Kink)
end

"""
Multi-particle trajectory in imaginary time.

`Configuration{T}` is a parametric type depending on the single-particle basis `T <: Orbital`.

**Fields**
- `occupations :: Set{T}`  -- orbitals which are occupied initially (`??=0`)
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

Base.:(==)(c1::Configuration, c2::Configuration) = (c1.occupations == c2.occupations) && (c1.kinks == c2.kinks)

"""
    drop_orbs(occ, orbs)

Return a copy of the occupation `occ` without the orbitals in `orbs`.
"""
drop_orbs(occ, ::Nothing) = occ
drop_orbs(occ, orbs) = filter(o -> !in(o,orbs), occ)

"""
    add_orbs(occ, orbs)

Return a copy of the occupation `occ` with the orbitals in `orbs` added.
"""
add_orbs(occ, ::Nothing) = occ
add_orbs(occ, orbs) = union(occ, orbs)

"""
    drop_kinks(ck, s)

Return a copy of the kinks in `ck` without the kinks in `s`.
"""
drop_kinks(ck, ::Nothing) = ck
drop_kinks(ck, s) = filter(k -> !in(k,s), ck)

"""
    add_kinks(ck::Kinks, s)

Return a copy of the kinks in `ck` with the kinks in `s` added.
"""
add_kinks(ck::Kinks, ::Nothing) = ck

function add_kinks(ck::Kinks, s)
    ck_copy = deepcopy(ck)
    add_kinks!(ck_copy,s)
    return ck_copy
end

"""
    drop_orbs!(occ, s)

Remove the orbitals in `s` from the occupation `occ`.
"""
drop_orbs!(occ, ::Nothing) = nothing
drop_orbs!(occ, s) = filter!(o -> !in(o,s), occ)

"""
    add_orbs!(occ, s)

Add the orbitals in `s` to the occupation `occ`.
"""
add_orbs!(occ, ::Nothing) = nothing
add_orbs!(occ, s) = union!(occ, s)

"""
   drop_kinks!(ck, s)

Remove the kinks in `s` from `ck`.
"""
drop_kinks!(ck, ::Nothing) = nothing
drop_kinks!(ck, s) = filter!(k -> !in(k,s), ck)

"""
   add_kinks!(ck::Kinks, s)

Add the kinks in `s` to `ck`.
"""
add_kinks!(ck::Kinks, ::Nothing) = nothing
function add_kinks!(ck::Kinks, s)
    for k in s
        insert!(ck, searchsortedfirst(ck, by=first, first(k)), first(k) => last(k))
    end
end


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
    ??(::ImgTime,::ImgTime)

return the length of the periodic interval between the two ::ImgTimes

the periodic distance of two imaginary times:

    `??(??1,??2) = ??2 - ??1` if `??2 >= ??1` and
    `??(??1,??2) = 1 - (??1 - ??2) = 1 + ??2 - ??1` else
"""
??(??1::ImgTime,??2::ImgTime) = ??1 > ??2 ? ImgTime(1) + ??2 - ??1 : ??2 - ??1

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
function excite(o::Set{T}, ??::T4{T}) where T
  @assert ( in(??.k, o) & in(??.l, o) ) "Kink ($(??.i),$(??.j),$(??.k),$(??.l)) cannot be applied: one or two of the annihilators $(??.k), $(??.l) is not occupied."
  @assert ( !in(??.i, o) & !in(??.j, o) ) "Kink ($(??.i),$(??.j),$(??.k),$(??.l)) cannot be applied: one or two of the creators $(??.i), $(??.j) is already occupied. (Pauli principle)"
  union(setdiff(o, Set([??.k,??.l])), Set([??.i, ??.j]))
end

"""
    excite(o::Set{T}, ??::Pair{ImgTime,T4{T}})

Apply a `T4`-kink to a set of basis states for a pair of a time and a `T4`-kink. This is useful for iterating over a `Kinks`-object when calculating the occupations at a specific time.
"""
excite(o::Set{T}, ??::Pair{Fixed{Int64,60},Union{T2{T}, T4{T}}}) where T = excite(o, last(??))
excite(o::Set{T}, ??::Pair{ImgTime,T4{T}}) where T = excite(o, last(??))

"""
    excite!(::Set{T}, i::T, j::T, k::T, l::T) where {T}
    excite!(o::Set{T}, ??::T4{T})

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

excite!(o::Set{T}, ??::T4{T}) where T = excite!(o, ??.i, ??.j, ??.k, ??.l)


"""
    occupations_at(o::Set{T}, kinks::Kinks{T})

Return the occupied orbitals after applying all kinks to initial occupation.
"""
function occupations_at(o::Set{T}, kinks::Kinks{T}) :: Set{T} where {T}
  foldl(excite, kinks; init=o)
end

"""
    occupations_at(c::Configuration, ??::ImgTime)

Return the occupied orbitals to the right of ??.
"""
function occupations_at(c::Configuration, ??::ImgTime)
  occupations_at(c.occupations, filter(x -> first(x) <= ??, c.kinks))
end

"""
    next(a, ??)

Return index of the first kink after to the given `??`. If there is no kink with `::ImgTime` larger than `??`, return index of first kink.
Throw a `DomainError` if an empty list is passed as first argument.

This function assumes that the kinks in `a` are ordered with respect to their imaginary time!
"""
function next(a, ??)
    if isempty(a)
        throw(DomainError(" next(x::Kinks, ??) is not supported for empty x. It makes no sense to get the element next to some element for an empty collection. "))
    end
    index = searchsortedfirst(a, ??, by=first) # find first element x ??? a with first(x) >= ??
    if index == length(a) + 1 # if there is no element x ??? a with first(x) >= ??
        return 1#  return index of the first element in a
    end
    # catch equality:
    # if there there is a kink at time ??, `searchsortedfirst` returns the index of that kink
    if ?? == first(a[index])
        index += 1 # increase index by one
    end
    # catch if kink at ?? is the last kink in a:
    # then the index after the index of the kink at ?? is length(a) + 1
    if index == length(a) + 1
        return 1 # return index of the first element
    else
        return index # return the index
    end
end

"""
     next_affecting(a::Any, orbs::Any, ??::Any)

Return the index of the closest kink to the right of ??
that affects one of the orbitals in orbs.
If no orbital in orbs is affected by any kink return 0.
"""
function next_affecting(a, orbs, ??)
    if isempty(a)
        return 0
    end
    i = next(a, ??)
    if any( orbs .??? ( last(a[i]) ,) ) # if any orbs affected by the kink at i
        return i
    end
    # otherwise, iterate forward periodically until first affecting kink is reached
    for _ = 1:length(a) - 1 # iterate until before the same kink it reached. note: this should not occur if o is affected by any kink since there should always be at least two kinks for each affected orbital due to periodicity (i.e. all creators have to pair up with a corresponding annihilator)
        i += 1 # increment
        if i == length(a) + 1 # periodicity: if last index is reached, start from beginning of collection
            i = 1
        end
        if any( orbs .??? ( last(a[i]) ,) )# if any orbs are affected by the kink at i
            return i # if o is affected by the kink at i
        end
    end
    return 0
end


"""
    prev(a, ??)

Return the index of the first kink earlier than `??::ImgTime`. If there is no such kink the last kink is returned.
Throw a `DomainError` if an empty list is passed as first argument.

This function assumes that the kinks in `x` are ordered with respect to their imaginary time!
"""
function prev(a, ??)
    if isempty(a)
        throw(DomainError(" prev(x::Kinks, ??) is not supported for empty x. It makes no sense to get the element previous to some element for an empty collection. "))
    end
    index = searchsortedlast(a, ??, by=first) # find last element x ??? a with first(x) <= ??
    if index == 0 # if there is no element x ??? a with first(x) <= ??
        return length(a)#  return index of the last element in a
    end
    # catch equality:
    # if there there is a kink at time ??, `searchsortedlast` returns the index of that kink
    if ?? == first(a[index])
        index -= 1 # reduce index by one
    end
    # catch if kink at ?? is the first kink in a:
    # then the index before the index of the kink at ?? is 0
    if index == 0
        return length(a) # return index of the last element
    else
        return index # return the index
    end
end

"""
     prev_affecting(a::Any, orbs::Any, ??::Any)

Return the index of the closest kink to the left of ??
that affects one of the orbitals in orbs.
If no orbital in orbs is affected by any kink return 0.
"""
function prev_affecting(a, orbs, ??)
    if isempty(a)
        return 0
    end
    i = prev(a, ??)
    if any( orbs .??? ( last(a[i]), ) ) # if any orbs are affected by the kink at i
        return i
    end
    # otherwise, iterate backwards periodically until first affecting kink is reached
    for _ = 1:length(a) - 1 # iterate until before the same kink it reached. note: this should not occur if o is affected by any kink since there should always be at least two kinks for each affected orbital due to periodicity (i.e. all creators have to pair up with a corresponding annihilator)
        i -= 1 # decrement
        if i == 0 # periodicity: if first index is reached, start from end of collection
            i = length(a)
        end
        if any( orbs .??? ( last(a[i]), ) ) # if any orbs are affected by the kink at i
            return i
        end
    end
    return 0
end


"""
     kinks_affecting_orbs(itr::Any, itr::Any)

Return all kinks in the first argument that affect one ore more of the orbitals the second argument.
The first itr is expected to contain tuples or at least types on which `last()` is defined
and `last(i)` is intended to contain a ::T4 for all elements `i ??? itr` or at least a type
which has fields with the names `i`, `j`, `k`, `l`.
"""
kinks_affecting_orbs(ck, os) = filter( p -> any( orbs(last(p)) .??? (os,) ), ck)
kinks_affecting_orbs(c::Configuration, os) = kinks_affecting_orbs(c.kinks, os)


"""
    adjacent_kinks(itr::Any, ::Any)

Return a tuple of the kink in `itr`
that is closest right of the time in the second argument.
This expects that methods for the functions `next()` and `prev()` are defined
for the given argument types.
Used for getting a tuple of the neighbouring kinks from some imaginary time.
"""
adjacent_kinks(ck::Any, ??::Any) = ck[prev(ck, ??)], ck[next(ck, ??)]


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
function adjacent_kinks_affecting_orbs(ck, os, ??)
    l = prev_affecting(ck, os, ??)

    if l == 0
        return (nothing,nothing)
    else
        return ck[l], ck[next_affecting(ck, os, ??)]
    end
end
adjacent_kinks_affecting_orbs(c::Configuration, os, ??) = adjacent_kinks_affecting_orbs(c.kinks, os, ??)


"""
     ??_prev_affecting(::Any, ::Any, ::Any)

Return the ImgTime of the closest kink to the left of ??
that affects one of the orbitals in os.
"""
function ??_prev_affecting(ck, os, ??)  # TODO: remove
    first(ck[prev_affecting(ck, os, ??)])
end


"""
     ??_next_affecting(::Any, ::Any, ::Any)

Return the ImgTime of the closest kink to the right of ??
that affects one of the orbitals in os.
"""
function ??_next_affecting(ck, os, ??)  # TODO: remove
    first(ck[next_affecting(ck, os, ??)])
end

"""
     ??_borders(::Kinks{T}, orbs, ::ImgTime) where {T <: Orbital}

Return a tuple of
the ImgTime of the closest kink to the right and
the ImgTime of the closest kink to the left of ??
that affect one of the orbitals in orbs.
If no orbital in os is affected by and kink from the collection in the first argument,
return a tuple of the interval bounds (ImgTime(0), ImgTime(1))."""
??_borders(ck::Kinks{T}, orbs, ??::ImgTime) where {T <: Orbital} = ImgTime(adjacent_kinks_affecting_orbs(ck, orbs, ??))  # TODO: remove
??_borders(c::Configuration{T}, orbs, ??::ImgTime) where {T <: Orbital} = ??_borders(c.kinks, orbs, ??)


"""
    isunaffected(Kinks, orb)

Return if an orbital is not affected by any kink.
"""
isunaffected(kinks, orb) = all(kink -> orb ??? last(kink), kinks)

Base.in(i, k::T2) = i == k.i || i == k.j
Base.in(i, k::T4) = i == k.i || i == k.j || i == k.k || i == k.l

"""
    in_open_interval(??, ??_first, ??_last)

Check if `??` lies in the open `ImgTime`-interval `(??_first,`??_last)`, assuming that `??_first` is the left border and `??_last` the right one.
"""
in_open_interval(??, ??_first, ??_last) = (??_first != ?? != ??_last != ??_first) && ((??_first < ??_last) == ((?? < ??_first) != (?? < ??_last)))

"""
    isunaffected_in_interval(kinks, orb, ??_first::ImgTime, ??_last::ImgTime)

Return if an orbital is not affected by any of the kinks from `ck` in the open interval `(??_first,??_last)`.
"""
isunaffected_in_interval(kinks, orb, ??_first::ImgTime, ??_last::ImgTime) = all(kink -> ((orb ??? last(kink)) || !in_open_interval(first(kink), ??_first, ??_last) ), kinks)


"""
    time_ordered_orbs(ck::Kinks{T}) where {T <: Orbital}

Get a list of orbitals that affect each kink in the time-ordering of the kinks and in the conventional ordering i, j, k, l.
"""
function time_ordered_orbs(ck::Kinks{T}) where {T <: Orbital}
    if isempty(ck)
        return Array{T,1}()
    else
        return collect( Iterators.flatten( orbs(k) for k in excitations(ck) ) )
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
    right_type_1_chain_length(ck, ??, count = 0)

Returns the length of the chain of type-1-entaglements starting with the Kink at ?? counting to the right.
"""
function right_type_1_chain_length(ck::Kinks, ??, counted_??s = [])
    next_kink = ck[next_affecting(ck, orbs(ck[??]), ??)]
    if is_type_1(ck[??], last(next_kink)) & !in(??, counted_??s)
        push!(counted_??s,??)
        return right_type_1_chain_length(ck, first(next_kink), counted_??s)
    else
        return length(counted_??s)
    end
end

"""
    left_type_1_chain_length(ck::Kinks, ??, counted_??s = [])

Returns the length of the chain of type-1-entaglements starting with the Kink at ?? counting to the left.
"""
function left_type_1_chain_length(ck::Kinks, ??, counted_??s = [])
    prev_kink = ck[prev(kinks_affecting_orbs(ck, orbs(ck[??])), ??)]
    if is_type_1(last(prev_kink), ck[??]) & !in(??, counted_??s)
        push!(counted_??s,??)
        return left_type_1_chain_length(ck, first(prev_kink), counted_??s)
    else
        return length(counted_??s)
    end
end

"""
    longest_type_1_chain_length(ck::Kinks)

Returns the longest chain of type-1-entaglements in ck.
"""
function longest_type_1_chain_length(ck::Kinks)
    longest_length = 0
    for (??, kink) in ck
        longest_length = max(right_type_1_chain_length(ck, ??),
                                left_type_1_chain_length(ck, ??), longest_length)
    end
    return longest_length
end

"""
    right_type_1_count(ck::Kinks)

Returns the number of kinks that are type 1 entagled with their right neighbour.
"""
function right_type_1_count(ck::Kinks)
    count = 0
    for (??, kink) in ck
        if is_type_1(kink, last(ck[next_affecting(ck, orbs(kink),??)]))
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
shuffle_annihilators(??::T4) = T4(shuffle_annihilators(??.i, ??.j, ??.k, ??.l)...)

"""
    shuffle_creators(::Any, ::Any, ::Any, ::Any)

Return a tuple where the first two arguments are randomly shuffled with equal probability 0.5
in comparison to the ordering of the input.
"""
shuffle_creators(i, j, k, l) = random_shuffle(i, j)..., k, l
shuffle_creators(??::T4) = T4(shuffle_creators(??.i, ??.j, ??.k, ??.l)...)

"""
    shuffle_indices(::Any, ::Any, ::Any, ::Any)

Return a tuple where both the first two arguments and the last two arguments
are randomly shuffled with equal probability 0.5 respectively
in comparison to the ordering of the input.
"""
shuffle_indices(i, j, k, l) = random_shuffle(i, j)..., random_shuffle(k, l)...
shuffle_indices(??::T4) = T4(shuffle_indices(??.i, ??.j, ??.k, ??.l)...)


"""
    kinks_from_periodic_interval(::Kinks, ??1, ??2)

return kinks with ?? ??? (??1,??2) if ??1 < ??2 and ?? ??? (??2,1) ??? (0,??1) if ??1 > ??2
"""
function kinks_from_periodic_interval(ck::Kinks, ??1, ??2)
    filter(x-> in_open_interval(first(x), ??1, ??2), ck)
end

"""
    times_from_periodic_interval(::Kinks, ::ImgTime, ::ImgTime)

return a list of all times of kinks with ?? ??? (??1,??2) if ??1 < ??2 or ?? ??? (??2,1) ??? (0,??1) if ??1 > ??2
in the periodic ordering suggested by the relation of the first time-argument ??1 to the second time-argument ??2
"""
function times_from_periodic_interval(ck::Kinks, ??1::ImgTime, ??2::ImgTime)
    times(filter(x-> in_open_interval(first(x), ??1, ??2), ck))
end
