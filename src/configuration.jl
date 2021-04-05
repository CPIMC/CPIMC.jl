export Configuration, Orbital, Kink, T2, T4, ImgTime

@doc raw"""
    T2{T}

Parametric type representing a 1-particle excitation by specifying
a transition from one state (annihilator) to another (creator).
In the occupation number representation, this reads

    `a^{\dagger}_i a_j`

with a creator orbital `i` and an annihilator orbital `j`.
The single-particle basis is represented by the type parameter T.
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

    `a^{\dagger}_i a^{\dagger}_j a_k a_l`

with creator orbitals `i`,`j` and an annihilator orbitals `k`, `l`.
The single-particle basis is represented by the type parameter T.
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
    const Kink{T}

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
const Kink{T} = Union{T2{T},T4{T}}

" outer constructor method to construct a T2 kink, inferring the type parameter from the arguments "
Kink(i,j) = T2(i,j)
" outer constructor method to construct a T4 kink, inferring the type parameter from the arguments "
Kink(i,j,k,l) = T4(i,j,k,l)
""" outer constructor method to extract a kink from a pair where the second element is a kink.
    This is useful for automatic conversion when looping over SortedDict{S,Kink{T}} """
Kink(p::Pair{S,T} where {T<:Kink} where {S}) = p[2]# first substitute S, then T

"""
    orbs(::T2)

return a set of all orbitals which are affected by a T2 kink
"""
orbs(x::T2) = Set([k.i, k.j])

"""
    orbs(x::T4)

return a set of all orbitals which are affected by a T4 kink """
orbs(x::T4) = Set([x.i, x.j, x.k, x.l])

"""
    creators(::T4)

return a set of the two creators which are affected by a T4 kink
"""
creators(x::T4) = Set([x.i, x.j])

"""
    annihilators(x::T4)

return a set of the two annihilators which are affected by a T4 kink
"""
annihilators(x::T4) = Set([x.k, x.l])

""" type alias for imaginary time
    FixedPointNumbers are used since these are stable for == and are thus stable as keys in Dict """
const ImgTime = Fixed{Int64,60}

" outer constructor method for a Tuple{ImgTime,ImgTime} from Tuple{Nothing,Nothing} which returns the bounds of the ImgTime interval as Tuple (ImgTime(0), ImgTime(1)) "
ImgTime(t::Tuple{Nothing,Nothing}) = (ImgTime(0), ImgTime(1))

"""
    ImgTime(::Pair{ImgTime,<:Kink})

outer constructor method for an ImgTime from Pair{ImgTime,<:Kink}
this returns the first argument of the Pair
useful for automatic conversion
"""
ImgTime(p::Pair{ImgTime,<:Kink}) = first(p)


" outer constructor method for a Tuple{ImgTime,ImgTime} from Tuple{Pair{ImgTime,<:Kink},Pair{ImgTime,<:Kink}} which returns the ImgTimes of the Pairs as Tuple{ImgTime,ImgTime} "
ImgTime(t::Tuple{Pair{ImgTime,<:Kink},Pair{ImgTime,<:Kink}}) = (first(first(t)), first(last(t)))

@doc raw"""
    Δ(::ImgTime,::ImgTime)

return the length of the periodic interval between the two ::ImgTimes

the periodic distance of two imaginary times:

    `Δ(τ1,τ2) = τ2 - τ1` if `τ2 >= τ1` and
    `Δ(τ1,τ2) = 1 - (τ1 - τ2) = 1 + τ2 - τ1` else
"""
Δ(τ1::ImgTime,τ2::ImgTime) = τ1 > τ2 ? ImgTime(1) + τ2 - τ1 : τ2 - τ1

" multi-particle trajectory using single particle states with type T "
mutable struct Configuration{T}
  " set of orbitals occupied at τ=0 "
  occupations :: Set{T}
  " excitations, using τ as an index "
  kinks :: SortedDict{ImgTime, Kink{T}}
end

" outer constructor method for a configuration with occupations given by o and kinks given by k. k can be anything from which a SortedDict can be constructed from. "
Configuration(o::Set{T}, k) where {T} = Configuration(o, SortedDict{ImgTime,Kink{T}}(Base.Forward, k) )
" outer constructor method for a configuration with occupations given by o and no kinks. "
Configuration(o::Set{T}) where {T} = Configuration(o, SortedDict{ImgTime,Kink{T}}(Base.Forward))
" outer constructor method for a configuriation with no occupations and kinks given by k. k can be anything from which a SortedDict can be constructed from. "
Configuration(k::SortedDict{ImgTime,<:Kink{T}}) where {T} = Configuration(Set{T}(), k)
" outer constructor method for a configuration with no occupations and kinks as given by varargs p... of Pair{ImgTime,<:Kink}, which are passed to the SortedDict constructor "
Configuration(p::Pair{ImgTime,<:Kink{T}}...) where {T} = Configuration(SortedDict{ImgTime,Kink{T}}(Base.Forward, p...))
" outer constructor method for a configuration with occupations given by o and kinks as given by varargs p... of Pair{ImgTime,<:Kink}, which are passed to the SortedDict constructor "
Configuration(o::Set{T}, p::Pair{ImgTime,<:Kink{T}}...) where {T} = Configuration(o, SortedDict{ImgTime,Kink{T}}(Base.Forward, p...))

" outer constructor method for empty Configurations{T} "
Configuration{T}() where T = Configuration(Set{T}())

" get single particle basis type "
basis(c::Configuration{T}) where T = T

" abstract type for single-particle basis states, implementation is required for each model "
abstract type Orbital end

@doc raw"""
    excite(::Set{T}, ::T4{T}) where T

'Apply a T4 kink to a set of basis states', i.e.
return a set of basis states where the states specified by the creators of the kink are added
and the states specified by the annihilators of the kink are dropped from the given set of basis states.
This has the physical meaning of a two-particle scattering event where
two (quasi-)particles in a many-body state change the single-particle states they occupy.
In occupation number representation this reads

    `a^{\dagger}_i a^{\dagger}_j a_k a_l |\{n\}\rangle`

for creator orbitals `i` and `j` and annihilator orbitals `k` and `l`.
This function assumes fermionic particle statistics,

    `a^{\dagger}_i a^{\dagger}_i = a_i a_i = 0` (Pauli principle)

i.e. the target (creator) states must not be occupied and
the initial (annihilator) states must be occupied in the given set of states.
"""
function excite(o::Set{T}, κ::T4{T}) where T
  @assert ( in(κ.k, o) & in(κ.l, o) ) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the annihilators $(κ.k), $(κ.l) is not occupied."
  @assert ( !in(κ.i, o) & !in(κ.j, o) ) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the creators $(κ.i), $(κ.j) is already occupied. (Pauli principle)"
  union(setdiff(o, Set([κ.k,κ.l])), Set([κ.i, κ.j]))
end

""" Apply a T4 kink to a set of basis states for a pair of a time and a T4-kink.
    This is useful for iteration of a SortedDict{ImgTime, T4{T}}."""
excite(o::Set{T}, κ::Pair{ImgTime,T4{T}}) where T = excite(o, last(κ))


"""
    excite!(::Set{T}, i::T, j::T, k::T, l::T) where {T}

Apply an excitation in-place to a set of basis states,
that is given by creating orbitals i, j
and annihilating the orbitals k, l
"""
function excite!(o::Set{T}, i::T, j::T, k::T, l::T) where {T}
    @assert (in(k, o) & in(l, o)) "Kink ($(i),$(j),$(k),$(l)) cannot be applied: one or two of the annihilators $(k), $(l) is not occupied. (Pauli-Principle)"
    @assert (!in(i, o) & !in(j, o)) "Kink ($(i),$(j),$(k),$(l)) cannot be applied: one or two of the creators $(i), $(j) is already occupied. (Pauli-Principle)"
    delete!(o, k)
    delete!(o, l)
    push!(o, i)
    push!(o, j)
end

" Apply a T4 kink in-place to a set of basis states. "
excite!(o::Set{T}, κ::T4{T}) where T = excite!(o, κ.i, κ.j, κ.k, κ.l)


" Return the occupied orbitals after applying all kinks to initial occupation. "
function occupations(o::Set{T}, kinks::SortedDict{ImgTime,Kink{T}}) :: Set{T} where {T}
  foldl(excite, kinks; init=o)
end

" Return the occupied orbitals to the right of τ ."
function occupations(c::Configuration, τ::ImgTime)
  occupations(c.occupations, filter(x -> first(x) <= τ, c.kinks))
end

"""
    next(::SortedDict{ImgTime,<:Kink}, ::ImgTime)

Return first kink after to given ::ImgTime. If there is no kink with ::ImgTime larger than τ, return first kink.
"""
function next(ck::SortedDict{ImgTime,<:Kink}, τ::ImgTime)
    toc = searchsortedafter(ck, τ)
    if toc == pastendsemitoken(ck)
        first(ck)
    else
        deref((ck, toc))
    end
end

"""
    prev(::SortedDict{ImgTime,<:Kink}, ::ImgTime)

Return first kink before given ::ImgTime. If there is no kink with ::ImgTime smaller than τ, return last kink.
"""
function prev(ck::SortedDict{ImgTime,<:Kink}, τ::ImgTime)
    toc = regress((ck, searchsortedfirst(ck, τ)))
    if toc == beforestartsemitoken(ck)
        last(ck)
    else
        deref((ck, toc))
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
    if isempty(κs)
        ImgTime(0.0)
    else
        ImgTime(prev(κs, τ))
    end
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
    if isempty(κs)
        ImgTime(1.0)
    else
        ImgTime(next(κs, τ))
    end
end

"""
     τ_borders(::SortedDict{ImgTime,<:Kink}, ::Set{T}, ::ImgTime) where {T <: Orbital}

Return a tuple of
the ImgTime of the closest kink to the right and
the ImgTime of the closest kink to the left of τ
that affect one of the orbitals in os.
If no orbital in os is affected by and kink from the collection in the first argument,
return a tuple of the interval bounds (ImgTime(0), ImgTime(1))."""
τ_borders(ck::SortedDict{ImgTime,<:Kink{T}}, os::Set{T}, τ::ImgTime) where {T <: Orbital} = ImgTime(adjacent_kinks_affecting_orbs(ck, os, τ))

τ_borders(c::Configuration{T}, os::Set{T}, τ::ImgTime) where {T <: Orbital} = τ_borders(c.kinks, os, τ)


" Return if an orbital is not affected by any kink. "
function isunaffected(ck::SortedDict{ImgTime,<:Kink{T}}, orbital::T) where {T<:Orbital}
    if isempty(ck)
        return true
    else
        return !in(orbital, union( orbs.( values(ck))... ) )
    end
end

" Return if an orbital is not affected by any of the kinks from ck in the open interval (τ_first,τ_last). "
function isunaffected_in_interval(ck::SortedDict{ImgTime,<:Kink{T}}, orbital::T, τ_first::ImgTime, τ_last::ImgTime) :: Bool where {T<:Orbital}
    @assert τ_first != τ_last
    if τ_first < τ_last
        for (τ_kink,kink) in ck
            if (τ_kink <= τ_first) | (τ_kink >= τ_last)
                continue
            elseif in(orbital, orbs(kink))
                return false
            end
        end
    else
        for (τ_kink,kink) in ck
            if ((τ_kink <= τ_first) & (τ_kink >= τ_last))
                continue
            elseif in(orbital, orbs(kink))
                return false
            end
        end
    end
    return true
end

## Method definitions for function drop
# This function is mostly used for calculating the changes proposed by an update
# in order to calculate proposal and/or acceptance probabilities.

" return the configuration c without dropping anything. "
drop(c::Configuration, n::Nothing) = c

" return a set with the orbital o dropped from oc. "
drop(oc::Set{T}, o::T) where {T <: Orbital} = setdiff(oc, Set([o]))
" return a configuration with the orbital o dropped from c.occupations. "
drop(c::Configuration{T}, o::T) where {T <: Orbital} = Configuration(drop(c.occupations, Set([o])),c.kinks)

" return a set with the orbitals in oc2 dropped from oc1. "
drop(oc1::Set{T}, oc2::Set{T}) where {T <: Orbital} = setdiff(oc1, oc2)
" return a configuration with the orbitals in oc dropped from c.occupations. "
drop(c::Configuration{T}, oc::Set{T}) where {T <: Orbital} = Configuration(drop(c.occupations, oc), c.kinks)

" return a SortedDict{ImgTime,<:Kink{<:Orbital}} with the pairs in ck2 dropped from ck1. "
drop(ck1::SortedDict{ImgTime,<:Kink{T}}, ck2::SortedDict{ImgTime,<:Kink{T}}) where {T <: Orbital} = setdiff(ck1, ck2)
" return a configuration with the pairs in ck dropped from c.kinks. "
drop(c::Configuration{T}, ck::SortedDict{ImgTime,<:Kink{T}}) where {T <: Orbital} = Configuration(c.occupations, drop(c.kinks, ck))

" return a SortedDict{ImgTime,<:Kink{<:Orbital}} with the pairs ps... dropped from ck. "
drop(ck::SortedDict{ImgTime,<:Kink{T}}, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital} = setdiff(ck, SortedDict(ps...))
" return a configuration with the pairs ps... dropped from c.kinks. "
drop(c::Configuration{T}, ps::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital} = Configuration(c.occupations, drop(c.kinks, ps...))

" return a configuration with occupations and kinks in c2 dropped from c1. "
drop(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital} = Configuration(drop(c1.occupations,c2.occupations), drop(c1.kinks,c2.kinks))


## Method definitions for function drop!
# This function is mostly for applying the changes determined in an MC Step Δ
# to the current configuration c in function apply_step!(c, Δ). cf. CPIMC.jl
# These methods change the first argument in-place.

" drop nothing "
function drop!(c::Configuration, n::Nothing)
    nothing
end

" drop a single orbital "
function drop!(c::Configuration{T}, o::T) where {T <: Orbital}
    delete!(c.occupations, o)
end

" drop all orbitals given in a set oc from c.occupations "
function drop!(c::Configuration{T}, oc::Set{T}) where {T <: Orbital}
    for o in oc
        drop!(c, o)
    end
end

" drop a single kink at time τ from c.kinks "
function drop!(c::Configuration, τ::ImgTime)
    delete!(c.kinks, τ)
end

" drop all kinks given in a SortedDict from c.kinks "
function drop!(c::Configuration{T}, ck::SortedDict{ImgTime,<:Kink{T}}) where {T <: Orbital}
    for (τ, k) in ck
        drop!(c, τ)
    end
end

" drop all kinks given by Varargs{Pair{ImgTime,<:Kink{T}}} from c.kinks "
function drop!(c::Configuration{T}, pk::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital}
    for (τ, k) in pk
        drop!(c, τ)
    end
end

" drop occupations and kinks given by second argument c2 from first argument c1 "
function drop!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}
    drop!(c1, c2.occupations)
    drop!(c1, c2.kinks)
end


## Method definitions for function add.
# This function is mostly used for calculating the changes proposed by an update
# in order to calculate proposal and/or acceptance probabilities.

" return the configuration c without adding anything "
add(c::Configuration, n::Nothing) = c

" return a set with the orbital o added to oc "
add(oc::Set{T}, o::T) where {T <: Orbital} = union(oc, Set([o]))
" return a configuration with the orbital o added to c.occupations "
add(c::Configuration{T}, o::T) where {T <: Orbital} = Configuration(add(c.occupations, o), c.kinks)

" return a set with the orbitals in oc2 added to oc1 "
add(oc1::Set{T}, oc2::Set{T}) where {T <: Orbital} = union(oc1, oc2)
" return a configuration with the orbitals in oc added to c.occupations "
add(c::Configuration{T}, oc::Set{T}) where {T <: Orbital} = Configuration(union(c.occupations, oc), c.kinks)

" return a SortedDict{ImgTime,<:Kink{<:Orbital}} with the pairs ps... added to ck "
add(ck::SortedDict{ImgTime,<:Kink{T}}, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital} = merge(ck, SortedDict(ps...))
" return a configuration with the pairs ps... added to c.kinks "
add(c::Configuration{T}, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital} = Configuration(c.occupations, add(c.kinks, ps...))
" return a configuration with the pairs in ck added to c.kinks "
add(c::Configuration{T}, ck::SortedDict{ImgTime,<:Kink{T}}) where {T <: Orbital} = add(c, ck...)

" return a configuration with occupations and kinks in c2 dropped from c1 "
add(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital} = Configuration(add(c1.occupations, c2.occupations), merge(c1.kinks, c2.kinks))


## Method definitions for function add!.
# This function is mostly for applying the changes determined in an MC Step Δ
# to the current configuration c in function apply_step!(c, Δ). cf. CPIMC.jl
# These methods change the first argument in-place.

" add nothing "
function add!(c::Configuration, n::Nothing)
    nothing
end

" add a single orbital "
function add!(c::Configuration{T}, o::T) where {T <: Orbital}
    push!(c.occupations, o)
end

" add all orbitals given in a set oc to c.occupations "
function add!(c::Configuration{T}, oc::Set{T}) where {T <: Orbital}
    union!(c.occupations, oc)
end

" add a single kink at time τ to c.kinks "
function add!(c::Configuration{T}, τ::ImgTime, k::Kink{T}) where {T <: Orbital}
    insert!(c.kinks, τ, k)
end

" add all kinks given by the pairs ps... to ck::SortedDict{ImgTime,<:Kink{<:Orbital}} "
add!(ck::SortedDict{ImgTime,<:Kink{T}}, ps::Pair{ImgTime, <:Kink{T}}...) where {T <: Orbital} = merge!(ck, SortedDict(ps...))
" add all kinks given by the pairs ps... to c.kinks "
add!(c::Configuration{T}, ps::Pair{ImgTime,<:Kink{T}}...) where {T <: Orbital} = merge!(c.kinks, SortedDict(ps))

" add all kinks given in a SortedDict to c.kinks "
function add!(c::Configuration{T}, ck::SortedDict{ImgTime,<:Kink{T}}) where {T <: Orbital}
    merge!(c.kinks, ck)
end

" add occupations and kinks given by second argument c2 from first argument c1 "
function add!(c1::Configuration{T}, c2::Configuration{T}) where {T <: Orbital}
    add!(c1, c2.occupations)
    add!(c1, c2.kinks)
end


""" get a list of orbitals that are affected by a kink in the conventional ordering i, j, k, l
    this corresponds to the matrix element w_ijkl in the same order """
time_ordered_orbs(x::T4) = [x.i, x.j, x.k, x.l]


""" get a list of orbitals that affect each kink in the time-ordering of the kinks and in the conventional ordering i, j, k, l """
function time_ordered_orbs(ck::SortedDict{ImgTime,<:Kink{T}}) where T time_ordered_orbs
    if isempty(ck)
        return Array{T,1}()
    else
        return vcat([time_ordered_orbs(k) for k in values(ck)]...)
    end
end

""" returns 1 or -1 depending on the order of all ladder operators as given by a list of orbitals """
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


""" returns 1 or -1 depending on the order of all ladder operators as given by the orbitals that affect each kink in the time-ordering of the kinks and in the conventional ordering i, j, k, l,
used in the sign estimator."""
ladder_operator_order_factor(ck::SortedDict{ImgTime,<:Kink}) = ladder_operator_order_factor(time_ordered_orbs(ck))




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
function right_type_1_chain_length(ck::SortedDict{ImgTime,<:Kink}, τ, counted_τs = [])
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
    left_type_1_chain_length(ck, τ, count = 0)

Returns the length of the chain of type-1-entaglements starting with the Kink at τ counting to the left.
"""
function left_type_1_chain_length(ck::SortedDict{ImgTime,<:Kink}, τ, counted_τs = [])
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
    longest_type_1_chain_length(ck)

Returns the longest chain of type-1-entaglements in ck.
"""
function longest_type_1_chain_length(ck::SortedDict{ImgTime,<:Kink})
    longest_length = 0
    for (τ, kink) in ck
        longest_length = max(right_type_1_chain_length(ck, τ),
                                left_type_1_chain_length(ck, τ), longest_length)
    end
    return longest_length
end

"""
    longest_type_1_chain_length(ck)

Returns the longest chain of type-1-entaglements in ck.
"""
function right_type_1_count(ck::SortedDict{ImgTime,<:Kink})
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
    kinks_from_periodic_interval(::SortedDict{ImgTime,<:Kink}, τ1, τ2)

return kinks with τ ∈ (τ1,τ2) if τ1 < τ2 and τ ∈ (τ2,1) ∪ (0,τ1) if τ1 > τ2
"""
function kinks_from_periodic_interval(ck::SortedDict{ImgTime,<:Kink}, τ1, τ2)
    if τ1 > τ2# interval is periodically continued
        filter(x -> ( τ1 < first(x) ) | ( first(x) < τ2 ), ck)
    else# this also catches τ1 == τ2
        filter(x -> τ1 < first(x) < τ2, ck)
    end
end

"""
    times_from_periodic_interval(::SortedDict{ImgTime,<:Kink}, ::ImgTime, ::ImgTime)

return a list of all times of kinks with τ ∈ (τ1,τ2) if τ1 < τ2 or τ ∈ (τ2,1) ∪ (0,τ1) if τ1 > τ2
in the periodic ordering suggested by the relation of the first time-argument τ1 to the second time-argument τ2
"""
function times_from_periodic_interval(ck::SortedDict{ImgTime,<:Kink}, τ1::ImgTime, τ2::ImgTime)
    if τ1 > τ2# interval is periodically continued
        vcat(
            collect(keys( filter(x -> τ1 < first(x), ck) )),
            collect(keys( filter(x -> first(x) < τ2, ck) ))
            )
    else# this also catches τ1 == τ2
        collect(keys( filter(x -> τ1 < first(x) < τ2, ck) ))
    end
end
