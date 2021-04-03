using DataStructures
using FixedPointNumbers

struct T2{T}
  " creator "
  i :: T

  " annihilator "
  j :: T
end

struct T4{T}
  " creator "
  i :: T
  j :: T

  " annihilator "
  k :: T
  l :: T
end

const Kink{T} = Union{T2{T},T4{T}}

" outer constructor method to construct a T2 kink, inferring the type parameter from the arguments "
Kink(i,j) = T2(i,j)
" outer constructor method to construct a T4 kink, inferring the type parameter from the arguments "
Kink(i,j,k,l) = T4(i,j,k,l)
""" outer constructor method to extract a kink from a pair where the second element is a kink.
    This is useful for automatic conversion when looping over SortedDict{S,Kink{T}} """
Kink(p::Pair{S,T} where {T<:Kink} where {S}) = p[2]# first substitute S, then T

" return a set of all orbitals which are affected by a T2 kink "
orbs(x::T2) = Set([k.i, k.j])

" return a set of all orbitals which are affected by a T4 kink "
orbs(x::T4) = Set([x.i, x.j, x.k, x.l])

" return a set of the two creators which are affected by a T4 kink "
creators(x::T4) = Set([x.i, x.j])

" return a set of the two annihilators which are affected by a T4 kink "
annihilators(x::T4) = Set([x.k, x.l])

""" type alias for imaginary time
    FixedPointNumbers are used since these are stable for == and are thus stable as keys in Dict """
const ImgTime = Fixed{Int64,60}

" outer constructor method for a Tuple{ImgTime,ImgTime} from Tuple{Nothing,Nothing} which returns the bounds of the ImgTime interval as Tuple (ImgTime(0), ImgTime(1)) "
ImgTime(t::Tuple{Nothing,Nothing}) = (ImgTime(0), ImgTime(1))

" outer constructor method for a Tuple{ImgTime,ImgTime} from Tuple{Pair{ImgTime,<:Kink},Pair{ImgTime,<:Kink}} which returns the ImgTimes of the Pairs as Tuple{ImgTime,ImgTime} "
ImgTime(t::Tuple{Pair{ImgTime,<:Kink},Pair{ImgTime,<:Kink}}) = (first(first(t)), first(last(t)))

@doc raw"""
    Δ(::ImgTime,::ImgTime)

return the periodic difference of two imaginary times
    `Δ(τ1,τ2) = τ1 - τ2` if `τ1 >= τ2` and
    `Δ(τ1,τ2) = 1 - (τ2 - τ1) = 1 + τ1 - τ2` else
"""
Δ(τ1::ImgTime,τ2::ImgTime) = τ1 < τ2 ? ImgTime(1) + τ1 - τ2 : τ1 - τ2

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

" apply a T4 kink to a set of basis states "
function excite(o::Set{T}, κ::T4{T}) where T
  @assert ( in(κ.k, o) & in(κ.l, o) ) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the annihilators $(κ.k), $(κ.l) is not occupied. (Pauli-Principle)"
  @assert ( !in(κ.i, o) & !in(κ.j, o) ) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the creators $(κ.i), $(κ.j) is already occupied. (Pauli-Principle)"
  union(setdiff(o, Set([κ.k,κ.l])), Set([κ.i, κ.j]))
end

""" Apply a T4 kink to a set of basis states for a pair of a time and a T4-kink.
    This is useful for iteration of a SortedDict{ImgTime, T4{T}}."""
excite(o::Set{T}, κ::Pair{ImgTime,T4{T}}) where T = excite(o, last(κ))

" Apply a T4 kink in-place to a set of basis states. "
function excite!(o::Set{T}, κ::T4{T}) where T
  @assert (in(κ.k, o) & in(κ.l, o)) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the annihilators $(κ.k), $(κ.l) is not occupied. (Pauli-Principle)"
  @assert (!in(κ.i, o) & !in(κ.j, o)) "Kink ($(κ.i),$(κ.j),$(κ.k),$(κ.l)) cannot be applied: one or two of the creators $(κ.i), $(κ.j) is already occupied. (Pauli-Principle)"
  delete!(o, κ.k)
  delete!(o, κ.l)
  push!(o, κ.i)
  push!(o, κ.j)
end

" Return the occupied orbitals after applying all kinks to initial occupation. "
function occupations(o::Set{T}, kinks::SortedDict{ImgTime,Kink{T}}) :: Set{T} where {T}
  reduce(excite, kinks; init=o)
end

" Return the occupied orbitals to the right of τ ."
function occupations(c::Configuration, τ::ImgTime)
  occupations(c.occupations, filter(x -> first(x) <= τ, c.kinks))
end

""" return first kink after to τ
    if there is no kink with ImgTime larger than τ, return first kink """
function next(ck::SortedDict{ImgTime,<:Kink}, τ::ImgTime)
    toc = searchsortedafter(ck, τ)
    if toc == pastendsemitoken(ck)
        first(ck)
    else
        deref((ck, toc))
    end
end

""" return first kink before τ
    if there is no kink with ImgTime smaller than τ, return last kink """
function prev(ck::SortedDict{ImgTime,<:Kink}, τ::ImgTime)
    toc = regress((ck, searchsortedfirst(ck, τ)))
    if toc == beforestartsemitoken(ck)
        last(ck)
    else
        deref((ck, toc))
    end
end

"""     kinks_affecting_orbs(c::Configuration{T}, os::Set{T}) where {T <: Orbital}
return all kinks in c which affect one ore more of the orbitals in os
"""
kinks_affecting_orbs(c::Configuration{T}, os::Set{T}) where {T <: Orbital} = filter( p -> any( Set([last(p).i,last(p).j,last(p).k,last(p).l]) .∈ (os,) ), c.kinks)

kinks_affecting_orbs(ck::SortedDict{ImgTime,<:Kink}, os::Set{T}) where {T <: Orbital} = filter( p -> any( Set([last(p).i,last(p).j,last(p).k,last(p).l]) .∈ (os,) ), ck)
"""     adjacent_kinks(ck::SortedDict{ImgTime,<:Kink}, τ::ImgTime) where {T <: Orbital}
return a tuple of imaginary times the closest kink to the right and
the closest kink to the left of τ that affect one of the orbitals in os
if there are no such kinks in ck, return (nothing,nothing)
"""
adjacent_kinks(ck::SortedDict{ImgTime,<:Kink}, τ::ImgTime) where {T <: Orbital} = prev(ck, τ), next(ck, τ)

"""    adjacent_kinks_affecting_orbs(c::Configuration{T}, os::Set{T}, τ::ImgTime) where {T <: Orbital}

return a tuple of
the closest kink to the right and the closest kink to the left of τ
that affect one of the orbitals in os if there are no such kinks in c, return (nothing,nothing) """
function adjacent_kinks_affecting_orbs(c::Configuration{T}, os::Set{T}, τ::ImgTime) where {T <: Orbital}
    k = kinks_affecting_orbs(c, os)
    if isempty(k)
        return (nothing,nothing)
    else
        return adjacent_kinks(k, τ)
    end
end

"""     τ_borders(c::Configuration{T}, os::Set{T}, τ::ImgTime) where {T <: Orbital}
return a tuple of
the ImgTime of the closest kink to the right and
the ImgTime of the closest kink to the left of τ
that affect one of the orbitals in os """
τ_borders(c::Configuration{T}, os::Set{T}, τ::ImgTime) where {T <: Orbital} = ImgTime(adjacent_kinks_affecting_orbs(c, os, τ))

" return if an orbital is not affected by any kink "
function isunaffected(ck::SortedDict{ImgTime,<:Kink{T}}, orbital::T) where {T<:Orbital}
    if isempty(ck)
        return true
    else
        return !in(orbital, union( orbs.( values(ck))... ) )
    end
end

" return if an orbital is not affected by any of the kinks from ck in the open interval (τ_first,τ_last). "
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


#TODO change Name to orbs_time_ordered
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


@doc raw"""
    wminus(i,j,k,l)

return the difference of the two-particle matrix element with the same but the last two indices transposed
antisymmetric difference of the two-particle matrix elements:
    `w^-_{ijkl} = w_{ijkl} - w_{ijlk}`
This is called the antisymmetrized two-particle matrix element.
This is an abbreviation for these terms arising in the
Slater-Condon rules for the calculation of the many-body matrix elements
via the one- and two-particle matrix elements of the underlying single-particle basis.
"""
wminus(i,j,k,l) = w(i,j,k,l) - w(i,j,l,k)


@doc raw"""
    Woffdiag_element(::Ensemble, ::Orbital, ::Orbital, ::Orbital, ::Orbital)

Return the offdiagonal many body matrix element of the interaction operator

    `\frac{1}{4} ( w_{ijkl} - w_{ijlk} ) ( \pm \braket{\tilde{\onvv} | \onvv^{ij}_{kl}} )`

for an excitation given by creating orbitals i,j and annihilating orbitals k, l
as given by the Slater-Condon rules
"""
function Woffdiag_element(e::Ensemble, i::Orbital, j::Orbital, k::Orbital, l::Orbital)
    # We sample with the weight of antisymmetrized matrix element but we do not restrict
    # the order of indices of our possible kinks. We therefor need an extra factor 1/4 in the weight-function
    1/4 * e.λ * wminus(i,j,k,l)
end

Woffdiag_element(e::Ensemble, kink::T4) = Woffdiag_element(e, kink.i, kink.j, kink.k, kink.l)

"""
    ΔWoffdiag_element(e::Ensemble, itr, itr)

return the change in the offdiagonal many body matrix element of the interaction operator
given by adding the kinks in the first iterable and removing the kinks in the second iterable

both arguments `itr` are required to be iterables containing kinks
"""
function ΔWoffdiag_element(e::Ensemble, add_kinks, drop_kinks)
    if isempty(drop_kinks)
        prod(Woffdiag_element(e, k) for k in add_kinks)
    elseif isempty(add_kinks)
        1.0 / prod(Woffdiag_element(e, k) for k in drop_kinks)
    else
        prod(Woffdiag_element(e, k) for k in add_kinks) / prod(Woffdiag_element(e, k) for k in drop_kinks)
    end
end

ΔWoffdiag_element(e::Ensemble, add_kinks::SortedDict{ImgTime,<:Kink}, drop_kinks::SortedDict{ImgTime,<:Kink}) = ΔWoffdiag_element(e, values(add_kinks), values(drop_kinks))


""" change in the diagonal interaction matrix element due to a change in the occupation occ
    in an periodic interval (τ1,τ2) where no kinks occur
    the change in the occupation is assumed to consist in a creation of two orbitals i, j
    and in the annihlation of two orbitals k, l """
function ΔW_diag(i, j, k, l, occ)
    @assert (i ∉ occ) & (j ∉ occ) "Calculation of the change in the many-body diagonal interaction matrix element: This function assumes that the first two orbitals\n\t $(i)\n and\n\t $(j) given are the creator orbitals and thus that they are not occupied in the given occupation\n\t $(occ). "
    # contributions due to mean field interactions of the annihilated orbitals
    Δ = sum( w(ν,k,k,ν) + w(ν,l,l,ν) for ν in drop(occ, Set([k,l])) )# interactions of mean field with k and l
    Δ += w(k,l,l,k)# interaction between k and l
    # contributions due to mean field interactions of the created orbitals
    # the annihilator orbitals k, l are not in the new occupation
    Δ -= sum( w(ν,i,i,ν) + w(ν,j,j,ν) for ν in drop(occ, Set([k,l])) )# interactions of mean field with i and j
    Δ -= w(i,j,j,i)# interaction between i and j
    return Δ
end

""" change in the diagonal interaction matrix element due to a change in the occupation occ
    in an periodic interval (τ1,τ2) where no kinks occur
    the change in the occupation is assumed to consist in a creation of one orbitals i
    and in the annihlation of one orbitals j """
function ΔW_diag(i, j, occ)
    @assert (i ∉ occ) "Calculation of the change in the many-body diagonal interaction matrix element: This function assumes that the first two orbitals\n\t $(i)\n and\n\t $(j) given are the creator orbitals and thus that they are not occupied in the given occupation\n\t $(occ). "
    # contributions due to mean field interactions of the annihilated orbitals
    Δ = sum( w(ν,j,j,ν) for ν in drop(occ, j) )# interactions of mean field with k
    # contributions due to mean field interactions of the created orbitals
    Δ -= sum( w(ν,i,i,ν) for ν in drop(occ, j) )
end

"""
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
return a list of all times of kinks with τ ∈ (τ1,τ2) if τ1 < τ2 or τ ∈ (τ2,1) ∪ (0,τ1) if τ1 > τ2
in the periodic ordering suggested by the relation of the first time-argument τ1 to the second time-argument τ2
"""
function times_from_periodic_interval(ck::SortedDict{ImgTime,<:Kink}, τ1::ImgTime, τ2::ImgTime)
    if τ1 > τ2# interval is periodically continued
        vcat([τ1],
            collect(keys( filter(x -> τ1 < first(x), ck) )),
            collect(keys( filter(x -> first(x) < τ2, ck) )),
            [τ2])
    else# this also catches τ1 == τ2
        vcat([τ1], collect(keys( filter(x -> τ1 < first(x) < τ2, ck) )), [τ2])
    end
end

### convention: all Δelement represent only the difference in the matrix elements will be used as exp(-Δ) for the weight change

""" change in the kinetic many body matrix element due to creating i, j and annihilating k,l """
ΔT_element(i,j,k,l) = energy(i) + energy(j) - energy(k) - energy(l)



"""
    ΔWdiag_element(c::Configuration, e::Ensemble, i, j, k, l, τ1, τ2)

calculate the change in the diagonal interaction many-body matrix element
due to a change in the occupations given by creating two orbitals i and j
and annihilating two orbitals k, l in the interval (τ1, τ2)
This interval may be periodically extended over the bounds (0,1) if τ1 > τ2,
i.e. the change in the occupation is considered for (τ1,1] ∪ [0,τ2) in that case.
We do not need to evaluate the diagonal interaction
between all obritals in all time-intervalls, but it is sufficient to evaluate
the full diagonal interaction with the occupations at the start of the Intervall
and then consider only contrubations of orbitals that are changed by kinks in the intervall.
"""# TODO: ΔWdiag_element
function ΔWdiag_element(c::Configuration, e::Ensemble, i, j, k, l, τ1, τ2)# TODO: assuming that i, j are creators and k, l are annihilators. Use Step instead ?
    @assert τ1 != τ2 " The diagonal interaction matrix element changes when kinks are added at different times and thus the occupations between the kinks are altered. It has no meaning to calculate this matrix element (or to add kinks) at equal times τ1=$(τ1), τ2=$(τ2). "
    τs = times_from_periodic_interval(c.kinks, τ1, τ2)
    Ks = kinks_from_periodic_interval(c.kinks, τ1, τ2)
    #Calculate Wdiag with occupation at the start of the Intervall.
    ΔWdiag_τ1_occ = ΔW_diag(i, j, k, l, occupations(c,τ1)) * Δ(τ2,τ1)
    if isempty(Ks)
        return e.λ * ΔWdiag_τ1_occ
    else
        #Calculate contrubutions to Wdiag of the orbitals changed by kinks in the intervall.
        ΔWdiag_kinks = sum( (ΔW_diag(i, j, k, l, creators(Ks[t1])) - ΔW_diag(i, j, k, l, annihilators(Ks[t1]))) * Δ(τ2,t1) for t1 in τs[2:end-1])
        return e.λ * (ΔWdiag_τ1_occ + ΔWdiag_kinks)
    end
end

"""
    ΔWdiag_element(c::Configuration, e::Ensemble, i, j, τ1, τ2)

calculate the change in the diagonal interaction many-body matrix element
due to a change in the occupations given by creating an orbital i
and annihilating an orbital j in the interval (τ1, τ2)
This interval may be periodically extended over the bounds (0,1) if τ1 > τ2,
i.e. the change in the occupation is considered for (τ1,1] ∪ [0,τ2) in that case.
We do not need to evaluate the diagonal interaction
between all obritals in all time-intervalls, but it is sufficient to evaluate
the full diagonal interaction with the occupations at the start of the Intervall
and then consider only contrubations of orbitals that are changed by kinks in the intervall.
"""
function ΔWdiag_element(c::Configuration, e::Ensemble, i, j, τ1, τ2)# TODO: assuming that i is creator and j is annihilator. Use Step instead ?
    @assert τ1 != τ2 " The diagonal interaction matrix element changes when kinks are added at different times and thus the occupations between the kinks are altered. It has no meaning to calculate this matrix element (or to add kinks) at equal times τ1=$(τ1), τ2=$(τ2). "
    τs = times_from_periodic_interval(c.kinks, τ1, τ2)
    Ks = kinks_from_periodic_interval(c.kinks, τ1, τ2)
    #Calculate Wdiag with occupation at the start of the Intervall.
    ΔWdiag_τ1_occ = ΔW_diag(i, j, occupations(c,τ1)) * Δ(τ2,τ1)
    if isempty(Ks)
        return e.λ * ΔWdiag_τ1_occ
    else
        #Calculate contrubutions to Wdiag of the orbitals changed by kinks in the intervall.
        ΔWdiag_kinks = sum( (ΔW_diag(i, j, creators(Ks[t1])) - ΔW_diag(i, j, annihilators(Ks[t1]))) * Δ(τ2,t1) for t1 in τs[2:end-1])
        return e.λ * (ΔWdiag_τ1_occ + ΔWdiag_kinks)
    end
end




"""
    sign_offdiagonal_product(c::Configuration)

sign of the product of the offdiagonal elements
used for the calculation the sign of the weight function
"""
function sign_offdiagonal_product(c::Configuration)
    sign_ofd_prod = 1
    for ϰ in values(c.kinks)
        sign_ofd_prod *= sign(wminus(ϰ.i, ϰ.j, ϰ.k, ϰ.l))
    end
    return sign_ofd_prod
end
