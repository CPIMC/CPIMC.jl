using DataStructures
import Base: union, union!

abstract type Basis end

struct T2{T}
  " start orbital "
  i :: T

  " end orbital "
  j :: T
end

struct T4{T}
  " start orbitals "
  i :: T
  j :: T

  " end orbitals "
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

" multi-particle trajectory using single particle states with type T "
mutable struct Configuration{T}
  " set of orbitals occupied at tau=0 "
  occupations :: Set{T}

  " excitations, using tau as an index "
  kinks :: SortedDict{Float64, Kink{T}}

  Configuration(s::Set{T}) where {T} = new{T}(s,SortedDict{Float64,Kink}(Base.Forward))
end

" outer constructor method for empty Configurations{T} "
Configuration{T}() where T = Configuration(Set{T}())

" apply a T4 kink to a set of basis states "
function kink(o::Set{T}, κ::T4{T}) where T
  union(setdiff(o, Set([κ.i,κ.j])), Set([κ.k, κ.l]))
end

" apply a T4 kink to a set of basis states for a pair of a time and a T4-kink "
function kink(o::Set{T}, κ::Pair{Float64,T4{T}}) where T
  union(setdiff(o, Set([κ[2].i,κ[2].j])), Set([κ[2].k, κ[2].l]))
end

" return the occupied orbitals after applying all kinks to initial occupation "
function occupation(o::Set{T}, kinks::SortedDict{Float64,Kink{T}}) :: Set{T} where {T}
  reduce(kink, c.kinks; init=c.occupations)
end

" return the occupied orbitals to the right of τ "
function occupation(c::Configuration{T}, τ::Float64) :: Set{T} where T
  occupation(c.occupations, filter!(x -> x[1] <= τ, c.kinks))
end

function diff!(c::Configuration, n::Nothing)
  nothing
end

function diff!(c::Configuration{T}, o::T) where {T <: Basis}
  delete!(c.occupations, o)
end

function diff!(c::Configuration{T}, k::Kink{T}) where {T}
  delete!(c.kinks, k)
end

function diff!(c::Configuration{T}, ck::SortedDict{Float64,Kink{T}}) where {T}
  # TODO: find in-place operation for SortedDict, e.g. iterate delete!(token)
  throw(Exception("no in-place operation defined."))
end

#TODO:
# function diff(c1::Configuration{T}, c2::Configuration{T}) :: Configuration{T} where {T}
#   body
# end

function diff!(c1::Configuration{T}, c2::Configuration{T}) where {T}
  if isempty(c1.kinks)
    setdiff!(c1.occupations, c2.occupations)
  else
    if isempty(c2.kinks)
      setdiff!(c1.occupations, c2.occupations)
    else
      setdiff!(c1.occupations, c2.occupations)
      setdiff!(c1.kinks, c2.kinks)# TODO: find in-place operations for SortedDict
    end
  end
end

function union!(c::Configuration, n::Nothing)
  nothing
end

function union!(c::Configuration{T}, o::T) where {T <: Basis}
  push!(c.occupations, o)
end

function union!(c::Configuration{Basis}, co::Set{Basis})
  union!(c.occupations, co)
end

function union!(c::Configuration{T}, τ::Float64, k::Kink{T}) where {T}
  insert!(c.kinks, τ, k)
end

function union!(c::Configuration{T}, ck::SortedDict{Float64,Kink{T}}) where {T}
  # TODO: find in-place operation for SortedDict, e.g. iterate delete!(token)
  throw(Exception("no in-place operation defined."))
end

# TODO:
# function union(c1::Configuration{T}, c2::Configuration{T}) :: Configuration{T} where {T}
#   Configuration(union(c1.occupations,c2.occupations), merge(c1.kinks, c2.kinks))
# end

function union!(c1::Configuration{T}, c2::Configuration{T}) where {T}
  union!(c1.occupations, c2.occupations)
  merge!(c1.kinks, c2.kinks)
end

" get single particle basis type "
basis(c::Configuration{T}) where T = T
