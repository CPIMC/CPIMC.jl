using DataStructures
import Base: union, union!

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

Kink(i,j) = T2(i,j)
Kink(i,j,k,l) = T4(i,j,k,l)


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

function diff!(c::Configuration, n::Nothing)
  nothing
end

function diff!(c::Configuration{T}, o::T) where {T <: Orbital}
  delete!(c.occupations, o)
end

function diff!(c::Configuration{T}, k::Kink{T}) where {T}
  delete!(c.kinks, k)
end

function diff!(c::Configuration{T}, ck::SortedDict{Float64,Kink{T}}) where {T}
  # TODO: find in-place operation for SortedDict, e.g. iterate delete!(token)
  throw(Exception("no in-place operation defined."))
end

# function diff(c1::Configuration{T}, c2::Configuration{T}) :: Configuration{T} where {T}
#   if isempty(c1.kinks)
#     return Configuration(setdiff(c1.occupations,c2.occupations))
#   else
#     if isempty(c2.kinks)
#       return Configuration(setdiff(c1.occupations,c2.occupations), c1.kinks)
#     else
#       return Configuration(setdiff(c1.occupations, c2.occupations), setdiff(c1.kinks, c2.kinks))
#     end
#   end
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

function union!(c::Configuration{T}, o::T) where {T <: Orbital}
  push!(c.occupations, o)
end

function union!(c::Configuration{Orbital}, co::Set{Orbital})
  union!(c.occupations, co)
end

function union!(c::Configuration{T}, k::Kink{T}) where {T}
  insert!(c.kinks, k)
end

function union!(c::Configuration{T}, ck::SortedDict{Float64,Kink{T}}) where {T}
  # TODO: find in-place operation for SortedDict, e.g. iterate delete!(token)
  throw(Exception("no in-place operation defined."))
end

# function union(c1::Configuration{T}, c2::Configuration{T}) :: Configuration{T} where {T}
#   Configuration(union(c1.occupations,c2.occupations), merge(c1.kinks, c2.kinks))
# end

function union!(c1::Configuration{T}, c2::Configuration{T}) where {T}
  union!(c1.occupations, c2.occupations)
  merge!(c1.kinks, c2.kinks)
end

" get single particle basis type "
basis(c::Configuration{T}) where T = T
