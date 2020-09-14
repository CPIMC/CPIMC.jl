## Calculate a mapping from orbital indices to orbital quantum numbers
## e.g. kvectors
## The units are so that the possible values of each k-vector component are integer numbers (i.e. divide by 2pi/boxarea compared to a.u.)

using DataFrames
import LinearAlgebra: dot

# wrapper for the sign
# struct Sign
#     sign :: Bool
# end
# TODO : Avoid fields with abstract type
# mutable struct Sign{T<:Bool}
    # sign :: T
# end

function get_index(q::Array{Int16,1}) :: Int64
    b = convert(Array{Int64,1}, q)
    b[1] << 32 | b[2] << 16 | b[3]
end

function get_vector(index) :: Array{Int16,1}
    convert(Array{Int16,1}, [ (index >>> 32) & 0xFFFF, (index >>> 16) & 0xFFFF, index & 0xFFFF ])
end

function get_energy(index)
    vec = get_vector(index)
    dot(vec,vec)
end

function get_Nb(emax::Int) :: Int
    kk = Int(floor(sqrt(emax)))
    nb :: Int = 0
    for x in -kk:kk
        for y in -kk:kk
            for z in -kk:kk
                if x*x + y*y + z*z <= emax
                    nb = nb + 1
                end
            end
        end
    end
    nb
end

function get_Emax(Nb::Int) :: Int
    ek :: Int = 0
    while !(get_Nb(ek) == Nb)
        ek += 1
        @assert ek < Nb "Basis size not found: ek = $ek >= $Nb = Nb. Maybe this basis size does not correspond to a full shell?"
    end
    ek
end

function get_basis(e::Ensemble)
    b :: Set{Int64} = Set{Int64}()
    nb = 0

    emax = get_Emax(e.cutoff)

    # estimate upper bound for maximum k-Component
    # kc :: Int = Int(ceil(sqrt(e.cutoff^(1.748626417827964) + 13*e.cutoff)))# cutoff radius (TODO: better estimate)
    kk = Int(floor(sqrt(emax)))
    for x in -kk:kk
        for y in -kk:kk
            for z in -kk:kk
                if x*x + y*y + z*z <= emax
                    indx = get_index(Array{Int16,1}([x,y,z]))
                    @assert !in(indx, b) "index $(indx) of array [$x,$y,$z] already in Set."
                    push!(b, indx)
                    nb += 1
                end
            end
        end
    end
    print("SET :\n")
    print(b)
    @assert length(b) == e.cutoff "Basis size $(length(b)) not equal cutoff=$(e.cutoff). Hitted $nb orbitals."
    b
end

function emptyOrbs(e::Ensemble, c::Configuration) :: Set{Int}
    setdiff!(c.occupations, get_basis(e))
end

### Estimators
function Ekin(e::Ensemble, c::Configuration)
    sum(get_energy(n,orblist) for n in c.occupations)
end

function occVec(c::Configuration)
    Dict(Symbol(ci) => 1 for ci in c.occupations)
    c.occupations
end

function particleNumber(c::Configuration)
  return length(c.occupations)
end

# function emptyOrbs(e::Ensemble, c::Configuration)
  # return filter(x -> !(x in c.occupations), 1:e.cutoff)
# end



### Units
function get_beta_internal(theta, N)
  return ((2*pi)^2)/(((6*(pi^2)*N)^(2/3))*theta)
end
