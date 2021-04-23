export shuffle_indices, add_remove_kink_chain

"""
    shuffle_indices(m::Model, e::Ensemble, c::Configuration)
Performs an Update that takes a Kink and randomly shuffles its annihlators with each other and its creators.
The Sets of creators and annihilators of the kinks stay the same the acceptance propability is therefore always 1.
"""
function shuffle_indices(m::Model, e::Ensemble, c::Configuration)
    if isempty(c.kinks)
        return 1.0, Step()
    end
    kink = rand(c.kinks)
    if rand() > 0.5
        # shuffle creators
        Δ = Step(nothing, (kink,), nothing, (first(kink) => T4(last(kink).j,last(kink).i,last(kink).k,last(kink).l),))
    else
        # shuffle annihilators
        Δ = Step(nothing, (kink,), nothing, (first(kink) => T4(last(kink).i,last(kink).j,last(kink).l,last(kink).k),))
    end
    apply_step!(c,Δ)
    return 1.0, Δ
end

""" perform a number of subsequent updates
    first perform a number of updates which add kinks
    second perform the same number of updates which remove kinks
    THIS DOES NOT WORK CORRECTLY RIGHT NOW"""
function add_remove_kink_chain(m::Model, e::Ensemble, c::Configuration)
    if isempty(c.kinks)
        return 1.0, Step()
    end
    add_single_kinks = [add_type_C, add_type_D, add_type_E]
    remove_single_kinks = [remove_type_C, remove_type_D, remove_type_E]
    chain_length = rand(1:9)
    acc_prob = 1.0
    step_list = Array{Step,1}()# TODO: it may be more efficient to pre-allocate this output and only set it here
    for i in  1:chain_length
        dv, Δ = rand(add_single_kinks)(m,e,apply_step(c,step_list))
        acc_prob *= dv
        if iszero(acc_prob)
            return 0.0, Step()
        end
        push!(step_list, Δ)
    end
    for i in  1:chain_length
        dv, Δ = rand(remove_single_kinks)(m,e,apply_step(c,step_list))
        acc_prob *= dv
        if iszero(acc_prob)
            return 0.0, Step()
        end
        push!(step_list, Δ)
    end
    return acc_prob, step_list
end

"""
    isuseful(c::Configuration, up::Function)
Function with the name isuseful checkt wether it is senseble to propose a spezific Update-Type.
If no explizit function is defined for the type spezific type of Update its called on, this function will be called, which just returns true.
"""
function isuseful(c::Configuration, up::Function)
    return true
end


function isuseful(c::Configuration, up::typeof(shuffle_indices))
    if isempty(c.kinks)
        return false
    else
        return true
    end
end
