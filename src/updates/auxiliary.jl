
function shuffle_indices(m::Model, e::Ensemble, c::Configuration)
    if isempty(c.kinks)
        return 1.0, Step()
    end
    kink = rand(c.kinks)
    if rand() > 0.5
        # shuffle creators
        Δ = Step(Configuration(kink), Configuration(first(kink) => T4(last(kink).j,last(kink).i,last(kink).k,last(kink).l)))
    else
        # shuffle annihilators
        Δ = Step(Configuration(kink), Configuration(first(kink) => T4(last(kink).i,last(kink).j,last(kink).l,last(kink).k)))
    end
    return 1.0, Δ
end

""" perform a number of subsequent updates
    first perform a number of updates which add kinks
    second perform the same number of updates which remove kinks """
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