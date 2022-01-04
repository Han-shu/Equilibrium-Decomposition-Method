using JSON, JuMP, Gurobi, Roots

struct sets
    GENS # g
    TIME_BLOCKS # t
    PROFILES # r
    SHIFTERS # s
end

struct parameters
    VOLL # B
    DURATION # L
    LOAD_FIX # D_t^{fix}
    LOAD_FLEX # D_t^{res}
    AVAILABILITY # A_{grt}
    FUEL_PRICE # C_{g}^{EN}
    LOAD_SHIFT # D_s
    α
    β
    SCEN_PROB # p_{frs}
    INV_COST # C_g^{INV}
end

mutable struct evaluation
    CAPACITY
    UNSERVED_DEM
    CONS_SURPLUS
    EXPECTED_UNSERVED_DEM
    CONS_total_SURPLUS
    model_SURPLUS
    Objective
end

function init_evaluation(param::Dict)
    CAPACITY = Dict(g => 0.0 for g in param["GENS"])
    UNSERVED_DEM = Dict((r,s) => 0.0 for r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
    CONS_SURPLUS = Dict((r,s)=>0.0 for r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
    EXPECTED_UNSERVED_DEM = 0.0
    CONS_total_SURPLUS = 0.0
    model_SURPLUS = 0.0
    Objective = 0.0

    return evaluation(CAPACITY, UNSERVED_DEM, CONS_SURPLUS, EXPECTED_UNSERVED_DEM, CONS_total_SURPLUS, model_SURPLUS, Objective)
end

function init_sets(param::Dict)
    GENS = param["GENS"]
    TIME_BLOCKS = 1:param["N_SEGMENTS"]
    PROFILES = 1:param["N_PROFILES"]
    SHIFTERS = 1:param["N_SHIFTERS"]
    return sets(GENS, TIME_BLOCKS, PROFILES, SHIFTERS)
end

function init_parameters(param::Dict)
    VOLL = param["VOLL"]
    DURATION = param["duration"]
    LOAD_FIX = Dict(t => (param["nominal_demand"][t]-param["FLEX_DEMAND_MAX"]) for t in 1:param["N_SEGMENTS"])
    LOAD_FLEX = param["FLEX_DEMAND_MAX"]
    AVAILABILITY = Dict((g, r, t) => param["availability"][g][string(r)][t] for g in param["GENS"], t in 1:param["N_SEGMENTS"], r in 1:param["N_PROFILES"])
    FUEL_PRICE = Dict(g => param["fuel_price"][g] for g in param["GENS"])
    LOAD_SHIFT = param["adder"]
    α = 0.9 #0.05
    β = 0.2 #0.4
    SCEN_PROB = Dict((r, s) => param["p_profile"][r] * param["p_adder"][s] for r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
    INV_COST = param["investment"]

    return parameters(VOLL, DURATION, LOAD_FIX, LOAD_FLEX,
    AVAILABILITY, FUEL_PRICE, LOAD_SHIFT, α, β, SCEN_PROB, INV_COST)
end


function init_input(input_path)
    param = JSON.parsefile(input_path)
    sets = init_sets(param)
    parameters = init_parameters(param)
    eval = init_evaluation(param)
    return sets, parameters, eval
end

function make_dispatch_model(set::sets, param::parameters)
    model = Model(Gurobi.Optimizer)

    @variable(model, capacity[set.GENS] >= 0)
    @variable(model, 0 <= dem_served_fix[set.TIME_BLOCKS, set.PROFILES, set.SHIFTERS])
    @variable(model, 0 <= dem_served_flex[set.TIME_BLOCKS, set.PROFILES, set.SHIFTERS])
    @variable(model, 0 <= prod[set.GENS, set.TIME_BLOCKS, set.PROFILES, set.SHIFTERS])
    @variable(model, 0 <= surplus_aux[set.PROFILES, set.SHIFTERS])
    @variable(model, var)

    @expression(model, surplus[r in set.PROFILES, s in set.SHIFTERS],
        (-sum(param.INV_COST[g] * capacity[g] for g in set.GENS)
        + sum(param.VOLL * param.DURATION[t] * (dem_served_fix[t,r,s] + dem_served_flex[t,r,s] - (0.5*(dem_served_flex[t,r,s]^2)/param.LOAD_FLEX)) for t in set.TIME_BLOCKS)
        - sum(param.FUEL_PRICE[g] * param.DURATION[t] * prod[g,t,r,s] for g in set.GENS, t in set.TIME_BLOCKS))
        )

    @objective(model, Max,
        ((1 - param.β) * (var - (1 / param.α)*(sum(param.SCEN_PROB[r,s]*surplus_aux[r,s] for r in set.PROFILES, s in set.SHIFTERS)))
        + param.β * (sum(param.SCEN_PROB[r,s] * surplus[r,s] for r in set.PROFILES, s in set.SHIFTERS))  )
        )

    @constraint(model, risk_set[r in set.PROFILES, s in set.SHIFTERS], (var - surplus[r,s] <= surplus_aux[r,s])
        )

    @constraint(model, energy_balance[t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS],
        (dem_served_fix[t,r,s] + dem_served_flex[t,r,s] + param.LOAD_SHIFT[s] == sum(prod[g,t,r,s] for g in set.GENS))
        )

    @constraint(model, max_prod[g in set.GENS, t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS],
        (prod[g,t,r,s] <= param.AVAILABILITY[g,r,t] * capacity[g]))

    @constraint(model, max_dem_served_fix[t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS],
        (dem_served_fix[t,r,s] <= param.LOAD_FIX[t]))

    @constraint(model, max_dem_served_flex[t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS],
        (dem_served_flex[t,r,s] <= param.LOAD_FLEX))

    optimize!(model)
    return model
end

function solve_dispatch(set, param)
    dispatch_model = make_dispatch_model(set, param)
    set_silent(dispatch_model)
    optimize!(dispatch_model)
    return dispatch_model
end


function solve_evaluation(set, param, evalu)
    model = solve_dispatch(set, param)
    evalu.CAPACITY = Dict(g => value(model[:capacity][g]) for g in set.GENS)
    for r in set.PROFILES, s in set.SHIFTERS
        evalu.UNSERVED_DEM[r,s] = (sum(param.LOAD_FIX[t] -  value(model[:dem_served_fix][t,r,s]) for t in set.TIME_BLOCKS))
    end

    evalu.EXPECTED_UNSERVED_DEM = sum(param.SCEN_PROB[r,s] * evalu.UNSERVED_DEM[r,s]  for r in set.PROFILES, s in set.SHIFTERS)
    evalu.CONS_total_SURPLUS = sum(param.SCEN_PROB[r,s] * evalu.CONS_SURPLUS[r,s]  for r in set.PROFILES, s in set.SHIFTERS)
    evalu.model_SURPLUS = sum(param.SCEN_PROB[r,s] * value(model[:surplus][r,s]) for r in set.PROFILES, s in set.SHIFTERS)
    evalu.Objective = objective_value(model)

    return evalu
end


INPUT_FPATH = "../ercot_data.json"
set, param, evalu = init_input(INPUT_FPATH)

evalu = solve_evaluation(set, param, evalu)

OUTPUT_FPATH = "../Gurobi_complete.json"
open(OUTPUT_FPATH, "w") do f
    JSON.print(f, evalu,4)
end
