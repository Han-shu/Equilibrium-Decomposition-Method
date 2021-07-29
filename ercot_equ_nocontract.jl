"""
Run the equilibrium algorithm without contract

Simplified from equilbrium model by setting vol_gen=0 and vol_cons=0

"""

using JuMP, Gurobi, JSON, Roots

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
    MARG_COST # C_{g}^{EN}
    LOAD_SHIFT # D_s
    α_GEN # α_g
    α_CONS # α_c
    β_GEN # β_g
    β_CONS # β_c
    γ # γ   Not used here
    SCEN_PROB # p_{frs}
    INV_COST # C_g^{INV}
end

mutable struct evaluation
    CAPACITY
    UNSERVED_DEM
    CONS_SURPLUS
    NET_REVENUE_UNIT
    IRR
    INC_IRR
end

mutable struct provisional_parameters
    ENERGY_PRICE # λ_{rst}
    OPER_PROFIT # π_{grst}
    CAPACITY # x_g
    DEM_SERVED_FIX # d_trs^{fix}
    DEM_SERVED_FLEX # d_trs^{flex}
end

function init_evaluation(param::Dict)
    CAPACITY = Dict(g => 0.0 for g in param["GENS"])
    UNSERVED_DEM = Dict((r,s) => 0.0 for r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
    CONS_SURPLUS = Dict((r,s)=>0.0 for r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
    NET_REVENUE_UNIT = Dict((g,r,s)=>0.0 for g in param["GENS"], r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
    IRR = Dict(g=>0.0 for g in param["GENS"])
    INC_IRR = 0.0
    return evaluation(CAPACITY, UNSERVED_DEM, CONS_SURPLUS, NET_REVENUE_UNIT, IRR, INC_IRR)
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
    AVAILABILITY = Dict((g,r,t) => param["availability"][g][string(r)][t] for g in param["GENS"], t in 1:param["N_SEGMENTS"], r in 1:param["N_PROFILES"])
    MARG_COST = Dict(g => param["fuel_price"][g] for g in param["GENS"])
    LOAD_SHIFT = param["adder"]
    α_GEN = Dict(g => 0.9 for g in param["GENS"])
    α_CONS = 0.9 #0.05
    β_GEN = Dict(g => 0.2 for g in param["GENS"])
    β_CONS = 0.2
    γ = 25.0
    SCEN_PROB = Dict((r, s) => param["p_profile"][r] * param["p_adder"][s] for r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
    INV_COST = param["investment"]

    return parameters(VOLL, DURATION, LOAD_FIX, LOAD_FLEX, AVAILABILITY, MARG_COST, LOAD_SHIFT,
    α_GEN, α_CONS, β_GEN, β_CONS, γ, SCEN_PROB, INV_COST)
end

function init_provisional_parameters(param::Dict)
    ENERGY_PRICE = Dict()
    OPER_PROFIT = Dict()
    CAPACITY = Dict("CCGT" => 18.078, "w_CCGT" =>  58.0875)
    # CAPACITY = Dict("CCGT" => 15.067, "w_CCGT" =>  61.4113)
    # CAPACITY = Dict("CCGT" => 5.27163, "w_CCGT" =>  74.2632)
    DEM_SERVED_FIX = Dict()
    DEM_SERVED_FLEX = Dict()

    for r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"]
        for t in 1:param["N_SEGMENTS"]
            ENERGY_PRICE[t,r,s] = 0.0
            DEM_SERVED_FIX[t,r,s] = 0.0
            DEM_SERVED_FLEX[t,r,s] = 0.0
            for g in param["GENS"]
                OPER_PROFIT[g,t,r,s] = 0.0
            end
        end
    end

    return provisional_parameters(ENERGY_PRICE, OPER_PROFIT, CAPACITY, DEM_SERVED_FIX, DEM_SERVED_FLEX)
end

function make_dispatch_model(set::sets, param::parameters, prov_param::provisional_parameters)
    dispatch_model = Model(Gurobi.Optimizer)

    @variable(dispatch_model, 0 <= dem_served_fix[set.TIME_BLOCKS, set.PROFILES, set.SHIFTERS])
    @variable(dispatch_model, 0 <= dem_served_flex[set.TIME_BLOCKS, set.PROFILES, set.SHIFTERS] <= param.LOAD_FLEX)
    @variable(dispatch_model, 0 <= prod[set.GENS, set.TIME_BLOCKS, set.PROFILES, set.SHIFTERS])

    @objective(dispatch_model, Max,
    (sum(sum(param.DURATION[t]*param.VOLL * (dem_served_fix[t,r,s] + dem_served_flex[t,r,s] - 0.5*(dem_served_flex[t,r,s]^2)/param.LOAD_FLEX) for t in set.TIME_BLOCKS)
    - sum(param.DURATION[t]*param.MARG_COST[g]*prod[g,t,r,s] for t in set.TIME_BLOCKS, g in set.GENS) for r in set.PROFILES, s in set.SHIFTERS)   )
    )

    @constraint(dispatch_model, energy_balance[t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS],
    (dem_served_fix[t,r,s] + dem_served_flex[t,r,s] + param.LOAD_SHIFT[s] == sum(prod[g,t,r,s] for g in set.GENS))
    )

    @constraint(dispatch_model, max_prod[g in set.GENS, t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS],
    (prod[g,t,r,s] <= param.AVAILABILITY[g,r,t] * prov_param.CAPACITY[g])
    )

    @constraint(dispatch_model, max_dem_served_fix[t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS],
    (dem_served_fix[t,r,s] <= param.LOAD_FIX[t]))

    return dispatch_model
end

function make_cons_model(set::sets, param::parameters, prov_param::provisional_parameters)
    cons_model = Model(Gurobi.Optimizer)

    @variable(cons_model, 0 <= scen_surplus_aux[set.PROFILES, set.SHIFTERS])
    @variable(cons_model, var_cons)

    @expression(cons_model, scen_surplus[r in set.PROFILES, s in set.SHIFTERS],
    (sum(param.DURATION[t] * param.VOLL * (prov_param.DEM_SERVED_FIX[t,r,s] + prov_param.DEM_SERVED_FLEX[t,r,s] - 0.5 * (prov_param.DEM_SERVED_FLEX[t,r,s]^2) / param.LOAD_FLEX) for t in set.TIME_BLOCKS)
    - sum(param.DURATION[t] * prov_param.ENERGY_PRICE[t,r,s] * (prov_param.DEM_SERVED_FIX[t,r,s] + prov_param.DEM_SERVED_FLEX[t,r,s] + param.LOAD_SHIFT[s]) for t in set.TIME_BLOCKS))
    )

    @objective(cons_model, Max,
    ((1 - param.β_CONS) * (var_cons - (1 / param.α_CONS) * (sum(param.SCEN_PROB[r,s] * scen_surplus_aux[r,s] for r in set.PROFILES, s in set.SHIFTERS)))
    + param.β_CONS * (sum(param.SCEN_PROB[r,s] * scen_surplus[r,s] for r in set.PROFILES, s in set.SHIFTERS)))
    ) # Mays AMPL model multiplies last term with sum of installed capacities.

    @constraint(cons_model, risk_set_cons[r in set.PROFILES, s in set.SHIFTERS],
    (var_cons - scen_surplus[r,s] <= scen_surplus_aux[r,s])
    )

    return cons_model
end

function make_gen_model(g, set::sets, param::parameters, prov_param::provisional_parameters)
    gen_model = Model(Gurobi.Optimizer)

    @variable(gen_model, 0 <= scen_profit_aux[set.PROFILES, set.SHIFTERS])
    @variable(gen_model, var_gen)

    @expression(gen_model, scen_profit[r in set.PROFILES, s in set.SHIFTERS],
    (-param.INV_COST[g] * prov_param.CAPACITY[g]
    + sum(prov_param.OPER_PROFIT[g,t,r,s] * prov_param.CAPACITY[g] for t in set.TIME_BLOCKS))
    ) # Mays AMPL model does not multiply by capacity and uses only marginal values. Multiplies OPER_PROFIT by duration here instead of in algorithm.
    # Q: Should OPER_PROFIT be multiplied with production instead of CAPACITY?

    @objective(gen_model, Max,
    ((1 - param.β_GEN[g]) * (var_gen - (1 / param.α_GEN[g]) * (sum(param.SCEN_PROB[r,s] * scen_profit_aux[r,s] for r in set.PROFILES, s in set.SHIFTERS)))
    + param.β_GEN[g] * (sum(param.SCEN_PROB[r,s] * scen_profit[r,s] for r in set.PROFILES, s in set.SHIFTERS)))
    ) # Mays AMPL model has β_CONS in second term and sums over all generators in first and second term.

    @constraint(gen_model, risk_set_gen[r in set.PROFILES, s in set.SHIFTERS],
    (var_gen - scen_profit[r,s] <= scen_profit_aux[r,s])
    )

    return gen_model
end


function init_input(input_path)
    param = JSON.parsefile(input_path)
    sets = init_sets(param)
    parameters = init_parameters(param)
    provisional_parameters = init_provisional_parameters(param)
    eval = init_evaluation(param)
    return sets, parameters, provisional_parameters, eval
end

function solve_dispatch(set, param, prov_param)
    dispatch_model = make_dispatch_model(set, param, prov_param)
    set_silent(dispatch_model)
    optimize!(dispatch_model)
    return dispatch_model
end

function solve_gen(set, param, prov_param)
    gen_model = Dict()
    for g in set.GENS
        gen_model[g] = make_gen_model(g, set, param, prov_param)
        set_silent(gen_model[g])
        optimize!(gen_model[g])
    end
    return gen_model
end

function solve_cons(set, param, prov_param)
    cons_model = make_cons_model(set, param, prov_param)
    set_silent(cons_model)
    optimize!(cons_model)
    return cons_model
end


function calc_gen_objective(set, param, gen_model)
    gen_objective = Dict()
    for g in set.GENS
        gen_objective[g] = ((1 - param.β_GEN[g]) * (value(gen_model[g][:var_gen]) - (1 / param.α_GEN[g]) * (sum(param.SCEN_PROB[r,s] * value(gen_model[g][:scen_profit_aux][r,s]) for r in set.PROFILES, s in set.SHIFTERS)))
        + param.β_GEN[g] * (sum(param.SCEN_PROB[r,s] * value(gen_model[g][:scen_profit][r,s]) for r in set.PROFILES, s in set.SHIFTERS)))
    end
    return gen_objective
end

function solve_equilibrium(set, param, provisional_parameters)
    finish = false
    avg_profit = nothing
    step = 1e-2
    iter = 1
    prov_param = provisional_parameters
    gen_objective = Dict(g => 0.0 for g in set.GENS)

    while !finish && iter < 50_000
        println("Iteration: ", iter)

        prov_param.CAPACITY = Dict(g => max(0, prov_param.CAPACITY[g] + step*gen_objective[g] / param.INV_COST[g]) for g in set.GENS)

        dispatch_model = solve_dispatch(set, param, prov_param)
        print(iter, "\t", prov_param.CAPACITY, "\n")
        # print(iter, "\t", shadow_price(dispatch_model[:energy_balance][1,1,1])/param.DURATION[1], "\t", shadow_price(dispatch_model[:max_prod]["CCGT",94,2,10]), "\n")
        print(iter, "\t", "ED Obj", "\t", objective_value(dispatch_model), "\n")

        for t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS
            prov_param.DEM_SERVED_FIX[t,r,s] = value(dispatch_model[:dem_served_fix][t,r,s])
            prov_param.DEM_SERVED_FLEX[t,r,s] = value(dispatch_model[:dem_served_flex][t,r,s])
            prov_param.ENERGY_PRICE[t,r,s] = shadow_price(dispatch_model[:energy_balance][t,r,s])/param.DURATION[t]
            for g in set.GENS
                prov_param.OPER_PROFIT[g,t,r,s] = shadow_price(dispatch_model[:max_prod][g,t,r,s])*param.AVAILABILITY[g,r,t]
            end
        end
        # print(iter, "\t", "OPER_PROFIT: ", prov_param.OPER_PROFIT[("CCGT", 94,2,10)],"\n")
        cons_model = nothing
        gen_model = nothing

        gen_model = solve_gen(set, param, prov_param)
        cons_model= solve_cons(set, param, prov_param)
        print("GEN Var value: ", " \t", value(gen_model["CCGT"][:var_gen]), "\t", value(gen_model["w_CCGT"][:var_gen]), "\n")
        print("CONS Var value: ", "\t", value(cons_model[:var_cons]), "\n")
        print("GEN Scen surplus: ", value(gen_model["CCGT"][:scen_profit][1,1]), "\n")
        print("CONS Scen surplus: ", value(cons_model[:scen_surplus][1,1]), "\n")
        print("GEN scen_profit_aux: ", "\t", value(gen_model["CCGT"][:scen_profit_aux][1,1]), "\n")
        print("CONS scen_surplus_aux: ", "\t", value(cons_model[:scen_surplus_aux][1,1]), "\n")
        gen_objective = calc_gen_objective(set, param, gen_model)
        print("GEN objective: ", " \t", gen_objective["CCGT"], " \t", gen_objective["w_CCGT"], "\n")
        print("GEN prod: ", value(dispatch_model[:prod]["CCGT",1,1,1]), "\n")
        print("GEN dem_served_flex: ", value(dispatch_model[:dem_served_flex][1,1,1]), "\n")
        max_risk_imbalance = maximum(map(x->abs(x), values(gen_objective)))
        println("Max risk imbalance: ", max_risk_imbalance)
        open("/Users/hanshu/Box/ERCOT/data/data7/nocon_max_risk_imbalance.json", "a") do f
            println(f, iter, " \t", max_risk_imbalance, " \t", gen_objective["CCGT"], " \t", gen_objective["w_CCGT"], 4)
        end

        if max_risk_imbalance < 1e-1
            finish = true
            avg_profit = Dict(g => sum(param.SCEN_PROB[r,s] * value(gen_model[g][:scen_profit_aux][r,s]) for r in set.PROFILES, s in set.SHIFTERS) for g in set.GENS)
        end

        iter += 1
    end

    return prov_param, avg_profit
end

function solve_evaluation(set, param, prov_param, evalu)
    gen_model = solve_gen(set, param, prov_param)
    cons_model = solve_cons(set, param, prov_param)
    dispatch_model = solve_dispatch(set, param, prov_param)
    for r in set.PROFILES, s in set.SHIFTERS
        evalu.UNSERVED_DEM[r,s] = (sum(param.LOAD_FIX[t] - prov_param.DEM_SERVED_FIX[t,r,s] for t in set.TIME_BLOCKS))
        evalu.CONS_SURPLUS[r,s] = value(cons_model[:scen_surplus][r,s])
        for g in set.GENS
            evalu.NET_REVENUE_UNIT[g,r,s] = sum(prov_param.OPER_PROFIT[g,t,r,s]  for t in set.TIME_BLOCKS)
        end
    end

    for g in set.GENS
        evalu.CAPACITY[g] = prov_param.CAPACITY[g]
        f(x) = sum(param.SCEN_PROB[r,s]*evalu.NET_REVENUE_UNIT[g,r,s] for r in set.PROFILES, s in set.SHIFTERS)*sum(1/((1+x)^n) for n in 1:20) - param.INV_COST[g]*sum(1/((1.04)^n) for n in 1:20)
        evalu.IRR[g] = find_zero(f, (0,5), Bisection())
    end
    h(x) = sum(param.SCEN_PROB[r,s]*(evalu.NET_REVENUE_UNIT["w_CCGT",r,s]-evalu.NET_REVENUE_UNIT["CCGT",r,s]) for r in set.PROFILES, s in set.SHIFTERS)*sum(1/((1+x)^n) for n in 1:20) - (param.INV_COST["w_CCGT"]-param.INV_COST["CCGT"])*sum(1/((1.04)^n) for n in 1:20)
    evalu.INC_IRR = find_zero(h, (-1,1), Bisection())
    return evalu
end


INPUT_FPATH = "../ercot_data.json"

set, param, prov_param, evalu = init_input(INPUT_FPATH)

prov_param, avg_profit = solve_equilibrium(set, param, prov_param)

evalu = solve_evaluation(set, param, prov_param, evalu)


OUTPUT_FPATH = "../out_no_contract.json"
open(OUTPUT_FPATH, "w") do f
    JSON.print(f, evalu,4)
end
