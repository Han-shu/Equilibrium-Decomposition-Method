    """
    Module to run the equilibrium algorithm of ercot paper

    """

    using JuMP, Gurobi, JSON, Roots

    struct sets
        GENS # g
        TIME_BLOCKS # t
        PROFILES # r
        SHIFTERS # s
        CONTRACTS # c
    end

    struct parameters
        VOLL # B
        DURATION # L
        LOAD_FIX # D_t^{fix}
        LOAD_FLEX # F_t^{res}
        AVAILABILITY # A_{grt}
        MARG_COST # C_{g}^{EN}
        LOAD_SHIFT # D_s
        α_GEN # α_g
        α_CONS # α_c
        β_GEN # β_g
        β_CONS # β_c
        γ # γ
        SCEN_PROB # p_{rs}
        INV_COST # C_g^{INV}
    end

    mutable struct provisional_parameters
        ITERATION
        ENERGY_PRICE # λ_{rst}
        OPER_PROFIT # π_{grst}
        PAYOUT # η^k_{rs}
        CONTRACT_PRICE # ϕ^k
        CAPACITY # x_q
        DEM_SERVED_FIX # d_t^{fix}
        DEM_SERVED_FLEX # d_t^{flex}
        VOL_CONS # v_c^k
        VOL_GEN # v_g^k
        VOL_MIN_GEN # \underline{v}_g^k
        VOL_MIN_CONS # \underline{v}_c^k
        VOL_MAX_GEN # \overline{v}_g^k
        VOL_MAX_CONS # \overline{v}_c^k
        GEN_objective
    end

    mutable struct evaluation
        CAPACITY
        UNSERVED_DEM
        CONS_SURPLUS
        NET_REVENUE_UNIT1
        NET_REVENUE_UNIT2
        AVE_REVENUE1
        AVE_REVENUE2
        IRR
        INC_IRR
    end

    function init_evaluation(param::Dict)
        CAPACITY = Dict(g=>0.0 for g in param["GENS"])
        UNSERVED_DEM = Dict((r,s) => 0.0 for r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
        CONS_SURPLUS = Dict((r,s)=>0.0 for r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
        NET_REVENUE_UNIT1 = Dict((g,r,s)=>0.0 for g in param["GENS"], r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
        NET_REVENUE_UNIT2 = Dict((g,r,s)=>0.0 for g in param["GENS"], r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
        AVE_REVENUE1 = Dict(g => 0.0 for g in param["GENS"])
        AVE_REVENUE2 = Dict(g => 0.0 for g in param["GENS"])
        IRR = Dict(g=>0.0 for g in param["GENS"])
        INC_IRR = 0.0
        return evaluation(CAPACITY, UNSERVED_DEM, CONS_SURPLUS, NET_REVENUE_UNIT1, NET_REVENUE_UNIT2, AVE_REVENUE1, AVE_REVENUE2, IRR, INC_IRR)
    end

    function init_sets(param::Dict)
        GENS = param["GENS"]
        TIME_BLOCKS = 1:param["N_SEGMENTS"]
        PROFILES = 1:param["N_PROFILES"]
        SHIFTERS = 1:param["N_SHIFTERS"]
        CONTRACTS = param["CONTRACTS"]
        return sets(GENS, TIME_BLOCKS, PROFILES, SHIFTERS, CONTRACTS)
    end

    function init_parameters(param::Dict)
        VOLL = param["VOLL"]
        DURATION = param["duration"]
        LOAD_FIX = Dict(t => param["nominal_demand"][t]-param["FLEX_DEMAND_MAX"] for t in 1:param["N_SEGMENTS"])
        LOAD_FLEX = param["FLEX_DEMAND_MAX"]
        AVAILABILITY = Dict((g, r, t) => param["availability"][g][string(r)][t] for g in param["GENS"], t in 1:param["N_SEGMENTS"], r in 1:param["N_PROFILES"])
        MARG_COST = Dict(g => param["fuel_price"][g] for g in param["GENS"])
        LOAD_SHIFT = param["adder"]
        α_GEN = Dict(g => 0.9 for g in param["GENS"])
        α_CONS = 0.9
        β_GEN = Dict(g => 0.2 for g in param["GENS"])
        β_CONS = 0.2
        γ = 50.0
        SCEN_PROB = Dict((r, s) => param["p_profile"][r] * param["p_adder"][s] for r in 1:param["N_PROFILES"], s in 1:param["N_SHIFTERS"])
        INV_COST = param["investment"]

        return parameters(VOLL, DURATION, LOAD_FIX, LOAD_FLEX, AVAILABILITY, MARG_COST, LOAD_SHIFT,
        α_GEN, α_CONS, β_GEN, β_CONS, γ, SCEN_PROB, INV_COST)
    end

    function init_provisional_parameters(set::sets)
        ITERATION = 1
        ENERGY_PRICE = Dict()
        OPER_PROFIT = Dict()
        PAYOUT = Dict()
        CONTRACT_PRICE = Dict(c => 0.0 for c in set.CONTRACTS)

        CAPACITY = Dict("CCGT" => 77.0449869359021, "w_CCGT" => 0.9472080921646084)

        DEM_SERVED_FIX = Dict()
        DEM_SERVED_FLEX = Dict()
        VOL_CONS = Dict(c => 0.0 for c in set.CONTRACTS)
        VOL_GEN = Dict((g, c) => 0.0 for g in set.GENS, c in set.CONTRACTS)

        for r in set.PROFILES, s in set.SHIFTERS
            for t in set.TIME_BLOCKS
                ENERGY_PRICE[t,r,s] = 0.0
                DEM_SERVED_FIX[t,r,s] = 0.0
                DEM_SERVED_FLEX[t,r,s] = 0.0

                for g in set.GENS
                    OPER_PROFIT[g,t,r,s] = 0.0
                end
            end

            for c in set.CONTRACTS
                PAYOUT[c,r,s] = 0.0
            end
        end

        VOL_MIN_GEN = Dict((g, c) => 0.0 for g in set.GENS, c in set.CONTRACTS)
        VOL_MIN_CONS = Dict(c => 0.0 for c in set.CONTRACTS)
        VOL_MAX_GEN = Dict((g, c) => 0.0 for g in set.GENS, c in set.CONTRACTS)
        VOL_MAX_CONS = Dict(c => 0.0 for c in set.CONTRACTS)
        GEN_objective = Dict(g => 0.0 for g in set.GENS)

        return provisional_parameters(ITERATION, ENERGY_PRICE, OPER_PROFIT, PAYOUT, CONTRACT_PRICE, CAPACITY, DEM_SERVED_FIX,
        DEM_SERVED_FLEX, VOL_CONS, VOL_GEN, VOL_MIN_GEN, VOL_MIN_CONS, VOL_MAX_GEN, VOL_MAX_CONS, GEN_objective)
    end

    function make_dispatch_model(set::sets, param::parameters, prov_param::provisional_parameters)
        dispatch_model = Model(Gurobi.Optimizer)

        @variable(dispatch_model, 0 <= dem_served_fix[set.TIME_BLOCKS, set.PROFILES, set.SHIFTERS])
        @variable(dispatch_model, 0 <= dem_served_flex[set.TIME_BLOCKS, set.PROFILES, set.SHIFTERS] <= param.LOAD_FLEX)
        @variable(dispatch_model, 0 <= prod[set.GENS, set.TIME_BLOCKS, set.PROFILES, set.SHIFTERS])

        @objective(dispatch_model, Max,
        (sum(sum(param.DURATION[t]*param.VOLL*(dem_served_fix[t,r,s] + dem_served_flex[t,r,s] - 0.5 * (dem_served_flex[t,r,s]^2) / param.LOAD_FLEX) for t in set.TIME_BLOCKS)
        - sum(param.DURATION[t]*param.MARG_COST[g]*prod[g,t,r,s] for t in set.TIME_BLOCKS, g in set.GENS) for r in set.PROFILES, s in set.SHIFTERS))
        )

        @constraint(dispatch_model, energy_balance[t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS],
        (dem_served_fix[t,r,s] + dem_served_flex[t,r,s] + param.LOAD_SHIFT[s] == sum(prod[g,t,r,s] for g in set.GENS))
        )

        @constraint(dispatch_model, max_prod[g in set.GENS, t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS],
        (prod[g,t,r,s] <= param.AVAILABILITY[g,r,t] * prov_param.CAPACITY[g])
        )

        @constraint(dispatch_model, max_dem_served_fix[t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS],
        (dem_served_fix[t,r,s] <= param.LOAD_FIX[t] ))

        return dispatch_model
    end

    function make_cons_model(set::sets, param::parameters, prov_param::provisional_parameters)
        cons_model = Model(Gurobi.Optimizer)

        @variable(cons_model, prov_param.VOL_MIN_CONS[c] <= vol_cons[c in set.CONTRACTS]  <= prov_param.VOL_MAX_CONS[c])
        @variable(cons_model, 0 <= scen_surplus_aux[set.PROFILES, set.SHIFTERS])
        @variable(cons_model, var_cons)

        @expression(cons_model, scen_surplus[r in set.PROFILES, s in set.SHIFTERS],
        (-sum(vol_cons[c] * (prov_param.CONTRACT_PRICE[c] - prov_param.PAYOUT[c,r,s]) for c in set.CONTRACTS)
        + sum(param.DURATION[t] * param.VOLL * (prov_param.DEM_SERVED_FIX[t,r,s] + prov_param.DEM_SERVED_FLEX[t,r,s] - 0.5 * (prov_param.DEM_SERVED_FLEX[t,r,s]^2) / param.LOAD_FLEX) for t in set.TIME_BLOCKS)
        - sum(param.DURATION[t] * prov_param.ENERGY_PRICE[t,r,s] * (prov_param.DEM_SERVED_FIX[t,r,s] + prov_param.DEM_SERVED_FLEX[t,r,s] + param.LOAD_SHIFT[s] ) for t in set.TIME_BLOCKS))
        )

        @objective(cons_model, Max,
        ((1 - param.β_CONS) * (var_cons - (1 / param.α_CONS) * (sum(param.SCEN_PROB[r,s] * scen_surplus_aux[r,s] for r in set.PROFILES, s in set.SHIFTERS)))
        + param.β_CONS * (sum(param.SCEN_PROB[r,s] * scen_surplus[r,s] for r in set.PROFILES, s in set.SHIFTERS)))
         # - 0.5 * param.γ * (sum((vol_cons[c] + sum(prov_param.VOL_GEN[g,c] for g in set.GENS))^2 for c in set.CONTRACTS)))
        ) # Mays AMPL model multiplies last term with sum of installed capacities.

        @constraint(cons_model, risk_set_cons[r in set.PROFILES, s in set.SHIFTERS],
        (var_cons - scen_surplus[r,s] <= scen_surplus_aux[r,s])
        )

        return cons_model
    end


    function make_gen_model(g, set::sets, param::parameters, prov_param::provisional_parameters)
        gen_model = Model(Gurobi.Optimizer)

        @variable(gen_model, prov_param.VOL_MIN_GEN[g,c] <= vol_gen[c in set.CONTRACTS] <= prov_param.VOL_MAX_GEN[g,c])
        @variable(gen_model, 0 <= scen_profit_aux[set.PROFILES, set.SHIFTERS])
        @variable(gen_model, var_gen)

        @expression(gen_model, scen_profit[r in set.PROFILES, s in set.SHIFTERS],
        (-param.INV_COST[g] * prov_param.CAPACITY[g]
        - sum(vol_gen[c] * (prov_param.CONTRACT_PRICE[c] - prov_param.PAYOUT[c,r,s]) for c in set.CONTRACTS)
        + sum(prov_param.OPER_PROFIT[g,t,r,s] * prov_param.CAPACITY[g] for t in set.TIME_BLOCKS))
        )

        @objective(gen_model, Max,
        ((1 - param.β_GEN[g]) * (var_gen - (1 / param.α_GEN[g]) * (sum(param.SCEN_PROB[r,s] * scen_profit_aux[r,s] for r in set.PROFILES, s in set.SHIFTERS)))
        + param.β_GEN[g] * (sum(param.SCEN_PROB[r,s] * scen_profit[r,s] for r in set.PROFILES, s in set.SHIFTERS))
        - 0.5 * param.γ * (sum((prov_param.VOL_CONS[c] + vol_gen[c] + sum(prov_param.VOL_GEN[q,c] for q in filter(x->x != g, set.GENS)))^2 for c in set.CONTRACTS)))
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
        evaluation = init_evaluation(param)
        return sets, parameters, evaluation
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

    function calc_vol_limits(set, prov_param)
        vol_max_cons = prov_param.VOL_MAX_CONS
        vol_min_gen = prov_param.VOL_MIN_GEN
        percentage = 1.0
        for c in set.CONTRACTS
            vol_max_cons[c] = percentage * sum(prov_param.CAPACITY[g] for g in set.GENS)
            for g in set.GENS
                vol_min_gen[g,c] = -percentage * prov_param.CAPACITY[g]
            end
        end
        return vol_max_cons, vol_min_gen
    end

    function calc_payout(set, param, prov_param, fuel_price = 40.0)
        payout = prov_param.PAYOUT
        for r in set.PROFILES, s in set.SHIFTERS
            payout["HEAT_OPTION",r,s] = sum(param.DURATION[t] * max(0, prov_param.ENERGY_PRICE[t,r,s] - fuel_price) for t in set.TIME_BLOCKS)
        end
        return payout
    end

    function calc_gen_objective(set, param, prov_param, gen_model)
        gen_objective = Dict()
        for g in set.GENS
            gen_objective[g] = ((1 - param.β_GEN[g]) * (value(gen_model[g][:var_gen]) - (1 / param.α_GEN[g]) * (sum(param.SCEN_PROB[r,s] * value(gen_model[g][:scen_profit_aux][r,s]) for r in set.PROFILES, s in set.SHIFTERS)))
            + param.β_GEN[g] * (sum(param.SCEN_PROB[r,s] * value(gen_model[g][:scen_profit][r,s]) for r in set.PROFILES, s in set.SHIFTERS)))
        end
        return gen_objective
    end


    function intermediate_prov_param(inte_param::Dict, set::sets)
        ITERATION = inte_param["ITERATION"]
        ENERGY_PRICE = Dict((t,r,s) => inte_param["ENERGY_PRICE"][string((t,r,s))] for t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS)
        OPER_PROFIT = Dict((g,t,r,s)=> inte_param["OPER_PROFIT"][string((g,t,r,s))] for g in set.GENS, t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS)
        PAYOUT = Dict((c,r,s)=>inte_param["PAYOUT"][string((c,r,s))] for c in set.CONTRACTS, r in set.PROFILES, s in set.SHIFTERS)
        CONTRACT_PRICE = Dict(c => inte_param["CONTRACT_PRICE"][c] for c in set.CONTRACTS)
        CAPACITY = Dict(g => inte_param["CAPACITY"][g] for g in set.GENS)
        DEM_SERVED_FIX = Dict((t,r,s)=> inte_param["DEM_SERVED_FIX"][string((t,r,s))] for t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS)
        DEM_SERVED_FLEX = Dict((t,r,s)=> inte_param["DEM_SERVED_FLEX"][string((t,r,s))] for t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS)
        VOL_CONS = Dict(c => inte_param["VOL_CONS"][c] for c in set.CONTRACTS)
        VOL_GEN = Dict((g, c) => inte_param["VOL_GEN"][string((g,c))] for g in set.GENS, c in set.CONTRACTS)
        VOL_MIN_GEN = Dict((g, c) => inte_param["VOL_MIN_GEN"][string((g,c))] for g in set.GENS, c in set.CONTRACTS)
        VOL_MIN_CONS = Dict(c => inte_param["VOL_MIN_CONS"][c] for c in set.CONTRACTS)
        VOL_MAX_GEN = Dict((g, c) => inte_param["VOL_MAX_GEN"][string((g,c))] for g in set.GENS, c in set.CONTRACTS)
        VOL_MAX_CONS = Dict(c => inte_param["VOL_MAX_CONS"][c] for c in set.CONTRACTS)
        GEN_objective = Dict(g => inte_param["GEN_objective"][g] for g in set.GENS)

        return provisional_parameters(ITERATION, ENERGY_PRICE, OPER_PROFIT, PAYOUT, CONTRACT_PRICE, CAPACITY, DEM_SERVED_FIX,
        DEM_SERVED_FLEX, VOL_CONS, VOL_GEN, VOL_MIN_GEN, VOL_MIN_CONS, VOL_MAX_GEN, VOL_MAX_CONS, GEN_objective)
    end

    function solve_equilibrium(set, param, initial)
        FPATH = "/Users/hanshu/Desktop/response/data/Diff Failure/9per/intermediate_equ.json"
        finish = false
        avg_profit = nothing
        step = 1e-3
        ϵ = 1

            if initial ==1
                prov_param = init_provisional_parameters(set)
                open(FPATH, "w") do f
                    JSON.print(f, prov_param,4)
                end
                open("/Users/hanshu/Desktop/response/data/Diff Failure/9per/max_risk_imbalance.json", "w") do f
                    println(f, "Running ercot_equilibrium.jl ... contract percentage = 1.0 ")
                end
                open("/Users/hanshu/Desktop/response/data/Diff Failure/9per/GEN_objective.json", "w") do f
                    println(f, "Running ercot_equilibrium.jl ... contract percentage = 1.0 ")
                end
            else
                inte_param = JSON.parsefile(FPATH)
                prov_param = intermediate_prov_param(inte_param, set)
            end

        while !finish && prov_param.ITERATION < 50_000
            println("Iteration: ", prov_param.ITERATION )

            prov_param.CAPACITY = Dict(g => max(0, prov_param.CAPACITY[g] + step*prov_param.GEN_objective[g]/param.INV_COST[g]) for g in set.GENS)
            print("GENERATION CAPACITY: ",  "\t", prov_param.CAPACITY, "\n")
            prov_param.VOL_MAX_CONS, prov_param.VOL_MIN_GEN = calc_vol_limits(set, prov_param)
            print("prov_param.VOL_MAX_CONS: ",  "\t", prov_param.VOL_MAX_CONS, "\n")
            print("prov_param.VOL_MIN_GEN: ",  "\t", prov_param.VOL_MIN_GEN, "\n")
            dispatch_model = solve_dispatch(set, param, prov_param)
            print(prov_param.ITERATION,  "\t", prov_param.CAPACITY, "\n")
            print(prov_param.ITERATION, "\t", "ED Obj", "\t", objective_value(dispatch_model), "\n")
            for t in set.TIME_BLOCKS, r in set.PROFILES, s in set.SHIFTERS
                prov_param.DEM_SERVED_FIX[t,r,s] = value(dispatch_model[:dem_served_fix][t,r,s])
                prov_param.DEM_SERVED_FLEX[t,r,s] = value(dispatch_model[:dem_served_flex][t,r,s])
                prov_param.ENERGY_PRICE[t,r,s] = shadow_price(dispatch_model[:energy_balance][t,r,s])/param.DURATION[t]
                for g in set.GENS
                    prov_param.OPER_PROFIT[g,t,r,s] = shadow_price(dispatch_model[:max_prod][g,t,r,s])*param.AVAILABILITY[g,r,t]
                end
            end

            prov_param.PAYOUT = calc_payout(set, param, prov_param)

            if prov_param.ITERATION  == 1
                prov_param.CONTRACT_PRICE = Dict(c => sum(param.SCEN_PROB[r,s] * prov_param.PAYOUT[c,r,s] for r in set.PROFILES, s in set.SHIFTERS) for c in set.CONTRACTS)
                prov_param.VOL_CONS = Dict(c => prov_param.VOL_MAX_CONS[c] / max(length(set.CONTRACTS)) for c in set.CONTRACTS)
                prov_param.VOL_GEN = Dict((g, c) => -prov_param.VOL_CONS[c] / length(set.GENS) for g in set.GENS, c in set.CONTRACTS)
            end

            cons_model = nothing
            gen_model = nothing

            imbalance = 100
            while imbalance > ϵ

                println("Iteration: ", prov_param.ITERATION , " Imbalance: ", imbalance)
                # println("Gen objective: ", prov_param.GEN_objective)
                print(prov_param.ITERATION,  " \t", prov_param.CAPACITY, " \n")
                # print(prov_param.ITERATION, " \t", shadow_price(dispatch_model[:energy_balance][1,1,1]),  " \t", shadow_price(dispatch_model[:max_prod]["CCGT",94,2,10]), " \n")
                gen_model = solve_gen(set, param, prov_param)
                prov_param.VOL_GEN = Dict((g, c) => value(gen_model[g][:vol_gen][c]) for g in set.GENS, c in set.CONTRACTS)
                # print("GEN Var value: ", " \t", value(gen_model["CCGT"][:var_gen]), "\t", value(gen_model["w_CCGT"][:var_gen]), "\n")
                println("Generator volume: ", prov_param.VOL_GEN, "\n")

                cons_model = solve_cons(set, param, prov_param)
                println("consumer volume lower bound: ", prov_param.VOL_MIN_CONS, "\n")
                println("consumer volume higher bound: ", prov_param.VOL_MAX_CONS, "\n")
                println("consumer volume: ", value(cons_model[:vol_cons]["HEAT_OPTION"]), "\n")
                prov_param.VOL_CONS = Dict(c => value(cons_model[:vol_cons][c]) for c in set.CONTRACTS)
                # println("CONS Var value: ", "\t", value(cons_model[:var_cons]), "\n")
                println("Consumer volume: ", prov_param.VOL_CONS, "\n")

                difference = Dict(c => (prov_param.VOL_CONS[c] + sum(prov_param.VOL_GEN[g,c] for g in set.GENS)) for c in set.CONTRACTS)

                prov_param.CONTRACT_PRICE = Dict(c => prov_param.CONTRACT_PRICE[c] + param.γ * difference[c] for c in set.CONTRACTS)

                imbalance = maximum(map(x->abs(x), values(difference)))
                prov_param.GEN_objective = calc_gen_objective(set, param, prov_param, gen_model)
                println("GEN objective: ", " \t", prov_param.GEN_objective["CCGT"], " \t", prov_param.GEN_objective["w_CCGT"], "\n")
                # println("GEN Scen surplus: ", value(gen_model["CCGT"][:scen_profit][1,1]), "\n")
                # println("CONS Scen surplus: ", value(cons_model[:scen_surplus][1,1]), "\n")
            end

            prov_param.GEN_objective = calc_gen_objective(set, param, prov_param, gen_model)
            max_risk_imbalance = maximum(map(x->abs(x), values(prov_param.GEN_objective)))
            gen_model = solve_gen(set, param, prov_param)
            cons_model = solve_cons(set, param, prov_param)

            # print("F CONS Var value: ", "\t", value(cons_model[:var_cons]), "\n")
            # print("F GEN Var value: ", " \t", value(gen_model["CCGT"][:var_gen]), "\t", value(gen_model["w_CCGT"][:var_gen]), "\n")
            print("F Generstor volume: ", prov_param.VOL_GEN, "\n")
            print("F Consumer volume: ", prov_param.VOL_CONS, "\n")
            print("F GEN objective: ", " \t", prov_param.GEN_objective["CCGT"], " \t", prov_param.GEN_objective["w_CCGT"], "\n")

            open("/Users/hanshu/Desktop/response/data/Diff Failure/9per/max_risk_imbalance.json", "a") do f
                println(f, prov_param.ITERATION, " \t", max_risk_imbalance)
            end
            open("/Users/hanshu/Desktop/response/data/Diff Failure/9per/GEN_objective.json", "a") do f
                println(f, prov_param.ITERATION, " \t", prov_param.GEN_objective["CCGT"],  " \t", prov_param.GEN_objective["w_CCGT"])
            end
            #
            # println("Max risk imbalance: ", max_risk_imbalance)

            # step = max_risk_imbalance * 0.5 * 10^-6

            # Mays AMPL implementation calculates max_risk_imbalance as gen_objective times capacity.
            if max_risk_imbalance < 20000
                step = 1e-2
                ϵ = 0.3
                if max_risk_imbalance < 1_000
                    ϵ = 0.1
                    step = 0.02
                    # if max_risk_imbalance < 3
                    #     ϵ = 1e-1
                    #     step = 0.05
                    if max_risk_imbalance < 1e-1
                        finish = true
                        avg_profit = Dict(g => sum(param.SCEN_PROB[r,s] * value(gen_model[g][:scen_profit_aux][r,s]) for r in set.PROFILES, s in set.SHIFTERS) for g in set.GENS)
                    end
                    # end
                end
            end
            open(FPATH, "w") do f
                JSON.print(f, prov_param,4)
            end
            prov_param.ITERATION  += 1
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
                evalu.NET_REVENUE_UNIT1[g,r,s]  = sum(prov_param.OPER_PROFIT[g,t,r,s]  for t in set.TIME_BLOCKS)
                evalu.NET_REVENUE_UNIT2[g,r,s] = (sum((prov_param.ENERGY_PRICE[t,r,s] - param.MARG_COST[g])*param.DURATION[t]*value(dispatch_model[:prod][g,t,r,s]) for t in set.TIME_BLOCKS)) / prov_param.CAPACITY[g]
            end
        end

        for g in set.GENS
            evalu.AVE_REVENUE1[g] = sum(param.SCEN_PROB[r,s]*evalu.NET_REVENUE_UNIT1[g,r,s] for r in set.PROFILES, s in set.SHIFTERS)
            evalu.AVE_REVENUE2[g] = sum(param.SCEN_PROB[r,s]*evalu.NET_REVENUE_UNIT2[g,r,s] for r in set.PROFILES, s in set.SHIFTERS)
        end

        for g in set.GENS
            evalu.CAPACITY[g] = prov_param.CAPACITY[g]
            f(x) = evalu.AVE_REVENUE1[g]*sum(1/((1+x)^n) for n in 1:20) - param.INV_COST[g]*sum(1/((1.04)^n) for n in 1:20)
            evalu.IRR[g] = find_zero(f, (-1,1), Bisection())
        end

        h(x) = sum(param.SCEN_PROB[r,s]*(evalu.NET_REVENUE_UNIT1["w_CCGT",r,s]-evalu.NET_REVENUE_UNIT1["CCGT",r,s]) for r in set.PROFILES, s in set.SHIFTERS)*sum(1/((1+x)^n) for n in 1:20) - (param.INV_COST["w_CCGT"]-param.INV_COST["CCGT"])*sum(1/((1.04)^n) for n in 1:20)
        evalu.INC_IRR = find_zero(h, (-1,1), Bisection())

        return evalu
    end


    INPUT_FPATH = "/Users/hanshu/Desktop/response/code/fail_9.json"

    set, param, evalu = init_input(INPUT_FPATH)

    prov_param, avg_profit = solve_equilibrium(set, param, 1)

    evalu = solve_evaluation(set, param, prov_param, evalu)

    OUTPUT_FPATH = "/Users/hanshu/Desktop/response/data/Diff Failure/9per/out_equ_1.0.json"
    open(OUTPUT_FPATH, "w") do f
        JSON.print(f, evalu,4)
    end
