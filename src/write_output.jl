function _write_gurobi_log(master, method_name, instance_name, pars, n, V_tilt, time_benders)
    if method_name == "ilp"
        open("ilp_log.txt", "a") do io
            println(io,"")
            run_time = solve_time(master)
            obj = objective_value(master)
            gap = relative_gap(master)
        
            for j in [instance_name, n, length(V_tilt), method_name, pars.transformation, run_time, obj, gap, now()]
                print(io, "$(j), ")
            end
        end
    end
end

function _write_log_bruteforce(instance_name, n, V_tilt, method_name, objval_3, objsol_3, objtime_3, objval_4, objsol_4, objtime_4, objval_5, objsol_5, objtime_5)
    open("bruteforce_log.txt", "a") do io
        println(io,"")
        for j in [instance_name, n, length(V_tilt), method_name, objval_3, objsol_3, objtime_3, objval_4, objsol_4, objtime_4, objval_5, objsol_5, objtime_5, now()]
            print(io, "$(j), ")
        end
        @info "Write brute force 2: DONE"
    end
end


function _write_log_benders(instance_name, n, V_tilt, method_name, lower_bound, global_ub, total_time, pars)
    open("benders_log.txt", "a") do io
        println(io,"")
        for j in [instance_name, n, length(V_tilt), method_name, pars.one_cut, pars.transformation, pars.split_sp0, total_time, lower_bound, global_ub, now()]
            print(io, "$(j), ")
        end
    end
end

function _write_log_bbc(instance_name, master, n, V_tilt, method_name, pars)
    open("bbc_log.txt", "a") do io
        println(io, "")
        
        run_time = solve_time(master)
        obj = objective_value(master)
        gap = relative_gap(master)

        for j in [instance_name, n, length(V_tilt), method_name, pars.one_cut, pars.transformation, pars.split_sp0, run_time, obj, gap, now()]
            print(io, "$(j), ")
        end
    end
end


function write_ouput(pars, name, result_dict0, MainPar)
    name = split(split(name,"/")[end], ".")[1]
    pars0 = Dict(key=>getfield(pars, key) for key ∈ fieldnames(MainPar))
    open("results/$(name)_$(result_dict0["algorithm"])_$(pars.alpha)_$(pars.one_cut)_$(pars.transformation).txt", "w") do io
        print(io, "name: $(name), ")
        for i in keys(pars0)
            print(io, "$(i): $(pars0[i]), ")
        end
        for i in keys(result_dict0)
            print(io, "$(i): $(result_dict0[i]), ")
        end
    end
end


