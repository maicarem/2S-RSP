function _write_gurobi_log(filename, master::Model, method_name)
    open("results/$(filename).txt", "w") do io
        println(io, "################ $(method_name) ###############")
        # Model initialization
        println(io, "Model Summary")
        println(io, master)
        # Model attribute
        println(io, "Model Attribute")
        obj = objective_value(master)
        gap = relative_gap(master)
        run_time = solve_time(master)
        
        name = ["ObjVal", "Gap", "Run time"]
        for (i,j) in zip(name, [obj, gap, run_time])
            println(io, "$(i)--> $(j)")
        end
    end
end