using DataFrames
using CSV
function map_output()
    directory = "results"
    combined_text = ""
    files = readdir(directory)

    for file in files
        file_path = joinpath(directory, file)
        content = read(file_path, String)
        println(content)
        if occursin("time_limit: 3600",content)
            combined_text *= content
            combined_text *= "\n"
        end
        
    end

    output_file = "combined_output.txt"
    write(output_file, combined_text)
    @info("Done")

    # Read the input file
    input_file = "combined_output.txt"
    input_text = read(input_file, String)

    # Split the input text into lines
    lines = split(input_text, "\n")

    # Create an empty DataFrame
    df = DataFrame()

    # Extract the keys from the first line
    keys = split(lines[1], ",")
    # keys = hcat(map(strip, keys))
    value_dict = Dict()
    
    for each_key in keys
        a = split(each_key, ":")[1]
        if a == "name"
            value_dict[a] = []
        else
            b = split(a, " ")[2]
            if !(b == "")
                value_dict[b] = []
            end
        end
    end
    @show value_dict
    total = 0
    # Extract the values from the remaining lines
    for each_line in lines[1:end-1]
        line0 = split(each_line, ",")
        keys0 = hcat(map(strip, line0))
        list_keys0 = [split(each_key0, ":")[1] for each_key0 in keys0]
        if "node_count" ∉ list_keys0
            push!(value_dict["node_count"], 0)
            # total += 1 
        end
        # result_dict["gap"] = (result_dict["upper_bound"] - result_dict["lower_bound"])/result_dict["upper_bound"]

        ub, lb = 0, 0

        for each_key0 in keys0
            if !(each_key0 == "")
                
                key0 = split(each_key0, ":")[1]
                val0 = split(each_key0, ":")[2]

                if key0 ∈ ["name","algorithm", "split_sp0", "benders","transformation", "one_cut", "uc", "timestamp"]
                    push!(value_dict[key0], val0)
                else
                    push!(value_dict[key0], parse(Float64, val0))
                end

                if key0 == "upper_bound"
                    ub = value_dict["upper_bound"][end]
                elseif key0 == "lower_bound"
                    lb = value_dict["lower_bound"][end]
                end
            end
        end
        
        if "gap" ∉ list_keys0
            push!(value_dict["gap"], (ub - lb)/ub)
        end
    end
    df = DataFrame(value_dict)
    CSV.write("combined_test_1.csv", df)
    # @show value_dict
    # @show vcat(DataFrame.(value_dict)...)
    # CSV.write("combined_test_1.csv", value_dict)
    # @info value_dict
    # df = DataFrame(value_dict)
    
end

map_output()