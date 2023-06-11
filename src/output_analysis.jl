using DataFrames
using CSV
function map_output()
    directory = "results"
    combined_text = ""
    files = readdir(directory)

    for file in files
        file_path = joinpath(directory, file)
        content = read(file_path, String)
        combined_text *= content
        combined_text *= "\n"
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

    for each_line in lines
        line0 = split(each_line, ",")
        keys0 = hcat(map(strip, line0))
        for each_key0 in keys0
            # @show each_key0
            if !(each_key0 == "")
                key0 = split(each_key0, ":")[1]
                val0 = split(each_key0, ":")[2]
                # @show typeof(val0)

                if key0 ∈ ["name","algorithm", "split_sp0", "benders","transformation", "one_cut", "uc", "timestamp"]
                    push!(value_dict[key0], val0)
                else
                    push!(value_dict[key0], parse(Float64, val0))
                end
            end
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