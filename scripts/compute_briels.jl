using DrWatson
@quickactivate

using Sausage: solve_velocity_briels
using Logging
using Dates
using Formatting
using DelimitedFiles

function fmt(level, _module, group, id, file, line)
    return (:blue, format("{:<23}:", string(Dates.now())), "")
end

function main()
    global_logger(ConsoleLogger(meta_formatter=fmt))

    eb = 5e5:1e5:4e6
    a = similar(eb)
    v = similar(eb)

    for (i, eb1) in enumerate(eb)
        @info "[$i / $(length(eb))] eb = $eb1"
        v[i], a[i] = solve_velocity_briels(eb1)        
        
        open(datadir("briels.dat"), append=true) do fout
            println(fout, "$(eb1)\t$(v[i])\t$(a[i])")
        end
        @show eb1, v[i], a[i]
    end    
end

isinteractive() || main()
