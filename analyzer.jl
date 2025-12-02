using Plots

Es = Float64[]
ERPA0_1s = Float64[]
ERPA1_1s = Float64[]
ERPA2_1s = Float64[]
ERPA0_2s = Float64[]
ERPA1_2s = Float64[]
ERPA2_2s = Float64[]
ERPA0_3s = Float64[]
ERPA1_3s = Float64[]
ERPA2_3s = Float64[]
ERPA0_4s = Float64[]
ERPA1_4s = Float64[]
ERPA2_4s = Float64[]

rs = 1.0:0.1:4.5

for r in rs
    erpa0 = false
    erpa1 = false
    erpa2 = false
    filename = "n2/PNOF7/n2_$r.out"
    println(filename)
    E = 0.0
    open(filename, "r") do f
        for line in eachline(f)
             if occursin("Final NOF", line)
                 println("✅ ENCONTRADO: $line)")
		 E = parse(Float64,split(line)[6])
 		 push!(Es, E)
             end
             if occursin("Excitation energies ERPA0", line)
		 readline(f)
		 readline(f)
		 readline(f)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA0_1s, E + exc)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA0_2s, E + exc)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA0_3s, E + exc)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA0_4s, E + exc)
             end
             if occursin("Excitation energies ERPA/", line)
                 readline(f)
                 readline(f)
                 readline(f)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA1_1s, E + exc)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA1_2s, E + exc)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA1_3s, E + exc)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA1_4s, E + exc)
             end
             if occursin("Excitation energies ERPA2", line)
                 readline(f)
                 readline(f)
                 readline(f)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA2_1s, E + exc)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA2_2s, E + exc)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA2_3s, E + exc)
		 line2 = readline(f)
		 exc = parse(Float64,split(line2)[4])/27.2114
		 push!(ERPA2_4s, E + exc)
             end

         end
     end
end

plot(rs, Es, marker=:circle)
plot!(rs, ERPA0_1s, marker=:circle)
plot!(rs, ERPA0_2s, marker=:circle)
plot!(rs, ERPA0_3s, marker=:circle)
plot!(rs, ERPA0_4s, marker=:circle)
