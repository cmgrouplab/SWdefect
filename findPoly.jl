using LightGraphs
using CSV
using DelimitedFiles 

result = Array{Union{Float64, Int64}}(undef,20,4)
k = 1
for pValue in 0.05:0.05:1
    data = CSV.read("Np_216_random_p_0.05-1/p_$(pValue)/1/positionP.csv"; header=false)
    result[k,1] = pValue
    numP = length(data[1])
    numPt = numP ÷ 2

    println("Num of P: $numP, Num of Pt: $numPt")

    g = SimpleGraph(numP+numPt)

    for p in 1:numP
        add_edge!(g, p, data[p, 4])
        add_edge!(g, p, data[p, 5] + numP)  # Label for Pt_i = Pt_i + numP
        add_edge!(g, p, data[p, 6] + numP)
    end

    g = DiGraph(g)

    cycles_below_6 = simplecycles_limited_length(g, 6)

# for c in cycles_below_6
#     if length(c) == 6
#         println(c)
#     end
# end

    for l in 4:6
        n = count(c -> length(c) == l, cycles_below_6)
        println("Cycle len $l count: $(n ÷ 2)")
        result[k,l-2] = n ÷ 2 
    end
    global k += 1
end

writedlm("Np_216_random_poly.csv", result, ',')