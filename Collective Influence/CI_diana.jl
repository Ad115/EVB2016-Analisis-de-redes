
Pkg.clone("https://github.com/scheinerman/SimpleGraphs.jl.git")
using SimpleGraphs

red = readdlm("red_final_enfermos.dat",'\t')
g = StringGraph()
for n in 2:Int64(size(red,1))
    add!(g,red[n,1],red[n,3])
end

function CI(g,n::Int64)
    #n < diameter(g) || error("n debe ser menor que el diámetro de la red")
    CI = zeros(NV(g))
    for v in 1:NV(g)
        vec = neighbors(g,vlist(g)[v])
        vec2 = vec
        vec = vcat(vec,vlist(g)[v])
        todos = vec
        distn = Array(Int64,0)
        for l in 1:(n-1)
            todos = union(vec,vec2)
            for m in vec2
                vec2 = union(vec2,neighbors(g,m))
            end
            distn = union(setdiff(todos,vec2),setdiff(vec2,todos))
        end
        Î´ = zeros(length(distn))
        for m in 1:length(distn)
            Î´[m] = (deg(g,distn[m]))-1
        end
        Ïƒ = sum(Î´)
        CI[v] = Ïƒ * (deg(g,vlist(g)[v])-1)
    end
    CI = hcat(vlist(g),CI)
end


collective = CI(g,2)
writedlm("CI_Final_network_4_0.99_DPI_0.1.dat", collective)
add!(g,red[n,1],red[n,2])
