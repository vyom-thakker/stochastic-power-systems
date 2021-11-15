using JuMP
using LinearAlgebra
#using DiffOpt, Plots
using GLPK
using ChainRulesCore
g_max = [1000, 1000]
g_min = [0, 300]
c_g = [50, 100]
c_g0 = [1000, 0]
c_w = 50
d = 1500
w_f = 200
function solve_ed(g_max, g_min, c_g, c_w, d, w_f)
    ed = Model(GLPK.Optimizer)
    @variable(ed, 0 <= g[i = 1:2] <= g_max[i])
    @variable(ed, 0 <= w <= w_f)
    @objective(ed, Min, c_g' * g + c_w * w)
    @constraint(ed, [i = 1:2], g[i] <= g_max[i])
    @constraint(ed, [i = 1:2], g[i] >= g_min[i]) 
    @constraint(ed, w <= w_f)
    @constraint(ed, sum(g) + w == d)
    optimize!(ed)
    return value.(g), value(w), w_f - value(w), objective_value(ed)
end
(g_opt, w_opt, ws_opt, obj) = solve_ed(g_max, g_min, c_g, c_w, d, w_f);

println("Dispatch of Generators: ", g_opt, " MW")
println("Dispatch of Wind: ", w_opt, " MW")
println("Wind spillage: ", w_f - w_opt, " MW")
println("\n")
println("Total cost: ", obj, "\$")

g_max = [1000, 1000]
g_min = [0,30]
c_g = [50, 100]
c_g0 = [1000, 0]
c_w = 1
d = [1500 1500]
w_f = 200
function solve_ed(g_max, g_min, c_g, c_w, d, w_f)
    ed = Model(GLPK.Optimizer)
    @variable(ed, 0 <= g[i=1:2,j=1:2] <= g_max[i])
    @variable(ed, 0 <= w[i=1:1,j=1:2] <= w_f)
    @objective(ed, Min, sum((g[i,j] *c_g[i,1]) + (c_w* w[1,j]) for j=1:2, i=1:2))
   # @objective(
    #    model,
     #   Min,
      #  sum((Cp[g] * p[g, t]) + (Cnl[g] * u[g, t]) for g in unit_codes, t in 1:n_periods),
    #)

    @constraint(ed,[i=1:2,j=1:2], g[i,j] <= g_max[i])
    @constraint(ed,[i=1:2,j=1:2], g[i,j] >= g_min[i]) 
    @constraint(ed,[i=1:1,j=1:2],w[i,j] <= w_f)
    @constraint(ed,[j=1:2],(g[1,j]+ g[2,j] + w[1,j] == d[1,j]))
    optimize!(ed)
    return JuMP.value.(g), JuMP.value.(w), objective_value(ed),(g[i,j] for i=1:2, j=1:2)
end
(g_opt, w_opt, obj,un) = solve_ed(g_max, g_min, c_g, c_w, d, w_f);

#display(un)
println("Dispatch of Generators: ",g_opt, " MW")
println("Dispatch of Wind: ", w_opt, " MW")
#println("Wind spillage: ", w_f - w_opt, " MW")
println("\n")
println("Total cost: ", obj, "\$")

g_max = [1000, 1000]
g_min = [0,300]
c_g = [50, 100]
c_g0 = [1000, 0]
c_w = 50
d = [1500 1000]
w_f = 200
function solve_ed(g_max, g_min, c_g, c_w, d, w_f)
    ed = Model(GLPK.Optimizer)
    @variable(ed, 0 <= g[i=1:2,j=1:2] <= g_max[i])
    @variable(ed, 0 <= w[i=1:1,j=1:2] <= w_f)
    @objective(ed, Min, sum(sum((g[i,j] *c_g[i,1]) for i=1:2)+ (c_w* w[1,j]) for j=1:2))
   # @objective(
    #    model,
     #   Min,
      #  sum((Cp[g] * p[g, t]) + (Cnl[g] * u[g, t]) for g in unit_codes, t in 1:n_periods),
    #)

    @constraint(ed,[i=1:2,j=1:2], g[i,j] <= g_max[i])
    @constraint(ed,[i=1:2,j=1:2], g[i,j] >= g_min[i]) 
    @constraint(ed,[i=1:1,j=1:2],w[i,j] <= w_f)
    @constraint(ed,[j=1:2],(g[1,j]+ g[2,j] + w[1,j] == d[1,j]))
    optimize!(ed)
    return JuMP.value.(g), JuMP.value.(w), objective_value(ed)
end
(g_opt, w_opt, obj) = solve_ed(g_max, g_min, c_g, c_w, d, w_f);

#display(un)
println("Dispatch of Generators: ",g_opt, " MW")
println("Dispatch of Wind: ", w_opt, " MW")
#println("Wind spillage: ", w_f - w_opt, " MW")
println("\n")
println("Total cost: ", obj, "\$")

g_max = [1000, 1000]
g_min = [300,300]
c_g = [100, 100]
c_g0 = [1000, 0]
c_w = 50
d = [1000 1000 1300 1000 1000 1000 1300]
w_f = 200
UTg=[4,0]
UDg=[4,0]
#u= Int8
function solve_ed(g_max, g_min, c_g, c_w, d, w_f)
    ed = Model(GLPK.Optimizer)
    @variable(ed, 0 <= w[i=1:1,j=1:7] <= w_f)
    @variable(ed,u[i=1:2,j=1:7], Bin)
    @variable(ed, 0 <= g[i=1:2,j=1:7] <= g_max[i])
    @objective(ed, Min, sum(sum(((g[i,j] *c_g[i,1])+(u[i,j]*c_g0[i,1])) for i=1:2)+ (c_w* w[1,j]) for j=1:7))
    @constraint(ed,[i=1:2,j=1:7], g[i,j] <= g_max[i]*u[i,j])
    @constraint(ed,[i=1:2,j=1:7], g[i,j] >= g_min[i]*u[i,j]) 
    @constraint(ed,[i=1:1,j=1:7],w[i,j] <= w_f)
    @constraint(ed,[j=1:7],(g[1,j]+ g[2,j] + w[1,j] == d[1,j]))
    #@constraint(ed,[i=1:2],u[i,0]==0)
    @constraint(ed,[j=2:7,i=1:2,s=j+1:min(7,j.+UTg[i,1]-1)],u[i,s]>=u[i,j]-u[i,j-1])
    @constraint(ed,[j=1:1,i=1:2,s=j+1:min(7,j.+UTg[i,1]-1)],u[i,s]>=u[i,j]-0)
    @constraint(ed,[j=2:7,i=1:2,s=j+1:min(7,j.+UDg[i,1]-1)],(1-u[i,s])>=u[i,j-1]-u[i,j])
    @constraint(ed,[j=1:1,i=1:2,s=j+1:min(7,j.+UDg[i,1]-1)],(1-u[i,s])>=0-u[i,j])
    optimize!(ed)
    return JuMP.value.(g), JuMP.value.(w), objective_value(ed),JuMP.value.(u)
end
(g_opt, w_opt, obj,u) = solve_ed(g_max, g_min, c_g, c_w, d, w_f);

#display(un)
println("Dispatch of Generators: ",g_opt, " MW")
println("Binary: ",u, " Bin")
println("Dispatch of Wind: ", w_opt, " MW")
#println("Wind spillage: ", w_f - w_opt, " MW")
println("\n")
println("Total cost: ", obj, "\$")

g_max = [1000, 1000]
g_min = [300,300]
c_g = [100, 100]
c_g0 = [1000, 0]
c_w = 50
d = [1000 1000 1300 1000 1000 1000 1300]
w_f = 200
UTg=[4,0]
UDg=[4,0]
#u= Int8
function solve_ed(g_max, g_min, c_g, c_w, d, w_f)
    ed = Model(GLPK.Optimizer)
    @variable(ed, 0 <= w[i=1:1,j=1:7] <= w_f)
    @variable(ed,u[i=1:2,j=1:7], Bin)
    @variable(ed, 0 <= g[i=1:2,j=1:7] <= g_max[i])
    @objective(ed, Min, sum(sum(((g[i,j] *c_g[i,1])+(u[i,j]*c_g0[i,1])) for i=1:2)+ (c_w* w[1,j]) for j=1:7))
    @constraint(ed,[i=1:2,j=1:7], g[i,j] <= g_max[i]*u[i,j])
    @constraint(ed,[i=1:2,j=1:7], g[i,j] >= g_min[i]*u[i,j]) 
    @constraint(ed,[i=1:1,j=1:7],w[i,j] <= w_f)
    @constraint(ed,[j=1:7],(g[1,j]+ g[2,j] + w[1,j] == d[1,j]))
    #@constraint(ed,[i=1:2],u[i,0]==0)
    @constraint(ed,[j=2:7,i=1:2,s=j+1:min(7,j.+UTg[i,1]-1)],u[i,s]>=u[i,j]-u[i,j-1])
    @constraint(ed,[j=1:1,i=1:2,s=j+1:min(7,j.+UTg[i,1]-1)],u[i,s]>=u[i,j]-0)
    @constraint(ed,[j=2:7,i=1:2,s=j+1:min(7,j.+UDg[i,1]-1)],(1-u[i,s])>=u[i,j-1]-u[i,j])
    @constraint(ed,[j=1:1,i=1:2,s=j+1:min(7,j.+UDg[i,1]-1)],(1-u[i,s])>=0-u[i,j])
    optimize!(ed)
    return JuMP.value.(g), JuMP.value.(w), objective_value(ed),JuMP.value.(u)
end
(g_opt, w_opt, obj,u) = solve_ed(g_max, g_min, c_g, c_w, d, w_f);

#display(un)
println("Dispatch of Generators: ",g_opt, " MW")
println("Binary: ",u, " Bin")
println("Dispatch of Wind: ", w_opt, " MW")
#println("Wind spillage: ", w_f - w_opt, " MW")
println("\n")
println("Total cost: ", obj, "\$")


