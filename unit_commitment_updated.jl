#Om vinayaka
#ISE-7300 Stochastic Process project

using JuMP
using LinearAlgebra
using GLPK

#Conventional generator Min and Max values
g_max = [20, 50, 100]
g_min = [5,20, 0]

#Conventional generator per MW cost
c_g = [27.7, 35.5, 10000]

#Conventional generator no load cost
c_g0 = [120, 100, 10000]

#Wind generator per MW cost
c_w = 10

#Solar generator per MW cost
c_s = 10

# Load demand for one day or 24 hours
d = [76.07 90.33 33.78 39.41 34.96 36.73 32.5 29.65 35.89 45.5 41.9 38.6 43.16 65.56 45.02 35.63 40.31 56.85 74.55 90 57.22 49.6 59.49 83.5]

# Wind power maximum value at each hour for 24 hours
w_f = [40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40]

#Solar power maximum value at each hour for 24 hours
s_f = [30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30]

#Minimum uptime and downtime for conventional generation
UTg=[2,4,0]
UDg=[2,4,0]

#Minimum up and down ramp rate for conventional generator
UR=[20,50,100]
DR=[20,50,100]

function solve_ed(g_max, g_min, c_g, c_w, d, w_f, s_f)
    ed = Model(GLPK.Optimizer)
    
#Wind and Solar generation variable for each hour    
    @variable(ed, 0 <= w[i=1:1,j=1:24] <= w_f[i,j])
    @variable(ed, 0 <= s[i=1:1,j=1:24] <= s_f[i,j])
    
#Unit commitment binary variable for conventional generator    
    @variable(ed,u[i=1:3,j=1:24], Bin)
    
#Conventional generation varaible for each hour     
    @variable(ed, 0 <= g[i=1:3,j=1:24] <= g_max[i])
    
#Objective Minimize generation production cost    
    @objective(ed, Min, sum(sum(((g[i,j] *c_g[i])+(u[i,j]*c_g0[i])) for i=1:3)+ (c_w* w[1,j])+(c_s*s[1,j]) for j=1:24))
    
#Conventional generation constraints    
    @constraint(ed,[i=1:3,j=1:24], g[i,j] <= g_max[i]*u[i,j])
    @constraint(ed,[i=1:3,j=1:24], g[i,j] >= g_min[i]*u[i,j]) 
    
#Wind and Solar generation constraints    
    @constraint(ed,[i=1:1,j=1:24],w[i,j] <= w_f[1,j])
    @constraint(ed,[i=1:1,j=1:24],s[i,j] <= s_f[1,j])
    
#Power balance equation constraint   
    @constraint(ed,[j=1:24],(sum(g[i,j] for i=1:3)+ w[1,j]+ s[1,j] == d[1,j]))

#Minimum uptime and downtime constraint   
    @constraint(ed,[j=2:24,i=1:3,s=j+1:min(24,j.+UTg[i,1]-1)],u[i,s]>=u[i,j]-u[i,j-1])
    @constraint(ed,[j=1:1,i=1:3,s=j+1:min(24,j.+UTg[i,1]-1)],u[i,s]>=u[i,j]-0)
    @constraint(ed,[j=2:24,i=1:3,s=j+1:min(24,j.+UDg[i,1]-1)],(1-u[i,s])>=u[i,j-1]-u[i,j])
    @constraint(ed,[j=1:1,i=1:3,s=j+1:min(24,j.+UDg[i,1]-1)],(1-u[i,s])>=0-u[i,j])

#Ramp rate up and down constraint
    @constraint(ed,[i=1:3,j=2:24], g[i,j]-g[i,j-1]<=UR[i])
    @constraint(ed,[i=1:3,j=1:1], g[i,j]-0<=UR[i])
    @constraint(ed,[i=1:3,j=2:24], g[i,j-1]-g[i,j]<=DR[i])
    @constraint(ed,[i=1:3,j=1:1], 0-g[i,j]<=DR[i])
    
    optimize!(ed)
    return JuMP.value.(g), JuMP.value.(w),JuMP.value.(s), objective_value(ed),JuMP.value.(u)
end
(g_opt, w_opt, s_opt, obj,u) = solve_ed(g_max, g_min, c_g, c_w, d, w_f, s_f);



#display(un)
println("Dispatch of Coventional Generators: ",g_opt, " MW")
println("\n")
println("Dispatch of Wind Generator: ",w_opt, " MW")
println("\n")
println("Dispatch of Solar Generator: ",s_opt, " MW")
println("\n")
println("Binary:",u, " Bin")
println("\n")
println("Total cost: ", obj, "\$")


