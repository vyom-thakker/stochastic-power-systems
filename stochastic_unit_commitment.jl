using JuMP
using LinearAlgebra
#using DiffOpt, Plots
using GLPK
using ChainRulesCore
using CSV
using DataFrames
using FileIO

#initial wind speed
w_p_i=[0 0 0 0 1 0 0 0]
#initial cloud night
s_p_i=[0 0 1 0 0 0 0 0]

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
d = [35.89 45.5 41.9 38.6 43.16 65.56 45.02 35.63 40.31 56.85 74.55 90 57.22 49.6 59.49 83.5 76.07 90.33 33.78 39.41 34.96 36.73 32.5 29.65]

# Wind power maximum value at each hour for 24 hours
w_f = [10 40]
w_pt_df = CSV.File("./wind_transitionMat.csv", header=0, delim=',') |> DataFrame 
w_pt=Matrix(w_pt_df)

#Solar power maximum value at each hour for 24 hours
s_f = [10 30]
s_pt_df = CSV.File("./solar_transitionMat.csv", header=0, delim=',') |> DataFrame
s_pt=Matrix(s_pt_df)


#Minimum uptime and downtime for conventional generation
UTg=[2,4,0]
UDg=[2,4,0]

#Minimum up and down ramp rate for conventional generator
UR=[20,50,100]
DR=[20,50,100]

function solve_ed(g_max, g_min, c_g, c_w, d, w_f, s_f)
    ed = Model(GLPK.Optimizer)
    
#Wind and Solar generation variable for each hour    
    @variable(ed, 0 <= w[i=1:1,j=1:24] <= w_f[1])
    @variable(ed, 0 <= s[i=1:1,j=1:24] <= s_f[1])
    
#Unit commitment binary variable for conventional generator    
    @variable(ed,u[i=1:3,j=1:24], Bin)
    
#Conventional generation varaible for each hour     
    @variable(ed, 0 <= g[i=1:3,j=1:24] <= g_max[i])
    
#Objective Minimize generation production cost    
    @objective(ed, Min, sum(sum(((g[i,j] *c_g[i])+(u[i,j]*c_g0[i])) for i=1:3)+ (c_w* w[1,j])+(c_s*s[1,j]) for j=1:24))
    
#Conventional generation constraints    
    @constraint(ed,[i=1:3,j=1:24], g[i,j] <= g_max[i]*u[i,j])
    @constraint(ed,[i=1:3,j=1:24], g[i,j] >= g_min[i]*u[i,j]) 
    @constraint(ed,[i=1:1,j=12:20],s[i,j] ==0) 
    
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


#stochastic problem

function wind_pow_corr(p_t)
    w_pow=[0 0.06 0.18 0.56 0.85 1 0.96 0.87];
    return w_f[2]*sum(w_pow.*p_t);
end

function solar_pow_corr(p_t)
    s_pow=[1 1 0.94 0.85 0.7 0.5 0.3 0.1];
    return s_f[2]*sum(s_pow.*p_t);
end

function solve_st(g_max, g_min, c_g, c_w, d, w_f, s_f,w_p_i,s_p_i)
    st = Model(GLPK.Optimizer)
    
#Wind and Solar generation variable for each hour    
    @variable(st, 0 <= w[i=1:1,j=1:24] <= wind_pow_corr(w_p_i*w_pt^j))
    @variable(st, 0 <= s[i=1:1,j=1:24] <= solar_pow_corr(s_p_i*s_pt^j))   
    @constraint(st,[i=1:1,j=12:20],s[i,j] ==0) 

#Unit commitment binary variable for conventional generator    
    @variable(st,u[i=1:3,j=1:24], Bin)
    
#Conventional generation varaible for each hour     
    @variable(st, 0 <= g[i=1:3,j=1:24] <= g_max[i])
    
#Objective Minimize generation production cost    
    @objective(st, Min, sum(sum(((g[i,j] *c_g[i])+(u[i,j]*c_g0[i])) for i=1:3)+ (c_w* w[1,j])+(c_s*s[1,j]) for j=1:24))
    
#Conventional generation constraints    
    @constraint(st,[i=1:3,j=1:24], g[i,j] <= g_max[i]*u[i,j])
    @constraint(st,[i=1:3,j=1:24], g[i,j] >= g_min[i]*u[i,j]) 
    
    
#Power balance equation constraint   
    @constraint(st,[j=1:24],(sum(g[i,j] for i=1:3)+ w[1,j]+ s[1,j] == d[1,j]))

#Minimum uptime and downtime constraint   
    @constraint(st,[j=2:24,i=1:3,s=j+1:min(24,j.+UTg[i,1]-1)],u[i,s]>=u[i,j]-u[i,j-1])
    @constraint(st,[j=1:1,i=1:3,s=j+1:min(24,j.+UTg[i,1]-1)],u[i,s]>=u[i,j]-0)
    @constraint(st,[j=2:24,i=1:3,s=j+1:min(24,j.+UDg[i,1]-1)],(1-u[i,s])>=u[i,j-1]-u[i,j])
    @constraint(st,[j=1:1,i=1:3,s=j+1:min(24,j.+UDg[i,1]-1)],(1-u[i,s])>=0-u[i,j])

#Ramp rate up and down constraint
    @constraint(st,[i=1:3,j=2:24], g[i,j]-g[i,j-1]<=UR[i])
    @constraint(st,[i=1:3,j=1:1], g[i,j]-0<=UR[i])
    @constraint(st,[i=1:3,j=2:24], g[i,j-1]-g[i,j]<=DR[i])
    @constraint(st,[i=1:3,j=1:1], 0-g[i,j]<=DR[i])
    
    optimize!(st)
    return JuMP.value.(g), JuMP.value.(w),JuMP.value.(s), objective_value(st),JuMP.value.(u)
end
(g_opt, w_opt, s_opt, obj,u) = solve_ed(g_max, g_min, c_g, c_w, d, w_f, s_f);

println("Dispatch of Conventional Generators: ",g_opt, " MW")
println("\n")
println("Dispatch of Wind Generator: ",w_opt, " MW")
println("\n")
println("Dispatch of Solar Generator: ",s_opt, " MW")
println("\n")
println("Binary:",u, " Bin")
println("\n")
println("Total cost: ", obj, "\$")

(g_opt, w_opt, s_opt, obj,u) = solve_st(g_max, g_min, c_g, c_w, d, w_f, s_f,w_p_i,s_p_i);
println("Dispatch of Coventional Generators: ",g_opt, " MW")
println("\n")
println("Dispatch of Wind Generator: ",w_opt, " MW")
println("\n")
println("Dispatch of Solar Generator: ",s_opt, " MW")
println("\n")
println("Binary:",u, " Bin")
println("\n")
println("Total cost: ", obj, "\$")
