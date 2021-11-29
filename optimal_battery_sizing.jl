#Optimal battery sizing
#Om Vinayaka

import Pkg; Pkg.add("CSV")
using JuMP
using LinearAlgebra
using GLPK
using Pkg
using CSV
using DataFrames; 

#Loading wind power,solar power and load data
dt=CSV.read("V:\\PHD\\Stochastic_Process\\Project_details\\Literature\\Battery_sizing_data.csv", DataFrame)

W_P=(dt[!,1]*1.0);S_P=(dt[!,2]*1.0);L_P=(dt[!,3]*0.1);

len=size(W_P,1)
D_P=zeros(len)

#Excess avalibaility or deficiency of renewable energy
for i in 1:len
   D_P[i]=round(W_P[i].+S_P[i].-L_P[i])
end


pi_z=zeros(w)
pi_m=zeros(w)
B_C=1


#Range of batteries
w=100;
for C in 1:w
function solve_opt(B_C)
 opt = Model(GLPK.Optimizer) 
 @variable(opt,0<=pi_s[i=1:1,j=1:(C.+1)])
 @variable(opt,0<=pi_z[i=1:1,j=1:(C)])   
 #@objective(opt, Min, (C*B_C))
   
  T_P=zeros(C.+1,C.+1)
  X=zeros(len)

#DTMC modeling and transition probability matrix        
   for i in 1:len 
    if i==1       
     X[i]=min(C,max(0,(C).+D_P[i]))
     a=convert(Int64,X[i])
     b=convert(Int64,(C))
     T_P[b.+1,a.+1]+=1
    else
     X[i]=min(C,max(0,X[i-1].+D_P[i]))
     a=convert(Int64,X[i])
     b=convert(Int64,X[i-1])
     T_P[b.+1,a.+1]+=1    
    end      
   end


#Verification for tranisition probability matrix
   sum=0;
   for i=1:(C+1)
    for j=1:(C+1)
      sum=sum.+T_P[i,j]
    end
   end

 T_PM=zeros(C+1,C+1)

   for i=1:(C+1)
    for j=1:(C+1)
     d=0
     for k=1:(C+1)
      d=d.+T_P[i,k]        
     end
    T_PM[i,j]=(T_P[i,j]./d)     
    end
  end


   sum1=zeros(C+1)
   for i=1:C+1   
    for j=1:C+1
     sum1[i]=sum1[i].+T_PM[i,j]
    end
   end
        
        
#Stationary probability matrix
        
  @constraint(opt,[j=1:C.+1],(pi_s[1,j]==(sum(pi_s[1,i].*T_PM[i,j] for i=1:(C.+1)))))
  @constraint(opt,[i=1:1],(sum(pi_s[i,j] for j=1:((C.+1))))==1)
  #@constraint(opt,[i=1:1,j=1:1],pi[i,j]<=1.0)
  #@constraint(opt,[i=1:1],(sum(pi[i,j] for j=((24*(C+1))-(C-1)):24*(C+1)))>=0.999)
        
 optimize!(opt)
    return JuMP.value.(pi_s),JuMP.value.(pi_s[1,1]),JuMP.value.(pi_s[1,C.+1])         
end
(pi_st,pi_z[C],pi_m[C]) = solve_opt(B_C);    
end 

opts=zeros(w)
optsize=0;
for C in 1:w
 if pi_z[C]<=0.1 
  if optsize<=C & optsize==0
  optsize=C-1          
  end                 
 end      
end


#optsize=min(optsi)
#println("wind power data: ",W_P, " MW")
#println("\n")

#println("solar power data: ",S_P, " MW")
#println("\n")

#println("load power data: ",L_P, " MW")
#println("\n")

#println("Diff power data: ",D_P, " MW")
#println("\n")

#println("X Values: ",X, " MW")
#println("\n")

#println("Sum Values: ",sum, " Check")
#println("\n")


#println("Sum Values: ",sum1, " Check")
#println("\n")


#println("Transistion Probability Matrix: ",T_PM, " Matrix")
#println("\n")

#println("Length of data: ",len)
#println("\n")

#println("Stationary Probability: ",pi_st)
#println("\n")

#println("Stationary Probability at state zero for all batteries: ",pi_z)
#println("\n")

println("Stationary Probability at state zero for all batteries: ",pi_z)
println("\n")

println("Stationary Probability at Max state for all batteries: ",pi_m)
println("\n")

println("optimal battery sizes: ",optsize)
println("\n")


#CSV.read("V:\\PHD\\Stochastic_Process\\Project_details\\Literature\\Battery_sizing_data_1.csv");


