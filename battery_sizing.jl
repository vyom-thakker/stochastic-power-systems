using JuMP
using LinearAlgebra
#using DiffOpt, Plots
using GLPK
using ChainRulesCore
using CSV
using DataFrames
using FileIO


# get data
df = CSV.File("./Battery_sizing_data.csv", header=0, delim=',') |> DataFrame 
