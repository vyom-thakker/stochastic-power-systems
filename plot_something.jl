using Plots
using FileIO

@userplot StackedArea

# a simple "recipe" for Plots.jl to get stacked area plots
# usage: stackedarea(xvector, datamatrix, plotsoptions)
@recipe function f(pc::StackedArea)
    x, y = pc.args
    n = length(x)
    y = cumsum(y, dims=2)
    seriestype := :shape

    # create a filled polygon for each item
    for c=1:size(y,2)
        sx = vcat(x, reverse(x))
        sy = vcat(y[:,c], c==1 ? zeros(n) : reverse(y[:,c-1]))
        @series (sx, sy)
    end
end

a = [1,1,1,1.5,2,3]
b = [0.5,0.6,0.4,0.3,0.3,0.2]
c = [2,1.8,2.2,3.3,2.5,1.8]
sNames = ["a","b","c"]
x = [2001,2002,2003,2004,2005,2006]

plotly()
stackedarea(x, [a b c], labels=reshape(sNames, (1,3)))

