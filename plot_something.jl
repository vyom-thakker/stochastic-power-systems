using DataFrames
using VegaLite, VegaDatasets
using FileIO
a = [1,1,1,1.5,2,3]
b = [0.5,0.6,0.4,0.3,0.3,0.2]
c = [2,1.8,2.2,3.3,2.5,1.8]
sNames = ["a","b","c"]
xLabels = [2001,2002,2003,2004,2005,2006]

df = DataFrame(year=xLabels, a=a, b=b, c=c)

df |> stack |> @vlplot(:area, x=:year, y={:value, stack=:zero}, color="variable:n") |> save("some_plot")
