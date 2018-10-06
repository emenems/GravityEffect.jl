using GravityEffect
using Test
using DataFrames
using Dates
import DelimitedFiles
import ResampleAndFit

# List of test files. Run the test from GravityEffect.jl folder
tests = ["simplebodies.jl",
		 "eopeffects.jl",
		 "atmacseffects.jl",
		 "eosteffects.jl",
		 "inpolygon.jl",
		 "layerresponse.jl",
		 "data2effect.jl"];
# Run all tests in the list
for i in tests
	include(i)
end
