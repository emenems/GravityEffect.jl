using GravityEffect
using Base.Test
using DataFrames

# List of test files. Run the test from GravityEffect.jl folder
tests = ["simplebodies_test.jl",
		 "eopeffects_test.jl",
		 "atmacseffects_test.jl",
		 "eosteffects_test.jl",
		 "inpolygon_test.jl",
		 "layerresponse_test.jl",
		 "data2effect_test.jl"];
# Run all tests in the list
for i in tests
	include(i)
end
println("End test!")
