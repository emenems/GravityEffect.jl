using GravityEffect
using Base.Test
using DataFrames

# List of test files. Run the test from GravityEffect.jl folder
tests = ["gravityeffects_test.jl",
		 "eopeffects_test.jl",
		 "atmacseffects_test.jl"];
# Run all tests in the list
for i in tests
	include(i)
end
println("End test!")
