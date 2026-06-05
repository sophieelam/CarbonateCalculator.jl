using Test
include("../../src/Calculator.jl")
using .Calculator

println("Running full test suite...")

@testset "Calculator.jl Suite" begin
    include("round_robin.jl")
end