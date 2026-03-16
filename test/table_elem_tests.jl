# test the table arithmetic implemented in julia

using Test
include("../src/minors_strategy.jl")

@testset verbose=true "TableElem" begin
    G = small_group(8, 3) # D8
    H = small_group(8, 4) # Q8
    Q = H
    l = 3
    vals, _ = working_values(G, Q, l)
    A, B = vals

    @testset "identity" begin
        @test A == A
        @test A != B
    end

    @testset "addition" begin
        @test A + B == B + A
        @test A - A == B - B
    end

    @testset "multiplication" begin
        @test A%B != B%A
        @test (A%B)%A == A%(B%A)
    end

    @testset "quasi-multiplication" begin
        @test A^2 == A * A
        @test A^3 == prod([A, A, A])
        @test A * B == A + B + A % B
        @test (A * B) * A == A * (B * A)
        @test A + A * B == A + (A * B)
        @test A * one(A) == A
        @test one(A) * A == A
        @test A * A^-1 == one(A)
        @test A^-1 * A == one(A)
        @test_throws AssertionError A^0
    end

end
