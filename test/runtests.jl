# test the main functions on the example D8 vs. Q8

# $ julia runtests.jl

using Test
include("../src/minors_strategy.jl")

include("./table_elem_tests.jl")

@testset verbose=true "MinorsStrategy" begin
    G = small_group(8, 3) # D8
    H = small_group(8, 4) # Q8
    Q = H
    l = 3
    F = free_group(:x, :y)
    x, y = gens(F)
    rels = [x^2, y^2]
    zero_rels, S = zero_relations(G, Q, l, rels)
    m = rank(G)
    n = rank(Q)
    mat = matrix(S, m, div(ngens(S), m), gens(S))[:, 1:n]
    d = det(mat)
    I = ideal(zero_rels)
    @test d^4 in I
    @test !(d^2 in I)
    @test radical_membership(d, I)

    @testset "distinguishing_quotients" begin
        C = small_group(8, 1) # C8
        @test describe(first(distinguishing_quotients(C, G))) == "C2 x C2"
        @test distinguishing_quotients(G, G) == PcGroup[]
    end

    @testset "splitting_layer" begin
        @test splitting_layer(G, H) == 2
        @test splitting_layer(G, H, deg=2) == 3
        @test_throws AssertionError splitting_layer(G, G)
    end

    @testset "all_minors_in_radical" begin
        @test !all_minors_in_radical(G, H)
        @test all_minors_in_radical(G, H, l=splitting_layer(G, H, deg=2))
        @test_throws AssertionError all_minors_in_radical(G, G)
    end

    @testset "practical_geometric_test" begin
        @test describe(practical_geometric_test(G, H, l=3)) == "Q8"
        @test practical_geometric_test(G, H, l=2) == nothing
    end

    @testset "dimension_modulo" begin
        @test dimension_modulo(G, 3) == 6
    end
end
