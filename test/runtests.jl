using ElectronicStructure
using Test

@testset "Atom" begin
    a = Atom(:H, (1, 2, 3))
    b = Atom(:H, (1, 2, 3))
    @test a == b
    @test a === b
    @test_throws MethodError copy(a)
    @test_throws ErrorException Atom(:NonExistentElement, (1, 2, 3))
end

@testset "Geometry" begin
    a = Atom(:H, (0.0, 0.0, 0.0))
    b = Atom(:H, (0.0, 0.0, 0.74))
    g = Geometry(a, b)
    g1 = copy(g)
    @test g == g1
    @test length(g) == 2
    @test g[1] == a
    @test g[2] == b
end

