using NewtonsMethod
using Test

@testset "(x-1)^3 Analytical" begin
    f(x)  = (x-1)^3
    df(x) = 3*(x-1)^2
    out = newton_root(f, df)
    @test isapprox(out.value, 1, rtol=1e-6)
    @test out.iter < 1_000
end

@testset "(x-1)^3 AutoDiff" begin
    f(x) = (x-1)^3
    out = newton_root_AD(f)
    @test isapprox(out.value, 1, rtol=1e-6)
    @test out.iter < 1_000
end

@testset "(x-π)^2 Analytical" begin
    f(x)  = (x-π)^2
    df(x) = 2*(x-π)
    out = newton_root(f, df)
    @test isapprox(out.value, π, rtol=1e-6)
    @test out.iter < 1_000
end

@testset "(e^x)-2 AutoDiff" begin
    f(x) = exp(x)-2
    out = newton_root_AD(f)
    @test isapprox(out.value, log(2), rtol=1e-6)
    @test out.iter < 1_000
end

@testset "(x-3)^5 AutoDiff" begin
    f(x) = (x-3)^5
    out = newton_root_AD(f)
    @test isapprox(out.value, 3, rtol=1e-6)
    @test out.iter < 1_000
end

