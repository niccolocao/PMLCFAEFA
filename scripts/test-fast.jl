using Test, Random, LinearAlgebra, Statistics

include(joinpath(@__DIR__, "..", "julia", "PMLRBridge", "src", "PMLRBridge.jl"))
const PML = PMLRBridge.PMLCFAEFA

function fixture(n::Int, p::Int, q::Int, K::Int, c::Int, seed::Int)
    Random.seed!(seed)
    X = rand(1:c, n, p)
    npairs = binomial(p, 2)
    pairs = Matrix{Int}(undef, 2, npairs)
    idx = 0
    for i in 1:p-1, j in i+1:p
        idx += 1
        pairs[1,idx] = i
        pairs[2,idx] = j
    end
    L = zeros(Float64, p, q)
    for i in 1:p
        L[i, mod1(i,q)] = 1.0
    end
    indices = [i + (j-1)*p for j in 1:K for i in j:p]
    tau1 = repeat(reshape(collect(range(-1.5, 1.5; length=c-1)), 1, c-1), p, 1)
    tau2 = repeat(reshape(collect(range(0.3, 1.4; length=c-3)), 1, c-3), p, 1)
    x = vcat(vec(tau1), fill(0.6, p), fill(0.1, q*(q-1)÷2),
             vec(tau2), fill(0.35, length(indices)),
             fill(1.2, p), fill(0.15, K))
    T = zeros(Float64, n, 2, npairs)
    for pair in 1:npairs, r in 1:n
        t = 0.1 + 0.8*rand()
        T[r,1,pair] = t
        T[r,2,pair] = 1-t
    end
    return (; X, T, pairs, L, indices, x, n, p, q, K, c)
end

function verify_case(z; gradient=true)
    W = PML.pml_pair_category_weights(z.X,z.T,z.pairs,z.c)
    @test size(W) == (z.c^2,2,binomial(z.p,2))
    @test all(isapprox(sum(W[:,k,j]),sum(z.T[:,k,j]); atol=1e-11,rtol=1e-12)
              for k in 1:2, j in axes(z.pairs,2))
    call_obj(w) = PML.fglm!(1.0,nothing,z.x,z.X,z.T,z.p,z.q,z.K,z.c,z.n,z.pairs,z.L,z.indices,w)
    slow_f = call_obj(nothing)
    fast_f = call_obj(W)
    @test isfinite(slow_f) && isfinite(fast_f)
    @test isapprox(slow_f,fast_f; atol=1e-7,rtol=1e-10)
    println("p=",z.p," n=",z.n," q=",z.q," K=",z.K,
            " objective baseline=",slow_f," fast=",fast_f,
            " abs_difference=",abs(slow_f-fast_f))
    if gradient
        slow_g = zeros(length(z.x))
        fast_g = zeros(length(z.x))
        PML.fglm!(nothing,slow_g,z.x,z.X,z.T,z.p,z.q,z.K,z.c,z.n,z.pairs,z.L,z.indices,nothing)
        PML.fglm!(nothing,fast_g,z.x,z.X,z.T,z.p,z.q,z.K,z.c,z.n,z.pairs,z.L,z.indices,W)
        @test all(isfinite,slow_g)
        @test all(isfinite,fast_g)
        max_diff = maximum(abs.(slow_g .- fast_g))
        @test all(isapprox.(slow_g,fast_g; atol=1e-7,rtol=1e-9))
        println("gradient max_abs_difference=",max_diff)
    end
    return W
end

@testset "ordinal pair aggregation vs original analytical M-step" begin
    z1 = fixture(70,8,2,1,5,32)
    verify_case(z1)
    z2 = fixture(65,12,3,2,5,33)
    verify_case(z2)
    z3 = fixture(55,42,5,1,5,34)
    verify_case(z3)
end

println("All numerical equivalence tests passed.")
