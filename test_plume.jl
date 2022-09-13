using Test
include("plume.jl")

@testset "Second order central" begin
    f(x) = @. cos(x);
    fpp(x) = @. -cos(x);

    errs = zeros(2);
    i = 1;
    for Δ in [0.1, 0.001]
        x = [-Δ, 0.0, Δ];
        approx = second_order_central(f.(x), 2, Δ);
        errs[i] = abs(approx - fpp(0.0));
        i = i + 1;
    end
    
    # Estimate order of method and ensure its two
    p = log(errs[1]/errs[2])/log(100.);

    @test approx(p) == 2.0;

    # Ensure all errors are small
    @test all(errs .< 1e-3);

end

