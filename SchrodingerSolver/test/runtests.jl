using Test
using SchrodingerSolver

@testset "SchrodingerSolver.jl Tests" begin
    
    x = collect(range(0.0, 1.0, length=1000))
    V = zeros(length(x)) # Potencjał równy 0 wewnątrz studni
    m = 1.0
    planck = 1.0
    
    E_analytic_1 = (pi^2 * planck^2) / (2.0 * m * 1.0^2)

    @testset "Test Metody FDM" begin
        energies_fdm, wavefuncs_fdm = solve_FDM(V, x, m, planck)
        
        @test isapprox(energies_fdm[1], E_analytic_1, rtol=0.01)
    end

    # ====================================================================
    # Odkomentować, jak będzie gotowy solve_Numerov
    # @testset "Test Metody Numerova" begin
    #     energies_num, wavefuncs_num = solve_Numerov(V, x, m, planck)
    #     
    #     @test isapprox(energies_num[1], E_analytic_1, rtol=0.001)
    # end

end