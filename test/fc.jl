@testset "Fourier and Chebyshev" begin
    Nx, Ny = 16, 32

    x = range(0,2π,length=Nx+1)[1:Nx]
    y = cos.(π*(0:Ny)/Ny)

    f = cos.(2*x).*(1 .- 2*transpose(y).^2)
    fhat = physical2fourier(f,1)
    fchat = fchebt(fhat,2)
    fhat2 = ifchebt(fchat,2)
    f2 = fourier2physical(fhat2,1)  
    
    @test norm(f-f2) < 1e-13

    f = cos.(2*transpose(x)).*(1 .- 2*y.^2)
    fhat = physical2fourier(f,2)
    fchat = fchebt(fhat,1)
    fhat2 = ifchebt(fchat,1)
    f2 = fourier2physical(fhat2,2)  
    
    @test norm(f-f2) < 1e-13

end