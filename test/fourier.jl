@testset "Utilities" begin 
    m,n = 46,58
    f = rand(ComplexF64,m,n)

    m2, n2 = m+24, n+36
    fpad = zeropad(f,(m2,n2))
    f2 = unzeropad(fpad,(m,n))
    @test f2 == f

    m, n = 13, 12
    f = rand(ComplexF64,m,n)
    m2, n2 = 15, 24
    fpad = zeropad(f,(m2,n2))
    f2 = unzeropad(fpad,(m,n))
    @test f2 == f  

    m2, n2 = m, 24
    fpad = zeropad(f,(m2,n2))
    f2 = unzeropad(fpad,(m,n))
    @test f2 == f  

    m2, n2 = 16, n
    fpad = zeropad(f,(m2,n2))
    f2 = unzeropad(fpad,(m,n))
    @test f2 == f  

    m2, n2 = m, n
    fpad = zeropad(f,(m2,n2))
    f2 = unzeropad(fpad,(m,n))
    @test f2 == f

end

@testset "Transforms"  begin
    Nx, Ny = 16, 32
    f = rand(Nx,Ny+1)
    fhat = physical2fourier(f,1)
    f2 = fourier2physical(fhat,1)
    @test norm(f-f2) < 1e-14

    Nx, Ny = 16, 32
    f = rand(Nx,Ny)
    fhat = physical2fourier(f,1)
    f2 = fourier2physical(fhat,1)

    @test norm(f-f2) < 1e-14

    fhat = physical2fourier(f,2)
    f2 = fourier2physical(fhat,2)

    @test norm(f-f2) < 1e-14

    fhat = physical2fourier(f)
    f2 = fourier2physical(fhat)

    @test norm(f-f2) < 1e-14


end

@testset "Derivatives" begin

    Nx, Ny = 128, 128
    x = range(0,2π,Nx+1)[1:Nx]
    y = range(0,2π,Ny+1)[1:Ny]

    ω = shearlayer(x,y)
    ω̂ = physical2fourier(ω)

    dωdx = fourierDx(ω̂)
    dωdx_ex = PseudoSpectral._dshearlayerdx(x,y)
    @test norm(dωdx - dωdx_ex) < 1e-3

    dωdy = fourierDy(ω̂)
    dωdy_ex = PseudoSpectral._dshearlayerdy(x,y)
    @test norm(dωdy - dωdy_ex) < 1e-3

end

@testset "Velocity from vorticity" begin
    Nx, Ny = 64, 64
    Δx, Δy = 2π/Nx, 2π/Ny
    x = range(0,2π,Nx+1)[1:Nx]
    y = range(0,2π,Ny+1)[1:Ny]
    ω = twovortex(x,y,distance=1.0)

    ω̂ = physical2fourier(ω)
    u, v = velocity(ω̂)
    ψ = streamfunction(ω̂)

    @test maximum(abs.(u)) ≈ 0.5267879828
    @test maximum(abs.(ψ)) ≈ 0.5286985881

    ψ̂ = PseudoSpectral._poisson_fourier_fourier(-ω̂)
    ψ̂[1,1] = 0

    ugradw_hat = PseudoSpectral._convective_derivative_fourier_fourier(ω̂,ψ̂)
    ugradw = fourier2physical(ugradw_hat)
    @test maximum(abs.(ugradw)) ≈ 4.5207103746

end