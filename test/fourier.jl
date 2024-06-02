@testset "Utilities" begin 
    m,n = 46,58
    f = rand(ComplexF64,46,58)

    m2, n2 = m+24, n+36
    fpad = zeropad(f,(m2,n2))
    f2 = unzeropad(fpad,(m,n))
    @test maximum(abs.(f2-f)) == 0.0
end

@testset "Derivatives" begin

    Nx, Ny = 128, 128
    x = range(0,2π,Nx+1)[1:Nx]
    y = range(0,2π,Ny+1)[1:Ny]

    ω = shearlayer(x,y)
    ω̂ = physical2fourier(ω)

    dωdx = spectralDx(ω̂)
    dωdx_ex = PseudoSpectral._dshearlayerdx(x,y)
    @test norm(dωdx - dωdx_ex) < 1e-3

    dωdy = spectralDy(ω̂)
    dωdy_ex = PseudoSpectral._dshearlayerdy(x,y)
    @test norm(dωdy - dωdy_ex) < 1e-3

end

