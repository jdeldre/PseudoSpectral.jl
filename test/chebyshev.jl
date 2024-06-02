@testset "DCT and DChebT" begin
    n = 128
    f = rand(n+1)
    f2 = idct1(dct1(f))
    @test norm(f-f2) < 1e-13

    f2 = dct1(idct1(f))
    @test norm(f-f2) < 1e-13

    f = rand(n+1)
    f2 = ifchebt(fchebt(f))
    @test norm(f-f2) < 1e-13

    f2 = fchebt(ifchebt(f))
    @test norm(f-f2) < 1e-13

end

@testset "Differentiation" begin
    n = 128
    x = cos.(π*(0:n)/n)
    c = rand()
    f = cos.(c*x)
    fp = fchebt(f);
    dfp = chebdcoeffs(fp)
    df = ifchebt(dfp)
    dfexact = -c*sin.(c*x)
    @test norm(df-dfexact) < 1e-11

    df = chebdiff(f)
    @test norm(df-dfexact) < 1e-11

    D1 = PseudoSpectral._construct_chebd1(n)
    @test norm(D1*f - df) < 1e-11

    d2f = chebdiff(df)
    D2 = PseudoSpectral._construct_chebd2(n,D1)
    @test norm(D2*f - d2f) < 1e-7


end

@testset "Solution of equations" begin

    Nx = 128
    x = cos.(π*(0:Nx)/Nx)
    uex = (1 .- x.^2)/2

    D1 = PseudoSpectral._construct_chebd1(Nx)
    D2 = PseudoSpectral._construct_chebd2(Nx,D1)

    A = copy(D2)
    A[1,:] .= 0
    A[1,1] = 1
    A[Nx+1,:] .= 0
    A[Nx+1,Nx+1] = 1

    b = zeros(Nx+1)
    b[2:Nx] .= -1

    u = A\b
    @test norm(u-uex) < 1e-9

    A[1,:] .= D1[1,:]
    uex = 2*(1 .- 0.25*(x .- 1).^2)
    u = A\b

    @test norm(u-uex) < 1e-9
end