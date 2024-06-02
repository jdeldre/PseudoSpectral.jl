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