@testset "Utilities" begin 
    m,n = 46,58
    f = rand(ComplexF64,46,58)

    m2, n2 = m+24, n+36
    fpad = zeropad(f,(m2,n2))
    f2 = unzeropad(fpad,(m,n))
    @test maximum(abs.(f2-f)) == 0.0
end