module PseudoSpectral

    using FFTW
    using UnPack
    using LinearAlgebra

    export zeropad, unzeropad, velocity, streamfunction, vorticity, physical2fourier,
            fourier2physical, VorticityFourierFourierNSCache, vorticity_ff_ns_step!, fourierDx, fourierDy,
            shearlayer, twovortex, randomfield

    export dct1, idct1, fchebt, ifchebt, chebdcoeffs, chebdiff, VorticityFourierChebyshevNSCache, vorticity_fc_ns_step!,
            randomfield_fourier_chebyshev

    abstract type AbstractPseudoSpectralCache end
    

    include("utilities.jl")
    include("diff.jl")
    include("chebyshev.jl")
    include("nonlinear.jl")
    include("fields.jl")
    include("navierstokes.jl")
    include("examples.jl")



end # module PseudoSpectral
