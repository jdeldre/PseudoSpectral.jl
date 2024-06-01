module PseudoSpectral

    using FFTW
    using UnPack

    export zeropad, unzeropad, velocity, streamfunction, vorticity, physical2fourier,
            fourier2physical, VorticityNSCache, vorticity_ns_step!

    abstract type AbstractPseudoSpectralCache end
    

    include("utilities.jl")
    include("diff.jl")
    include("nonlinear.jl")
    include("fields.jl")
    include("navierstokes.jl")



end # module PseudoSpectral
