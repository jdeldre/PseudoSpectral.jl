struct VorticityNSCache{T,Nx,Ny} <: AbstractPseudoSpectralCache
    Re :: Real
    Δt :: Real
    ksq :: Array{T,2}
    ugradw :: Array{Complex{T},2}
end

function VorticityNSCache(Nx::Int,Ny::Int,Re::Real,Δt::Real;dtype=Float64)
    kx, ky, ksq = _get_ksq_shifted(Nx,Ny)
    ugradw = zeros(Complex{dtype},Nx,Ny)
    return VorticityNSCache{dtype,Nx,Ny}(Re,Δt,ksq,ugradw)
end

function vorticity_ns_step!(what::Array{T},params::VorticityNSCache{Nx,Ny}) where {T<:ComplexF64,Nx,Ny}
    @unpack Re, Δt, ksq, ugradw = params
    ugradw .= _convective_derivative_fourier_fourier(what)
    what .= (what .- Δt*ugradw)./(1 .+ Δt/Re*ksq)

end