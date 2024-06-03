# Example initial conditions

"""
    shearlayer(x,y;[thickness=0.2, perturb = 0.02])

Return a smooth shear layer with a thickness `thickness` (0.2 by default)
and a sinusoidal perturbation amplitude `perturb` from flat. The wavelength
of the perturbation is 2π
"""
function shearlayer(x::AbstractVector,y::AbstractVector;thickness=0.2,perturb=0.02)
    return -1/thickness*sech.((y' .- π .- perturb*sin.(x))/thickness).^2
end

function _dshearlayerdy(x::AbstractVector,y::AbstractVector;thickness=0.2,perturb=0.02)
    return 2/thickness^2*sech.((y' .- π .- perturb*sin.(x))/thickness).^2 .* tanh.((y' .- π .- perturb*sin.(x))/thickness)
end

function _dshearlayerdx(x::AbstractVector,y::AbstractVector;thickness=0.2,perturb=0.02)
    return -2/thickness^2*perturb*cos.(x).*sech.((y' .- π .- perturb*sin.(x))/thickness).^2 .* tanh.((y' .- π .- perturb*sin.(x))/thickness)
end


function twovortex(x::AbstractVector,y::AbstractVector; radius = 0.2, distance = 1.0, circulation = 1.0)
    amp = circulation/(π*radius^2)
    return amp*exp.(-((x .- 0.5*distance .- π).^2 .+ (y' .- π).^2)/radius^2) .+ 
           amp*exp.(-((x .+ 0.5*distance .- π).^2 .+ (y' .- π).^2)/radius^2)
end


function randomfield(x::AbstractVector,y::AbstractVector; factor = 10.0)
    Nx, Ny = length(x), length(y)
    f = zeros(Nx,Ny)
    fhat = physical2fourier(f)
    Nxmax, Nymax = Nx÷4, Ny÷4
    for i in 2:Nxmax, j in 2:Nymax
        f0 = exp(im*2π*rand())
        fhat[i,j] = f0
        fhat[Nx-i+2,Ny-j+2] = conj(f0)
    end
    for i in Nx:-1:Nx-Nxmax+2, j in 2:Nymax
        f0 = exp(im*2π*rand())
        fhat[i,j] = f0
        fhat[Nx-i+2,Ny-j+2] = conj(f0)
    end
    
    fhat .*= factor*sqrt(Nx*Ny)
    
    f = fourier2physical(fhat)
end

function _randomfield_fourier!(fhat::Array{T};factor=1.0) where T<:Complex
    Nx, Ny = size(fhat)
    fill!(fhat,complex(0.0))
    Nxmax, Nymax = Nx÷4, Ny÷4
    for i in 2:Nxmax, j in 2:Nymax
        f0 = exp(im*2π*rand())
        fhat[i,j] = f0
        fhat[Nx-i+2,Ny-j+2] = conj(f0)
    end
    for i in Nx:-1:Nx-Nxmax+2, j in 2:Nymax
        f0 = exp(im*2π*rand())
        fhat[i,j] = f0
        fhat[Nx-i+2,Ny-j+2] = conj(f0)
    end
    return factor*fhat
end