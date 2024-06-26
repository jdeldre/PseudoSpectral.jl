# Example initial conditions

"""
    shearlayer(x,y;[thickness=0.2, perturb = 0.02])

Return a smooth shear layer with a thickness `thickness` (0.2 by default)
and a sinusoidal perturbation amplitude `perturb` from flat. The wavelength
of the perturbation is 2π
"""
function shearlayer(x::AbstractVector,y::AbstractVector;thickness=0.2,perturb=0.02)
    return _shearlayer.(x,y',Ref(thickness),Ref(perturb))
end
#function shearlayer(x::AbstractVector,y::AbstractVector;thickness=0.2,perturb=0.02)
#    return -1/thickness*sech.((y' .- π .- perturb*sin.(x))/thickness).^2
#end


function _dshearlayerdy(x::AbstractVector,y::AbstractVector;thickness=0.2,perturb=0.02)
    return _dshearlayerdy.(x,y',Ref(thickness),Ref(perturb))
end

function _dshearlayerdx(x::AbstractVector,y::AbstractVector;thickness=0.2,perturb=0.02)
    return _dshearlayerdx.(x,y',Ref(thickness),Ref(perturb))
end

function _shearlayer(x,y,thickness,perturb)
    return -1/thickness*sech((y-π-perturb*sin(x))/thickness)^2
end

function _dshearlayerdy(x,y,thickness,perturb)
    arg = (y - π - perturb*sin(x))/thickness
    return 2/thickness^2*sech(arg)^2*tanh(arg)
end

function _dshearlayerdx(x,y,thickness,perturb)
    arg = (y - π - perturb*sin(x))/thickness
    return -2/thickness^2*perturb*cos(x)*sech(arg)^2*tanh(arg)
end


## Two-vortex problem

"""
    twovortex(x,y;[radius=0.2, distance = 1.0, circulation = 1.0])

Return a set of two identical vortices of radius `radius` and circulation `circulation`,
separated by distance `distance` along the x axis and centered at the origin.
"""
function twovortex(x::AbstractVector,y::AbstractVector; radius = 0.2, distance = 1.0, circulation = 1.0)
    return _twovortex.(x,y',Ref(radius),Ref(distance),Ref(circulation))
end

function _twovortex(x,y,radius,distance,circulation)
    amp = circulation/(π*radius^2)
    return amp*exp(-((x - 0.5*distance - π)^2 + (y - π)^2)/radius^2) + 
           amp*exp(-((x + 0.5*distance - π)^2 + (y - π)^2)/radius^2)
end


"""
    randomfield(x,y;[factor=1.0])

Return a random field, uniformly distributed in wavenumbers from the lowest (non-zero)
wavenumber on the domains `x` and `y` up to half the Nyquist frequency (length(x)/4
and length(y)/4). The optional `factor` multiplies all non-zero Fourier components by
this factor.
"""
function randomfield(x::AbstractVector,y::AbstractVector; factor = 1.0)
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

function randomfield_fourier_chebyshev(x::AbstractVector,y::AbstractVector; factor = 10.0)
    Nx, Ny = length(x), length(y)-1
    f = zeros(Nx,Ny+1)
    fhat = physical2fourier(f,1)
    Nxmax = Nx÷4
    for i in 2:Nxmax, j in 2:Ny
        f0 = exp(im*2π*rand())
        fhat[i,j] = f0
        fhat[Nx-i+2,j] = conj(f0)
    end
    for i in Nx:-1:Nx-Nxmax+2, j in 2:Ny
        f0 = exp(im*2π*rand())
        fhat[i,j] = f0
        fhat[Nx-i+2,j] = conj(f0)
    end
    
    fhat .*= factor*sqrt(Nx*Ny)
    
    f = fourier2physical(fhat,1)
end