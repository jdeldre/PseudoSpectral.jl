"""
    fourier2physical(fhat) -> f

Return the physical field associated with the Fourier-space representation `fhat`.
This function takes the inverse FFT, but also shifts the components in the array
to undo what is done in `physical_to_fourier`
"""
fourier2physical(fhat::Array{T}) where T<: Complex = real(fftshift(ifft(fhat)))

"""
    physical2fourier(f) -> fhat

Return the Fourier-space representation of the physical field `f`.
This function first shifts the array and then takes the Fourier transform.
"""
physical2fourier(f::Array{T}) where T<: Real = fft(fftshift(f))

"""
    velocity(what::Array{ComplexF64}) -> Array{Real}, Array{Real}

Return a tuple of the u and v components of the velocity field, given
the Fourier representation of the vorticity field `what`.
"""
function velocity(what::Array{T}) where T <: ComplexF64
    uhat, vhat = _velocity_from_vorticity_fourier(what)
    return fourier2physical(uhat), fourier2physical(vhat)
end

"""
    streamfunction(what::Array{ComplexF64}) -> Array{Real}

Return the streamfunction field, given
the Fourier representation of the vorticity field `what`.
"""
function streamfunction(what::Array{T}) where T <: ComplexF64
    psihat = _streamfunction_fourier(what)
    return fourier2physical(psihat)
end

"""
    vorticity(what::Array{ComplexF64}) -> Array{Real}

Return the vorticity field, given
the Fourier representation of the vorticity field `what`.
"""
function vorticity(what::Array{T}) where T <: ComplexF64
    return fourier2physical(what)
end

