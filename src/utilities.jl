# This should be doable in place with view I think?

"""
    zeropad(f::Array,(mpad::Int,npad::Int))

Return an array in which the 2d Fourier transform data in `f` have been padded
with zeros to `(mpad,npad)`. 
"""
function zeropad(f::Array{T},padded::Tuple{Int,Int}) where {T}
    
    m, n = size(f)
    mmod, nmod = padded
    
    fshift = fftshift(f) # put in -N/2 -> N/2-1 form
    mshift = (mmod-m)÷2
    nshift = (nmod-n)÷2
    
    fmod = zeros(T,mmod,nmod)
    copyto!(fmod,CartesianIndices((1:m,1:n)),fshift,CartesianIndices(fshift))
    fmod .= circshift(fmod,(mshift,nshift))
    return fftshift(fmod)
end

"""
    unzeropad(fpad::Array,(m::Int,n::Int))

Return an array in which the padded 2d Fourier transform data in `fpad` have had the
extra padding removed, reduced to an array of size `(m,n)`. 
"""
function unzeropad(f::Array{T},unpadded::Tuple{Int,Int}) where {T}
    m, n = size(f)
    mmod, nmod = unpadded
    
    fshift = fftshift(f) # put in -N/2 -> N/2-1 form
    mshift = (mmod-m)÷2  # these should be negative
    nshift = (nmod-n)÷2
    
    fmod = zeros(T,mmod,nmod)
    fshift .= circshift(fshift,(mshift,nshift))
    copyto!(fmod,CartesianIndices(fmod),fshift,CartesianIndices((1:mmod,1:nmod)))
    return fftshift(fmod)
end


function _get_ksq_shifted(Nx::Int,Ny::Int)
    kx = fftshift(-Nx÷2:Nx÷2-1)
    ky = fftshift(-Ny÷2:Ny÷2-1)
    ksq = kx.^2 .+ transpose(ky).^2
    return kx, ky, ksq
end