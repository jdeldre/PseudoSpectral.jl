
"compute df/dx in Fourier space"
function _dfhatdx(fhat)
    Nx, _ = size(fhat) 
    kx = -Nx÷2:Nx÷2-1
    _dfhatdx_kshifted(fhat,fftshift(kx))
end

"compute df/dy in Fourier space"
function _dfhatdy(fhat)
    _, Ny = size(fhat)
    ky = -Ny÷2:Ny÷2-1
    _dfhatdy_kshifted(fhat,fftshift(ky))
end

"compute df/dx in Fourier space, using kx wavenumbers that have already been fft-shifted"
function _dfhatdx_kshifted(fhat,kx_shifted)
    return im*kx_shifted.*fhat
end

"compute df/dy in Fourier space, using ky wavenumbers that have already been fft-shifted"
function _dfhatdy_kshifted(fhat,ky_shifted)
    return im*transpose(ky_shifted.*transpose(fhat))
end