function _convective_derivative_fourier_fourier(what::Array)

    m, n = size(what)
    mpad, npad = ceil(Int,3*m/2), ceil(Int,3*n/2)
    nfact = (mpad/n)*(npad/n)

    wpad = zeropad(what,(mpad,npad))

    dwdx_pad = nfact*ifft(_dfhatdx_fourier(wpad))
    dwdy_pad = nfact*ifft(_dfhatdy_fourier(wpad))

    uhat_pad, vhat_pad = _velocity_from_vorticity_fourier(wpad)
    upad = nfact*ifft(uhat_pad)
    vpad = nfact*ifft(vhat_pad)
    
    ugradw_hat = unzeropad(fft(upad.*dwdx_pad .+ vpad.*dwdy_pad)/nfact,(m,n))

end

"Compute convective derivative for Fourier Galerkin in x, Chebyshev collocation in y"
function _convective_derivative_fourier_chebyshev(what::Array)

    m, n = size(what)
    mpad, npad = ceil(Int,3*m/2), n
    nfact = (mpad/n)*(npad/n)

    wpad = zeropad(what,(mpad,npad))

    dwdx_pad = nfact*ifft(_dfhatdx_fourier(wpad))
    dwdy_pad = nfact*ifft(_dfhatdy_fourier(wpad)) # instead, use Chebyshev diff operator

    uhat_pad, vhat_pad = _velocity_from_vorticity_fourier(wpad) # need to solve this differently
    upad = nfact*ifft(uhat_pad)
    vpad = nfact*ifft(vhat_pad)
    
    ugradw_hat = unzeropad(fft(upad.*dwdx_pad .+ vpad.*dwdy_pad)/nfact,(m,n))

end

function _velocity_from_vorticity_fourier(what::Array)
    psihat = _streamfunction_fourier(what)
    psihat[1,1] = 0
    return _velocity_from_streamfunction_fourier(psihat)
end

function _velocity_from_streamfunction_fourier(psihat::Array)
    return _u_velocity_from_streamfunction_fourier(psihat), 
           _v_velocity_from_streamfunction_fourier(psihat)
end

_u_velocity_from_streamfunction_fourier(psihat::Array) = _dfhatdy_fourier(psihat)
_v_velocity_from_streamfunction_fourier(psihat::Array) = -_dfhatdx_fourier(psihat)


function _streamfunction_fourier(what::Array)
    m, n = size(what)
    _,_,ksq = _get_ksq_shifted(m,n)
    ksq[1,1] = 1
    return what./ksq
end