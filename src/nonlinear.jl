function _convective_derivative_fourier_fourier(what::Array,psihat::Array)

    m, n = size(what)
    mpad, npad = ceil(Int,3*m/2), ceil(Int,3*n/2)
    nfact = (mpad/m)*(npad/n)

    wpad = zeropad(what,(mpad,npad))
    psipad = zeropad(psihat,(mpad,npad))

    dwdx_pad = nfact*ifft(_dfhatdx_fourier(wpad))
    dwdy_pad = nfact*ifft(_dfhatdy_fourier(wpad))

    uhat_pad = _u_velocity_from_streamfunction_fourier(psipad)
    vhat_pad = _v_velocity_from_streamfunction_fourier(psipad)

    upad = nfact*ifft(uhat_pad)
    vpad = nfact*ifft(vhat_pad)
    
    ugradw_hat = unzeropad(fft(upad.*dwdx_pad .+ vpad.*dwdy_pad)/nfact,(m,n))

end

"Compute convective derivative for Fourier Galerkin in x, Chebyshev collocation in y"
function _convective_derivative_fourier_chebyshev(what::Array,psihat::Array,D1::Matrix)

    m, n = size(what)
    mpad, npad = ceil(Int,3*m/2), n
    nfact = mpad/m

    wpad = zeropad(what,(mpad,npad))
    psipad = zeropad(psihat,(mpad,npad))

    dwdx_pad = nfact*ifft(_dfhatdx_fourier(wpad),1)
    dwdy_pad = nfact*ifft(_dfhatdy_chebyshev(wpad,D1),1)

    uhat_pad = _u_velocity_from_streamfunction_chebyshev(psipad,D1)
    vhat_pad = _v_velocity_from_streamfunction_fourier(psipad)

    upad = nfact*ifft(uhat_pad,1)
    vpad = nfact*ifft(vhat_pad,1)
    
    ugradw_hat = unzeropad(fft(upad.*dwdx_pad .+ vpad.*dwdy_pad,1)/nfact,(m,n))

end

function _velocity_from_vorticity_fourier(what::Array)
    psihat = _poisson_fourier_fourier(-what)
    psihat[1,1] = 0
    return _velocity_from_streamfunction_fourier(psihat)
end

function _velocity_from_streamfunction_fourier(psihat::Array)
    return _u_velocity_from_streamfunction_fourier(psihat), 
           _v_velocity_from_streamfunction_fourier(psihat)
end

_u_velocity_from_streamfunction_fourier(psihat::Array) =  _dfhatdy_fourier(psihat)
_v_velocity_from_streamfunction_fourier(psihat::Array) = -_dfhatdx_fourier(psihat)

_u_velocity_from_streamfunction_chebyshev(psihat::Array,D1::Matrix) = psihat*transpose(D1)

function _poisson_fourier_fourier(fhat::Array)
    m, n = size(fhat)
    _,_,ksq = _get_ksq_shifted(m,n)
    ksq[1,1] = 1
    return -fhat./ksq
end


function _poisson_fourier_chebyshev_dirichlet(what::Array,gplus_hat::Vector,gminus_hat::Vector,D2::Matrix)
    m = size(what,1)
    kx, kxsq = _get_ksq_shifted_1d(m)

    psihat = zero(what)

    k = 1
    psihat[k,:] = _helmholtz_chebyshev_dirichlet(-what[k,:],kxsq[k],gplus_hat[k],gminus_hat[k],D2)

    for k in 2:m÷2
        psihat[k,:] .= _helmholtz_chebyshev_dirichlet(-what[k,:],kxsq[k],gplus_hat[k],gminus_hat[k],D2)
        psihat[m-k+2,:] .= conj(psihat[k,:])
    end
    return psihat

end 

"Solves the equation u'' - σ*u = f on Chebyshev points, subject to u(1) = gplus, u(-1) = gminus"
function _helmholtz_chebyshev_dirichlet(f::Vector,sigma::Real,gplus::Complex,gminus::Complex,D2::Matrix)
    
    n = length(f)-1
    A = D2 - sigma*I
    A[1,:] .= 0
    A[1,1] = 1
    
    A[n+1,:] .= 0
    A[n+1,n+1] = 1
    
    b = f
    b[1] = gplus
    b[n+1] = gminus
    
    u = A\b
    return u
end