struct VorticityFourierFourierNSCache{T,Nx,Ny} <: AbstractPseudoSpectralCache
    Re :: Real
    Δt :: Real
    ksq :: Array{T,2}
    ugradw_np1 :: Array{Complex{T},2}
    ugradw_n :: Array{Complex{T},2}
    ugradw_nm1 :: Array{Complex{T},2}

end

function VorticityFourierFourierNSCache(Nx::Int,Ny::Int,Re::Real,Δt::Real;dtype=Float64)
    kx, ky, ksq = _get_ksq_shifted(Nx,Ny)
    ugradw_np1 = zeros(Complex{dtype},Nx,Ny)
    ugradw_n = zeros(Complex{dtype},Nx,Ny)
    ugradw_nm1 = zeros(Complex{dtype},Nx,Ny)

    return VorticityFourierFourierNSCache{dtype,Nx,Ny}(Re,Δt,ksq,ugradw_np1,ugradw_n,ugradw_nm1)
end

function vorticity_ff_ns_step!(what::Array{T},params::VorticityFourierFourierNSCache{Nx,Ny}) where {T<:ComplexF64,Nx,Ny}
    @unpack Re, Δt, ksq, ugradw_nm1, ugradw_n = params

    # Solve for streamfunction
    psihat = _poisson_fourier_fourier(-what)
    psihat[1,1] = 0

    ugradw_nm1 .= ugradw_n
    ugradw_n .= _convective_derivative_fourier_fourier(what,psihat)

    # Implicit/Explicit Euler 
    #what .= (what .- Δt*ugradw_n)./(1 .+ Δt/Re*ksq)

    # Adams/Bashforth and Crank/Nicolson
    what .= ((1 .- 0.5*Δt/Re*ksq).*what .- 1.5*Δt*ugradw_n .+ 0.5*Δt*ugradw_nm1)./(1 .+ 0.5*Δt/Re*ksq)


end

#= Fourier-Chebyshev solver =#

struct VorticityFourierChebyshevNSCache{T,Nx,Ny} <: AbstractPseudoSpectralCache
    Re :: Real
    Δt :: Real
    ksq :: Array{T,1}
    ugradw_np1 :: Array{Complex{T},2}
    ugradw_n :: Array{Complex{T},2}
    ugradw_nm1 :: Array{Complex{T},2}
    fhat :: Array{Complex{T},2}
    gplus_hat :: Array{Complex{T},1}
    gminus_hat :: Array{Complex{T},1}
    hplus_hat :: Array{Complex{T},1}
    hminus_hat :: Array{Complex{T},1}
    D1 :: Matrix{T}
    D2 :: Matrix{T}
end

function VorticityFourierChebyshevNSCache(Nx::Int,Ny::Int,Re::Real,Δt::Real,gplus::Vector,gminus::Vector,hplus::Vector,hminus::Vector;dtype=Float64)
    kx, ksq = _get_ksq_shifted_1d(Nx)
    ugradw_np1 = zeros(Complex{dtype},Nx,Ny+1)
    ugradw_n = zeros(Complex{dtype},Nx,Ny+1)
    ugradw_nm1 = zeros(Complex{dtype},Nx,Ny+1)
    fhat = zeros(Complex{dtype},Nx,Ny+1)

    D1 = _construct_chebd1(Ny)
    D2 = _construct_chebd2(Ny,D1)

    gplus_hat = physical2fourier(gplus)
    gminus_hat = physical2fourier(gminus)
    hplus_hat = physical2fourier(hplus)
    hminus_hat = physical2fourier(hminus)

    return VorticityFourierChebyshevNSCache{dtype,Nx,Ny}(Re,Δt,ksq,ugradw_np1,ugradw_n,ugradw_nm1,fhat,
                                                        gplus_hat, gminus_hat, hplus_hat, hminus_hat, D1,D2)
end
 
function VorticityFourierChebyshevNSCache(Nx::Int,Ny::Int,Re::Real,Δt::Real,gplus::Real,gminus::Real,hplus::Real,hminus::Real;dtype=Float64)
    VorticityFourierChebyshevNSCache(Nx,Ny,Re,Δt,gplus*ones(Nx),gminus*ones(Nx),hplus*ones(Nx),hminus*ones(Nx);dtype=dtype)
end

function vorticity_fc_ns_step!(what,psihat,params::VorticityFourierChebyshevNSCache)
    @unpack Re, Δt, ksq, ugradw_np1, ugradw_n, ugradw_nm1, fhat, gplus_hat, gminus_hat, hplus_hat, hminus_hat, D1, D2 = params

    # call convective derivative to set ugradw_n and set ugradw_nm1 to previous one
    ugradw_nm1 .= ugradw_n
    ugradw_n .= _convective_derivative_fourier_chebyshev(what,psihat,D1)

    # This is to use AB for the non-linear term, but in a CN setting
    ugradw_np1 .= 2*ugradw_n .- ugradw_nm1 # This should be just ugradw_np1 .= ugradw_n at first step

    # build fhat
    #ϵ, θ = 1, 0.5 # Crank-Nicolson
    kfact = ksq .- 2*Re/Δt
    fhat .= kfact.*what - what*transpose(D2) + Re*ugradw_np1 + Re*ugradw_n    
    
    sigma0 = 2*Re/Δt

    # advance with _vorticity_and_streamfunction_chebyshev
    what_np1, psihat_np1 = _vorticity_and_streamfunction_chebyshev(fhat,sigma0,ksq,gplus_hat,gminus_hat,hplus_hat,hminus_hat,D1,D2)
    what .= what_np1
    psihat .= psihat_np1

end

function _vorticity_and_streamfunction_chebyshev(fhat::Array,sigma0::Number,ksq::AbstractVector,gplus_hat::Vector,gminus_hat::Vector,hplus_hat::Vector,hminus_hat::Vector,D1::Matrix,D2::Matrix)

    m = size(fhat,1)
    n = size(fhat,2) - 1

    what = zero(fhat)
    psihat = zero(fhat)

    # Run through all of the wavenumbers
    k = 1
    what_k, psihat_k = _vorticity_and_streamfunction_chebyshev_k(fhat[k,:],ksq[k]+sigma0,ksq[k],gplus_hat[k],gminus_hat[k],hplus_hat[k],hminus_hat[k],D1,D2)
    what[k,:] .= what_k
    psihat[k,:] .= psihat_k

    for k in 2:m÷2
        what_k, psihat_k = _vorticity_and_streamfunction_chebyshev_k(fhat[k,:],ksq[k]+sigma0,ksq[k],gplus_hat[k],gminus_hat[k],hplus_hat[k],hminus_hat[k],D1,D2)
        what[k,:] .= what_k
        psihat[k,:] .= psihat_k

        what[m-k+2,:] .= conj(what[k,:])
        psihat[m-k+2,:] .= conj(psihat[k,:])

    end
    return what, psihat

end

function _vorticity_and_streamfunction_chebyshev_k(f::Vector,sigma_k::Number,ksq::Number,gplus::Number,gminus::Number,hplus::Number,hminus::Number,D1::Matrix,D2::Matrix)
    n = length(f) - 1

    # P̃ problem
    wplus, wminus = 0,0
    w̃, ψ̃ = _general_vorticity_and_streamfunction_chebyshev(f,sigma_k,ksq,wplus,wminus,gplus,gminus,D2)

    fzero = zero(f)

    # P1 problem
    wplus, wminus = 1,0
    w1, ψ1 = _general_vorticity_and_streamfunction_chebyshev(fzero,sigma_k,ksq,wplus,wminus,0,0,D2)

    # P2 problem
    wplus, wminus = 0,1
    w2, ψ2 = _general_vorticity_and_streamfunction_chebyshev(fzero,sigma_k,ksq,wplus,wminus,0,0,D2)

    dψ1 = D1*ψ1
    dψ2 = D1*ψ2
    dψ̃ = D1*ψ̃

    M = [dψ1[1] dψ2[1]; dψ1[n+1] dψ2[n+1]]
    E = [hplus-dψ̃[1], hminus-dψ̃[n+1]]
    ξ = M\E

    return w̃ + ξ[1]*w1 + ξ[2]*w2, ψ̃ + ξ[1]*ψ1 + ξ[2]*ψ2
end


function _general_vorticity_and_streamfunction_chebyshev(f::Vector,sigma_k::Number,ksq::Number,
                                                wplus::Number,wminus::Number,
                                                psiplus::Number,psiminus::Number,D2::Matrix)
    w = _helmholtz_chebyshev_dirichlet(f,sigma_k,wplus,wminus,D2)
    psi = _helmholtz_chebyshev_dirichlet(-w,ksq,psiplus,psiminus,D2)
    return w, psi
end