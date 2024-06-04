for fcn in (:energy_fourier,:two_point_correlation_fourier,:tke,:restress)
    @eval $fcn(what::Array) = $fcn(_velocity_from_vorticity_fourier(what)...)
end

for fcn in (:dissrate,)
    @eval $fcn(what::Array,a) = $fcn(_velocity_from_vorticity_fourier(what)...,a)
end

function dissrate(uhat::Array,vhat::Array,nu::Real)
    m, n = size(uhat)
    kx, ky, ksq = _get_ksq_shifted(m,n)
    Ek = energy_fourier(uhat,vhat)
    return 2*nu*sum(ksq.*Ek)
end

restress(uhat::Array,vhat::Array) = sum(two_point_correlation_fourier(uhat,vhat))

tke(uhat::Array,vhat::Array) = sum(energy_fourier(uhat,vhat))


function energy_fourier(uhat::Array,vhat::Array)
    R = two_point_correlation_fourier(uhat,vhat)
    return [real(0.5*Rk[1,1] + 0.5*Rk[2,2]) for Rk in R]
end

function two_point_correlation_fourier(uhat::Array,vhat::Array)
    Ruu = _two_point_correlation_fourier_uu(uhat,vhat)
    Rvv = _two_point_correlation_fourier_vv(uhat,vhat)
    Ruv = _two_point_correlation_fourier_uv(uhat,vhat)

    return [[Ruu_k Ruv_k; Ruv_k Rvv_k] for (Ruu_k,Rvv_k,Ruv_k) in zip(Ruu,Rvv,Ruv)]
end


function _two_point_correlation_fourier_uu(uhat::Array,vhat::Array)
    m,n = size(uhat)
    uhat.*conj(uhat)/(m*n)
end

function _two_point_correlation_fourier_vv(uhat::Array,vhat::Array)
    m,n = size(uhat)
    vhat.*conj(vhat)/(m*n)
end

function _two_point_correlation_fourier_uv(uhat::Array,vhat::Array)
    m,n = size(uhat)
    uhat.*conj(vhat)/(m*n)
end