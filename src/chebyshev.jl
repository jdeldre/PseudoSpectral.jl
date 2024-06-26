

"Normalized type-I discrete cosine transforms (on the Chebyshev extrema points)"
dct1(f;kwargs...) = FFTW.r2r(f,FFTW.REDFT00;kwargs...)/(length(f)-1)
dct1(f,dims;kwargs...) = FFTW.r2r(f,FFTW.REDFT00,dims;kwargs...)/(size(f,dims)-1)

"Normalized type-I inverse discrete cosine transforms (on the Chebyshev extrema points)"
idct1(f;kwargs...) = 0.5*FFTW.r2r(f,FFTW.REDFT00;kwargs...)
idct1(f,dims;kwargs...) = 0.5*FFTW.r2r(f,FFTW.REDFT00,dims;kwargs...)

"""
    fchebt(f::AbstractVector[,dim::Integer]) -> AbstractVector

Given data `f` on the extrema Chebyshev points, evaluate
the Chebyshev transform. The optional argument `dim` specifies
the dimension along which to take the transform.
"""
function fchebt(f::AbstractVector{T}) where {T<:Number}
    N = length(f) - 1
    a = dct1(f)
    _adjust_fchebt!(a,N,Val(0))
end

function fchebt(f::AbstractArray{T},dim) where {T<:Number}
    N = size(f,dim) - 1
    a = dct1(f,dim)
    _adjust_fchebt!(a,N,Val(dim))
end

function _adjust_fchebt!(a,N,::Val{0})
    a[1] *= 0.5
    a[N+1] *= 0.5
    return a
end

function _adjust_fchebt!(a,N,::Val{1})
    a[1,:] *= 0.5
    a[N+1,:] *= 0.5
    return a
end

function _adjust_fchebt!(a,N,::Val{2})
    a[:,1] *= 0.5
    a[:,N+1] *= 0.5
    return a
end

"""
    ifchebt(f::AbstractVector[,dim]) -> AbstractVector

Given data `f` on the extrema Chebyshev points, evaluate
the inverse Chebyshev transform. The optional argument `dim` specifies
the dimension along which to take the inverse transform.
"""
function ifchebt(a::AbstractVector{T}) where {T<:Number}
    N = length(a)-1
    f = copy(a)
    _adjust_ifchebt!(f,N,Val(0))
    idct1(f)
end

function ifchebt(a::AbstractArray{T},dim) where {T<:Number}
    N = size(a,dim) - 1
    f = copy(a)
    _adjust_ifchebt!(f,N,Val(dim))
    idct1(f,dim)
end

function _adjust_ifchebt!(f,N,::Val{0})
    f[1] *= 2
    f[N+1] *= 2
    return f
end

function _adjust_ifchebt!(f,N,::Val{1})
    f[1,:] *= 2
    f[N+1,:] *= 2
    return f
end

function _adjust_ifchebt!(f,N,::Val{2})
    f[:,1] *= 2
    f[:,N+1] *= 2
    return f
end

"""
    chebdcoeffs(a::AbstractVector) -> AbstractVector

Given the Chebshev representation of a set of data, return
the Chebshev representation of the derivative.
"""
function chebdcoeffs(a::AbstractVector{T}) where {T<:Number}

    N = length(a)-1
    da = zero(a)
    da[N] = float(2N)*a[N+1]
    for n in N-1:-1:1
        da[n] = da[n+2] + float(2n)*a[n+1]
    end
    da[1] *= 0.5
    return da

end

"""
    chebdiff(f::AbstractVector) -> AbstractVector

Given a set of data on the Chebshev extrema points, return
the derivative of the data, computed via a Chebyshev transform
"""
function chebdiff(f::AbstractVector)
    a = fchebt(f)
    _zero_small_values!(a)
    da = chebdcoeffs(a)
    return ifchebt(da)
end

"""
    chebint(f[,dim]) -> Real

Given a set of data on the Chebshev extrema points, return
the integral (quadrature) of the data
"""
function chebint(f::AbstractVector)
    n = length(f) - 1
    wts = sin.(π*(0:n)/n)
    c = ones(n+1)
    c[1] = c[n+1] = 0.5
    return sum(c.*f.*wts)*π/n

end

function chebint(f::AbstractArray,dim::Int)
    n = size(f,dim) - 1
    wts = sin.(π*(0:n)/n)
    c = ones(n+1)
    c[1] = c[n+1] = 0.5
    _chebint(c,wts,f,n,Val(dim))
end

function _chebint(c,wts,f,n,::Val{1})
    sum(c.*wts.*f,dims=1)*π/n
end

function _chebint(c,wts,f,n,::Val{2})
    vec(sum(transpose(c.*wts.*transpose(f)),dims=2)*π/n)
end

function _zero_small_values!(a)
    a[abs.(a).<eps()] .= 0.0
end


function _construct_chebd1(N::Int)

    c = ones(Float64,N+1)
    c[1] = c[N+1] = 2
    D1 = zeros(N+1,N+1)
    for j in 1:N+1, k in 1:N+1
        j == k && continue
        kj = j+k
        xjxk = 2*sin(0.5π*(kj-2)/N)*sin(0.5π*(k-j)/N)
        D1[j,k] = c[j]/c[k]*(-1)^kj/xjxk
        D1[j,j] -= D1[j,k]
    end
    return D1
end

function _construct_chebd2(N,D1)
    D2 = D1*D1

    for j in 1:N+1
        D2[j,j] = 0
        for k in 1:N+1
            j == k && continue
            D2[j,j] -= D2[j,k]
        end
    end
    return D2
end

function _set_dirichlet_at_yplus1!(A)
end