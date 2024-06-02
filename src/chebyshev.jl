dct1(f;kwargs...) = FFTW.r2r(f,FFTW.REDFT00;kwargs...)/(length(f)-1)
dct1(f,dims;kwargs...) = FFTW.r2r(f,FFTW.REDFT00,dims;kwargs...)/(length(f)-1)

idct1(f;kwargs...) = 0.5*FFTW.r2r(f,FFTW.REDFT00;kwargs...)
idct1(f,dims;kwargs...) = 0.5*FFTW.r2r(f,FFTW.REDFT00,dims;kwargs...)

"""
    fchebt(f::AbstractVector) -> AbstractVector

Given data `f` on the extrema Chebyshev points, evaluate
the Chebyshev transform.
"""
function fchebt(f::AbstractVector{T}) where {T<:Real}
    N = length(f)-1

    a = dct1(f)
    a[1] *= 0.5
    a[N+1] *= 0.5
    return a

end

"""
    ifchebt(f::AbstractVector) -> AbstractVector

Given data `f` on the extrema Chebyshev points, evaluate
the inverse Chebyshev transform.
"""
function ifchebt(a::AbstractVector{T}) where {T<:Real}
    N = length(a)-1

    f = copy(a)
    f[1] *= 2
    f[N+1] *= 2
    f .= idct1(f)
   
    return f
end

"""
    chebdcoeffs(a::AbstractVector) -> AbstractVector

Given the Chebshev representation of a set of data, return
the Chebshev representation of the derivative.
"""
function chebdcoeffs(a::AbstractVector{T}) where {T<:Real}

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