module ConditionalMean

using Base: indices1, tail
export condmean, condvar, condstd
export nanmean, nanstd

"""
`M = condmean(A, cond, region)` calculates the mean value of A along the dimensions listed in `region`, ignoring values when cond is not true, or values equal to cond if cond is a number.
`M = condmean(f, A, cond, region)` calculates the mean of f(A).
"""
#condmean{T<:Real}(A::AbstractArray{T}, cond::Function, region=1) = _condmean(A, T, region)
# vectorizing mean(A[:]) is trivial, so keep default region definition of 1
condmean{T<:Number}(             A::AbstractArray{T}, cond::Function=(x->true), region=1; Km1::Bool=false) = _condmean(identity, A, T, cond      , region, Km1)
condmean{T<:Number}(             A::AbstractArray{T}, cond::Number            , region=1; Km1::Bool=false) = _condmean(identity, A, T, x->x!=cond, region, Km1)
condmean{T<:Number}(f::Function, A::AbstractArray{T}, cond::Function=(x->true), region=1; Km1::Bool=false) = _condmean(f       , A, T, cond      , region, Km1)
condmean{T<:Number}(f::Function, A::AbstractArray{T}, cond::Number            , region=1; Km1::Bool=false) = _condmean(f       , A, T, x->x!=cond, region, Km1)

# redirect _condmean to _condmeanoffset for Array{Float32} to improve precision
#_condmean{T<:Float32}(f::Function, A::AbstractArray, ::Type{T}, cond, region, Km1::Bool=false) = _condmeanoffset(f, A, T, cond, region, Km1)

# # general condition cond generalizes to all Number types, but behavior maybe unpredictable for some
# function _condmean{T<:Number}(f::Function, A::AbstractArray, ::Type{T}, cond, region, Km1::Bool=false)
#     sz = Base.reduced_dims(A, region)
#     K = zeros(Int, sz)
#     S = zeros(eltype(A), sz)
#     condsum!(S, K, A, f, cond)
#     if Km1
#        S./(K-1)
#     else
#        S./K
#     end
# end

# add offset to improve precision for floating point numerics, behavior unpredictable for other types
function _condmeanoffset{T<:Number}(f::Function, A::AbstractArray, ::Type{T}, cond, region, Km1::Bool=false)
    sz = Base.reduced_indices(A, region)
    K = zeros(Int, sz)
    S = zeros(eltype(A), sz)
    P = zeros(eltype(A), sz)
    condsumoffset!(S, P, K, A, f, cond)
    if Km1
        (S + P.*K) ./ (K-1)
    else
        S./K + P
    end
end

# always redirect _condmean to _condmeanoffset to improve precision
_condmean=_condmeanoffset

using Base: check_reducedims, reducedim1, safe_tail
using Base.Broadcast: newindex

"""
    condsum!(S, K, A, f, cond)

Compute the sum `S` and number of contributing values `K` for
reductions of the array `A` over dimensions. `S` and `K` must have
identical indices, and either match `A` or have singleton-dimensions for
the dimensions that are being summed over. Only values x with cond(x)==true
are included in the tallies of `S` and `K`.

Note that the mean is just S./K.
"""
function condsum!{T,N}(S, K, A::AbstractArray{T,N}, f, cond)
    check_reducedims(S, A)
    isempty(A) && return S, K
    indices(S) == indices(K) || throw(DimensionMismatch("S and K must have identical indices"))

    indsAt, indsSt = safe_tail(indices(A)), safe_tail(indices(S))
    keep, Idefault = _newindexer(indsAt, indsSt)
    if reducedim1(S, A)
        # keep the accumulators as a local variable when reducing along the first dimension
        i1 = first(indices1(S))
        @inbounds for IA in CartesianRange(indsAt)
            IS = newindex(IA, keep, Idefault)
            s, k = S[i1,IS], K[i1,IS]
            for i in indices(A, 1)
                tmp = A[i, IA]
                if cond(tmp)
                    s += f(tmp)
                    k += 1
                end
            end
            S[i1,IS], K[i1,IS] = s, k
        end
    else
        @inbounds for IA in CartesianRange(indsAt)
            IS = newindex(IA, keep, Idefault)
            for i in indices(A, 1)
                tmp = A[i, IA]
                if cond(tmp)
                    S[i, IS] += f(tmp)
                    K[i, IS] += 1
                end
            end
        end
    end
    S, K
end
if VERSION < v"0.6.0-dev.693"
    _newindexer(shape, inds) = Base.Broadcast.newindexer(shape, inds)
else
    _newindexer(shape, inds) = Base.Broadcast.shapeindexer(shape, inds)
end

# add offset to improve precision for floating point numerics, behavior unpredictable for other types
"""
    condsumoffset!(S, P, K, A, f, cond)

Compute the adjusted sum `S` and number of contributing values `K` for
reductions of `f` mapped on array `A` over dimensions.
`S`, `P`, and `K` must have identical indices, and either match `A`
or have singleton-dimensions for the dimensions that are being summed
over. Only values of A with cond(A[i])==true are included in the tallies
of `S` and `K`. `P` is the first good value in each A set of `A` values
to be summed.

By subtracting off the first good value `P` accumulated to each sum,
condsumoffset! improves precision in case the values of A have a nearly
constant offset. The improved precision is beneficial if `K`×`P` is so large
that interesting variations of A are truncated in the sum S. However, precision
of this implementation will be worse than condsum! if the first good value is an
outlier farther from zero than the mean of the values of `A` summed.

The value mean computed in `condmeanprecise` is S./K+P.

See also `condsum!` and `condmeanoffset`.
"""
function condsumoffset!{T,N}(S, P, K, A::AbstractArray{T,N}, f, cond)
  #loops over dimensions of A _then_ determines dimensions of reduced S
    check_reducedims(S, A)
    isempty(A) && return S, K
    indices(S) == indices(K) || throw(DimensionMismatch("S and K must have identical indices"))

    indsAt, indsSt = safe_tail(indices(A)), safe_tail(indices(S))
    keep, Idefault = _newindexer(indsAt, indsSt)
    if reducedim1(S, A)
        # keep the accumulators as a local variable when reducing along the first dimension
        i1 = first(indices1(S))
        @inbounds for IA in CartesianRange(indsAt)
            IS = newindex(IA, keep, Idefault)
            s, p, k = S[i1,IS], P[i1,IS], K[i1,IS]
            #isfirst = true
            for i in indices(A, 1)
                tmp = A[i, IA]
                if cond(tmp)
                    if k==0 #isfirst
                        p = f(tmp)
                        #isfirst = false
                    else
                        s += f(tmp)-p
                    end
                    k += 1
                end
            end
            S[i1,IS], P[i1,IS], K[i1,IS] = s, p, k
        end
    else
        @inbounds for IA in CartesianRange(indsAt)
            IS = newindex(IA, keep, Idefault)
            #isfirst = true
            for i in indices(A, 1)
                tmp = A[i, IA]
                if cond(tmp)
                    if K[i, IS]==0 #isfirst
                        P[i, IS] = f(tmp)
                        #isfirst = false
                    else
                        S[i, IS] += f(tmp)-P[i, IS]
                    end
                    K[i, IS] += 1
                end
            end
        end
    end
    S, P, K
end
#=
if testing for k==0 is slow, then an alternative that works in many cases for
nearly constant offsets across the whole of A would be to just subtract off the
first good value in all of A. This would be good for surface pressure, but not
good for a significant gradient of A, like atmospheric pressure from the surface
to the stratosphere, where p varies by an order of magnitude.
Another solution, that I don't think will be faster, is to make a boolean array
`isfirst` like `K` but private to condsumoffset that stores whether the offset
has been recorded in `P`.
=#

# NOT TESTED
# precise, but requires two passes of condmean through A
# instead could write a version of condsum that accumulates multiple outputs of f, e.g. f(x)=[x, x*x]
"""
`condvar(A, cond, region; m=knownmean, corrected=true)` gives the unbiased
estimate of the variance of A (normalized by N-1) along dims in region.
`corrected=false` gives the variance normalized by N.
"""
condvar(A,cond,region; m=condmean(A,cond,region), corrected::Bool=true) =       condmean(x->abs(x)*abs(x), broadcast(-,A,m), cond, region, Km1=corrected)
condstd(A,cond,region; m=condmean(A,cond,region), corrected::Bool=true) = sqrt.(condmean(x->abs(x)*abs(x), broadcast(-,A,m), cond, region, Km1=corrected))
# corrected=true normalizes variance by K-1

"nanmean(X,r) mean of X long dimension(s) r, ignoring NaNs."
nanmean(X,r)=condmean(X,!(isnan),r)
"nanstd(X,r) standard deviation of X along dimension(s) r, ignoring NaNs."
nanstd(X,r)=sqrt.(condmean(x->abs(x)*abs(x),X,!(isnan),r))
end # module

