using ConditionalMean
using Base.Test

# use all(isapprox.(x,y)) elementwise isapprox of arrays

A=zeros(5,5)
A[:]=1:25
# all finite values

# write your own tests here
@test condmean(A[:],isfinite)[1] == mean(A)


f(x)=2x
c(x)=x>5
@test condmean(f,A,c)[2:5] == condmean(f,A,c,1)[2:5]
@test isnan(condmean(f,A,c)[1])
@test condmean(f,A,c,2) == mean(2*A,2)+5

# use a missing value
mv=-999.0
B=zeros(5,5)
B[:]=A[:]
B[3,:]=mv
@test condmean(B,mv,1) == mean(A,1)
B[:]=A[:]
B[2,:]=mv
@test condmean(B,mv,1) == mean(A,1)+0.25

# test using NaN and mapping the function through the mean
B[:]=A[:]
B[3,:]=NaN
# easy to predict for a linear operator
@test all(isapprox.(condmean(x->2*x,B,isfinite,1), 2*[3.0 8.0 13.0 18.0 23.0]))

# test the offsetting version
@show condmean(A,isfinite,2)
@show mean(A,2)
@test condmean(A,isfinite,2) == mean(A,2)

# test precision issues; first no missing values
B=Array{Float64}(A)
B[:]=A[:]+1e6
@show residual64 = (condmean(B[:],NaN)[1] - 1e6 - 13.0)
B=Array{Float32}(B)
@show residual32 = (condmean(B[:],NaN)[1] - 1e6 - 13.0)
if     !(residual64 ≈ Float64(0.0))
  warn("ConditionalMean: condmean imprecise for double precision.")
elseif !(residual32 ≈ Float32(0.0))
  warn("ConditionalMean: condmean imprecise for single precision.")
end
# require to be precise at double precision
@test residual64 ≈ Float64(0.0)

# test precision with missing value
B[3,:]=Float32(mv)
@test condmean(B[:],mv)[1] - 1e6 ≈ 13.0
@test all(isapprox.(condmean(B,mv,1) - 1e6, [3.0 8.0 13.0 18.0 23.0]))
# with NaN
B[3,:]=Float32(NaN)
@test all(isapprox.(condmean(B,isfinite,1) - 1e6, [3.0 8.0 13.0 18.0 23.0]))
