using ConditionalMean
using Base.Test

A=zeros(5,5)
A[:]=1:25
# all finite values

# write your own tests here
@test condmean(A,isfinite) == mean(A)
@test condmean(A,true,2) == mean(A,2)

f=2x
c=x>5
@test condmean(f,A,c) == condmean(f,A,c,1)
@test condmean(f,A,c,2) == mean(2*A,2)

# use a missing value
mv=-999.0
B=A
B[3,:]=mv
@test condmean(B,mv,1) == mean(A,1)

# test precision issues
B=map(Float32,A+1e6)
@test condmean(B[:],true) - 1e6 == 13.0
B[3,:]=mv
@test (condmean(B,true) .> 799803) .== [false, true, true, true, true].'
@test isequal(condmean(B[:],mv)[1] - 1e6, 13.0)
@test isequal(condmean(B,mv,1) - 1e6, [3.0 8.0 13.0 18.0 23.0])
B[3,:]=NaN
@test isequal(condmean(B,isfinite,1) - 1e6, [3.0 8.0 13.0 18.0 23.0])

# test mapping the function
# easy to predict for a linear operator
@test isequal(condmean(x->2*x,B,isfinite,1) - 2e6, 2*[3.0 8.0 13.0 18.0 23.0])
