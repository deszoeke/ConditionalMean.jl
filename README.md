# ConditionalMean

Provides general methods for summing and averaging values that meet a specific condition, over various dimensions of Julia arrays. A useful example is to average finite values and ignore missing values.

```condmean(A, cond, region)```
returns the mean value of A along the dimensions listed in `region`, ignoring values when `cond` is `false`, or ignoring values equal to `cond` if `cond` is a number.

An inner function may be supplied. For example,
`condmean(f, A, cond, region)`
calculates the mean of f(A).

The method of averaging over `region` is copied from the [JuliaImages/Images](https://github.com/JuliaImages/Images.jl) package. For multiple dimensions in `region`, it averages *all* data over the specified dimensions, regardless of the order of looping over dimensions. I.e.
```condmean(X,NaN,[1,2])``` averages all non-NaN data with equal weight along dimensions 1 and 2
while
```condmean(condmean(X,NaN,1),NaN,2)``` averages along dimension 1, then averages that result along dimension 2.

## Missing values

Julia inherits different standards and functions for handling missing values from many languages. Missing values may be indicated by
1. `missing`
2. DataFrames `NA` type
3. IEEE `NaN`
4. `Int64(-999)` or some other value
5. any negative value

The desired behavior for handing missing values is almost always either
1. The result of any calculation involving any missing value(s) is a missing value.
or
2. Perform the calculation ignoring missing values. MATLAB's `nanmean()` does this.

You might want to read data formatted with one missing value convention, process it with methods using another, and maybe even plot data with methods using a third convention. `ConditionalMean` gives you tools to ignore missing values minimizing or eliminating the need to convert custom types at each step.

### Examples
```
nanmean(X,r) = condmean(X, !(isnan), r)
```
averages all values that are not `NaN` over the region `r`. These examples might suit better, depending on how the data is formatted.
```
condmean(X, isfinite, r)  # the default
condmean(X, x -> x!=-999) # use a generic function to filter out missing values equal to -999
condmean(X, -999)         # another method: pass the missing value as a number to the condition argument
```

You can also call `condmean` on a function of the array. For example,
```
f(x) = abs(x)*abs(x)
nanvar(X,r) = condmean( f ,X,!(isnan),r)
```
squares the values of `X` before accumulating them, thus computing the variance.

## Numerical precision

I often need to compute small anomalies of 1000s of `Floats`, each of which is near a large mean offset. Summing many large numbers blunts the computed mean's numerical precision, so that it is less precise than the anomalies.

`ConditionalMean` improves numerical precision when averaging large arrays of large floating point numbers by subtracting a representative value before accumulating the sum. It adds the offset back to the mean at the end. The offset is the first valid datum in the sum. This is faster than preconditioning the input data by subtracting an approximation of the mean beforehand, and usually adequate. However, if the first value is large and not representative of the mean, it will make the precision worse. I have never encountered this problem, but if you do, you need to condition the input data set by subtracting its mean, or increase the numerical precision, e.g. from `Float32` to `Float64`.
