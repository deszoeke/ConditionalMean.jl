# ConditionalMean

Provides general methods for summing and averaging values that meet a specific condition, over various dimensions of Julia arrays. A useful example is to average finite values and ignore missing values.

## Missing values

Julia inherits functionality and different standards for missing values from many languages. Examples:
1. `missing`
2. DataFrames `NA` type
3. IEEE `NaN`
4. `Int64(-999)`
5. any negative value

The desired behavior for handing missing values is almost always one of 
1. The result of any calculation involving any missing value(s) is a missing value.
2. Perform the calculation ignoring missing values. MATLAB's `nanmean()` does this.

You might want to read data formatted with one missing value convention, process it with methods using another, and maybe even plot data with methods using a third convention. `ConditionalMean` gives you tools to ignore missing values minimizing or eliminating the need to convert custom types at each step.

### Examples
```
nanmean(X,r) = condmean(X, !(isnan), r)
```
averages all values that are not `NaN` over the region `r`.

You can also call `condmean` on a function of the array. For example,
```
f(x) = abs(x)*abs(x)
nanvar(X,r) = condmean( f ,X,!(isnan),r)
```
squares the values of `X` before accumulating them, thus computing the variance.
