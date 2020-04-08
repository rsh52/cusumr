# cusumr

`cusumr` is a R package dedicated to analyzing binary data sets for cumulative summation. There are two main functions in this package:

- `cusumr()`
- `cusumr_plot()`

Both functions have defaults set for acceptable rates, type 1 and type 2 error rates, whether or not a learning period should be required in the analysis, and if resets to the data should be applied based on decision limits in the monitoring phase.

### Example:

df <- c(0,0,0,1,1,0,1,1,0,0,0,0,1,0,0,1)

cusumr(df)
cusumr_plot(df)
