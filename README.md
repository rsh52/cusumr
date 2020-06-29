# cusumr

## Overview
`cusumr` is a R package dedicated to analyzing binary data sets for cumulative summation. There are two main functions in this package:

- `cusumr()`
- `cusumr_plot()`

Both functions have defaults set for acceptable rates, type 1 and type 2 error rates, whether or not a learning period should be required in the analysis, and if resets to the data should be applied based on decision limits in the monitoring phase.

## Installation

The recommended method for installing `cusumr` is the following:

`devtools::install_github(repo = "rsh52/cusumr")`

### Example

```
df <- c(0,0,0,1,1,0,1,1,0,0,0,0,1,0,0,1)

cusumr(df)

cusumr_plot(df)
```

![cusumr_plot Output]("/images/cusumr_plot.png")

## Reference
The use of the Cusum technique in the assessment of trainee competence in new procedures. Int J Qual Health Care. 2000 Oct;12(5):433-8. [url](https://www.ncbi.nlm.nih.gov/pubmed/11079224)

Cumulative sum (CUSUM) assessment and medical education: a square peg in a round hole. Anaesthesia, 2011, 66, pages 243-254. [url](https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-2044.2011.06692.x)

An application of the learning curve cumulative summation test to evaluate training for endotracheal intubation in emergency medicine. Emerg Med J, 2015;32:291..294. [url](https://www.ncbi.nlm.nih.gov/pubmed/24154942)
