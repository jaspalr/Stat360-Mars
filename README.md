# Stat360-Mars
R: Multivariate Adaptive Regression Splines (MARS)    const macros = { "\\\\R": "\\\\textsf{R}", "\\\\code": "\\\\texttt"}; function processMathHTML() { var l = document.getElementsByClassName('reqn'); for (let e of l) { katex.render(e.textContent, e, { throwOnError: false, macros }); } return; } 

mars

R Documentation

Multivariate Adaptive Regression Splines (MARS)
-----------------------------------------------

### Description

Fit Friedman's Multivariate Adaptive Regression Splines (MARS) model. https://github.com/jaspalr/Stat360-Mars

### Usage

    first create a mars control object using mars.control then
    mars(y~.,data, marscontrol)
    

### Arguments

`formula`

an regression formula

`data`

a data frame containing the data

`control`

a mars control object (see mars.control to build one)

### Details

This function creates a Multivariate Adaptive Regression Splines model for set of data with a forward step and backward step based using Friedman's method

### Value

a mars class containing the call, input formula, y values, coefficients, and data frame containing the split points, and direction for a given coefficient

### Author(s)

Jaspal Raman, Ben Shires Nakamura, Jessica Kim

### References

Multivariate Adaptive Regression Splines

### See Also

anova.mars, print.mars, plot.mars, predict.mars, print.mars, summary.mars

### Examples

     #Example 1: Comparing Relative humidity to wind and temperature
    #data obtained from https://archive-beta.ics.uci.edu/dataset/162/forest+fires
    df <- read.csv("testfiles/forestfires.csv",header=TRUE)
    df <- subset(df, select = c(wind, temp, RH))
    
    m <- mars(RH~.,df)
    plot.mars(m)
    print.mars(m)
    # Example  2: Comparing Concrete compressive strength to many other variables
    #data obtained from https://archive-beta.ics.uci.edu/dataset/165/concrete+compressive+strength
    df2 <- read.csv("testfiles/Concrete_Data.csv",header=TRUE)
    names(df2) <- c("Cement","Blast.Furnace.Slag","Fly.Ash","Water",
                   "Superplasticizer","Coarse.Aggregate","Fine.Aggregate.",
                  "Age.(day)", "Concrete.compressive.strength")
    df2 <- head(df2,500)
    m <- mars(Concrete.compressive.strength~.,df2)
    print.mars(m)
    anova.mars(m)
    # Example  3: Comparing facebook likes to comments and shares
    #data obtained from https://archive-beta.ics.uci.edu/dataset/368/facebook+metrics
    df3 <- read.delim("testfiles/dataset_Facebook.csv",sep=";",header=TRUE)
    df3 <- subset(df3, select = c(like, comment, share))
    df3 <- head(df3,400)
    newdata <- tail(df3,100)
    names(newdata) <- names(df3)
    m <- mars(like~., df3, mars.control(10,trace=TRUE))
    predict.mars(m, newdata)
    anova.mars(m)
