library(devtools)
load_all("mars")
#please set working directory with setwd("file path here")
#ex. setwd("C:/Users/bensh/OneDrive/Documents/GitHub/Stat360-Mars/")


# Test 1: Comparing Relative humidity to wind and temperature
#data obtained from https://archive-beta.ics.uci.edu/dataset/162/forest+fires
df <- read.csv("testfiles/forestfires.csv",header=TRUE)
df <- subset(df, select = c(wind, temp, RH))

m <- mars(RH~.,df)

plot.mars(m)
print.mars(m)

# Test 2: Comparing Concrete compressive strength to many other variables
#data obtained from https://archive-beta.ics.uci.edu/dataset/165/concrete+compressive+strength
df2 <- read.csv("testfiles/Concrete_Data.csv",header=TRUE)
names(df2) <- c("Cement","Blast.Furnace.Slag","Fly.Ash","Water",
                "Superplasticizer","Coarse.Aggregate","Fine.Aggregate.",
                "Age.(day)", "Concrete.compressive.strength")
df2 <- head(df2,500)

m <- mars(Concrete.compressive.strength~.,df2)

print.mars(m)

anova.mars(m)



# Test 3: Comparing facebook likes to comments and shares
#data obtained from https://archive-beta.ics.uci.edu/dataset/368/facebook+metrics
df3 <- read.delim("testfiles/dataset_Facebook.csv",sep=";",header=TRUE)

df3 <- subset(df3, select = c(like, comment, share))

df3 <- head(df3,400)
newdata <- tail(df3,100)
names(newdata) <- names(df3)

m <- mars(like~., df3, mars.control(10,trace=TRUE))

predict.mars(m, newdata)

anova.mars(m)

