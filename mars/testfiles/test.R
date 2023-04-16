library(mars)
#please set working directory with setwd("file path here")
#ex. setwd("C:/Users/bensh/OneDrive/Documents/GitHub/Stat360-Mars/mars/testfiles")
#data obtained from https://archive-beta.ics.uci.edu/dataset/162/forest+fires
df <- read.csv("./forestfires.csv",header=TRUE)
df <- subset(df, select = c(wind, temp, RH))

m <- mars(RH~.,df)

print(m)

anova.mars(m)

plot.mars(m)


#data obtained from https://archive-beta.ics.uci.edu/dataset/165/concrete+compressive+strength

df2 <- read.csv("./Concrete_Data.csv",header=TRUE)
names(df2) <- c("Cement","Blast.Furnace.Slag","Fly.Ash","Water",
                "Superplasticizer","Coarse.Aggregate","Fine.Aggregate.",
                "Age.(day)", "Concrete.compressive.strength")
df2 <- head(df2,500)

m <- mars(Concrete.compressive.strength~.,df2)

print(m)

predict.mars(m, newdata)

anova.mars(m)

plot.mars(m)



#data obtained from https://archive-beta.ics.uci.edu/dataset/368/facebook+metrics
df3 <- read.delim("./dataset_Facebook.csv",sep=";",header=TRUE)

df3 <- subset(df3, select = c(Paid, like, comment, share))

df3 <- head(df3,400)
newdata <- tail(df3,100)
names(newdata) <- names(df3)

m <- mars(Paid~., df3, mars.control(10,trace=TRUE))

print(m)

predict.mars(m, newdata)

manova.mars(m)

plot.mars(m)

