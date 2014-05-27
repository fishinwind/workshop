library(ggplot2)
library(plyr)
library(dplyr)

dfx <- read.table('misc/data/expr-geno-covs.txt', header=TRUE)

# calculate some simple stats with plyr and dplyr
# with plyr
ddply(dfx, .(condition, genotype), summarise, mean.age = mean(age), count = n())

# with dplyr
grouped <- group_by(dfx, condition, genotype)
summarize(grouped, count = n(), mean.age = mean(age))