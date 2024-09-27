####### Script to identify spine clusters and spine types #######
## 0. required packages
## 1. spine morphology - spine types + normalisation
## 2. spine number & distribution - cluster + normalisation
## 3. example genearlised linear model taking nested data structure into account 


#####
### 0. required packages
install.packages("readxl", dependencies = TRUE)
install.packages("tidyverse", dependencies = TRUE)

install.packages("reshape2", dependecies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("generics", dependencies = TRUE)

install.packages("lme4", dependencies = TRUE)
install.packages("nlme", dependencies = TRUE)
install.packages("emmeans", dependencies = TRUE)

install.packages("scatterplot3d")
install.packages("cluster")
install.packages("factoextra")
install.packages("corrplot")

# Load packages
library(tidyverse)
library(readxl)

library(reshape2)
library(dplyr)
library(generics)

library(lme4) # recommeneded to run generalised linear models with nested data structure
library(nlme) #nonlinear mixed-effects models
library(emmeans)

library(scatterplot3d) # making 3D scatter plot plot 
library(cluster) # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(corrplot) # making correlation plots for PCA analysis

rm(list = ls())
graphics.off()


# there are two independent variables in this data set reactivation status and condition
# condition (CFC vs CE)
# reactivation status (cFos+RFP+ vs cFos-RFP+)
# !please adapt code to your experimental design! 

#adust working directory to location in which data files are stored
setwd("~/Documents/Data_Compiler_Example/Processed_data/Output")


#####
### 1. spine morphology - spine types + normalisation
rm(list = ls())
graphics.off()

# read data file
df <- read_excel("morph_data.xlsx")

df <- rename(df, spine_mean_head = `Spine Part Mean Diameter Head`)
df <- rename(df, spine_max_head = `Spine Part Max Diameter Head`)
df <- rename(df, spine_volume_head = `Spine Part Volume Head`)
df <- rename(df, spine_mean_neck = `Spine Part Mean Diameter Neck`)
df <- rename(df, spine_max_neck = `Spine Part Max Diameter Neck`)
df <- rename(df, spine_volume_neck = `Spine Part Volume Neck`)
df <- rename(df, spine_length = `Spine Length`)
df <- rename(df, spine_volume = `Spine Volume`)

# create new variable based on ratio of spine head and neck volume
df$head_neck_ratio <- df$spine_volume_head / df$spine_volume_neck

# create empty vector for for loop
df$spine_type <- vector(mode = "numeric", length = length(df$spine_length))

# classify each spine as stubby/mushroom/thin/unknown type based on specified criteria
for (i in 1:length(df$spine_length)) {
  if ((df$spine_length[i] <= 1 & df$head_neck_ratio[i] < 1.2) | df$spine_max_neck[i] == 0) {
    df$spine_type[i] <- "stubby"
  } else if (df$spine_length[i] <= 5 & df$head_neck_ratio[i] >= 1.2) {
    df$spine_type[i] <- "mushroom"
  } else if (df$spine_length[i] <= 5 & df$head_neck_ratio[i] < 1.2) {
    df$spine_type[i] <- "thin"
  } else {
    df$spine_type[i] <- "unknown"
  }
}

# factorise variable for later analyses
df$spine_type <- as.factor(df$spine_type)


# counts the number of spines on each dendrite
n_spines_dendrite <- df %>% count(dendrite, sort = TRUE)
names(n_spines_dendrite)[names(n_spines_dendrite) == "n"] <- "n_spines_type"
df <- merge(df, n_spines_dendrite, by = c("dendrite"), all = TRUE)

# counts the number of spine types on each dendrite
df1 <- df %>% 
  group_by(animal, image, dendrite, condition, reactivation_status,  n_spines_type) %>% 
  count(spine_type)
names(df1)[names(df1) == "n"] <- "frequency"

# calculates the ratio of each spine type to the total amount of spines on each dendrite
df1$relative_frequ <- df1$frequency / df1$n_spines_type


#####
### 2. spine number & distribution - cluster + normalisation
## sort data according to dendrite ID and to distance from dendrite beginning point 
rm(list = ls())
graphics.off()

# read data file
df <- read_excel("distr_dendrite_data.xlsx")

df <- rename(df, spine_attachment = `Spine Attachment Pt Distance`)
df <- rename(df, spine_attachment_diameter = `Spine Attachment Pt Diameter`)
df <- rename(df, spine_density = `Dendrite Spine Density`)
df <- rename(df, dendrite_length = `Dendrite Length`)

df <- df[
  with(df, order(dendrite, spine_attachment)),
]


## calculate interspine distance by calculating distance between spine and spine before in df
df$lag_spine_attachment <- dplyr::lag(df$spine_attachment, n = 1) # creates new variable with spine attachment ahead of spine
df$interspine_distance <- df$spine_attachment - df$lag_spine_attachment # calculates interspine distance between current spine and spine before
df$interspine_distance <- ifelse(df$interspine_distance < 0, NA, df$interspine_distance) #if distance is negative, it must be a previous dendrite -> set distance to NA

## calculate interspine distance to two spines before spine in df
df$lag_spine_attachment2 <- dplyr::lag(df$spine_attachment, n = 2) # creates new variable with spine attachment two spines ahead of spine

df$distance_to_twoneighbours_before <- (df$spine_attachment - df$lag_spine_attachment2) # calculates interspine distance between spine and TWO spines before

df$cluster <- ifelse(df$distance_to_twoneighbours_before <= 1.5 & df$distance_to_twoneighbours_before >= 0, 1, 0) # if interspine distance between spine and two spines before is larger than 0 and smaller than 1.5, these three spines are defined as a cluster, otherwise not
df$cluster[is.na(df$cluster)] <- 0 

df$cluster_before <- dplyr::lag(df$cluster, n = 1) # creates new variable with cluster information before current spine
df$cluster_before[is.na(df$cluster_before)] <- 0

df_omit <- na.omit(df)


###cluster
count_overall <- 0
count_size <- 0
i <- 0

## total number of clusters
# counts the total amount of clusters across the entire df 
for (i in 1:length(df_omit$cluster)) {
  if (df_omit$cluster[i] == 1 & df_omit$cluster_before[i] == 0) {
    count_overall <- count_overall + 1
  } else {
    count_overall <- count_overall
  }
}

## size of all clusters together
# 
for (i in 1:length(df_omit$cluster)) {
  if (df_omit$cluster[i] == 1 & df_omit$cluster_before[i] == 0) {
    count_overall <- count_overall + 1
    count_size <- count_size + 1
  } else if ((df_omit$cluster[i] == 1 & df_omit$cluster_before[i] == 1)) {
    count_overall <- count_overall
    count_size <- count_size + 1
  } else {
    count_overall <- count_overall
    count_size <- 0
  }
  df_omit$count_size[i] <- count_size
}

## size of each cluster
local_maxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}


df_omit$maxima <- 0
df_omit$maxima[local_maxima(df_omit$count_size)] <- df_omit$count_size[local_maxima(df_omit$count_size)] ## extracts the largest number of spines within a cluster, to determine size of cluster 
df_omit$cluster_size <- ifelse(df_omit$maxima > 0, df_omit$maxima + 2, 0) # each cluster contains of 3 spines, therefore cluster size is 1 + 2


# counts the number of spines on each dendrite
n_spines_dendrite <- df_omit %>% count(dendrite, sort = TRUE) 
names(n_spines_dendrite)[names(n_spines_dendrite) == "n"] <- "n_spines_cluster"
df_omit <- merge(df_omit, n_spines_dendrite, by = c("dendrite"), all = TRUE)


#counts the number of clusters consisting of 3 spines across each dendrite
df2 <- df_omit %>%
  group_by(animal, image, dendrite,  condition, reactivation_status, n_spines_cluster) %>% 
  count(cluster) 
names(df2)[names(df2) == "n"] <- "frequency"

df2$animal <- as.factor(df2$animal)
df2$cluster <- as.factor(df2$cluster)

# adding correct lables to cluster variable
levels(df2$cluster)[levels(df2$cluster)=="0"] <- "no cluster" 
levels(df2$cluster)[levels(df2$cluster)=="1"] <- "cluster"

# counts the number of clustered and non-clustered spines along each dendrite
df2_w <- dcast(df2, animal + dendrite + image + condition + reactivation_status + n_spines_cluster ~ cluster, value.var = "frequency")
df2_w <- df2_w %>% replace(is.na(df2_w), 0)

# calculates ratio fo clustered to non-clustered spines along each dendrite
df2_w$ratio <- (df2_w$cluster / df2_w$`no cluster`) 



#####
### 3. example genearlised linear model taking nested data structure into account 
# taking df2_w from part 2 of this script

df2_w$ratio <- df2_w$ratio + 0.5 # Gamma distribution does not work on values = 0 -> add a constant to all values

hist(df2_w$ratio) # investigate distribution of dependent variable
stat.desc(df2_w$ratio, basic = F, norm = T) # check whether dependent variable is normally distributed
describeBy(ratio ~ condition + reactivation_status, data = df2_w) # get descriptive statistics for each group

# fixed intercept model
cluster_model0 <- glm(ratio ~ 1, data = df2_w, family = Gamma(link="inverse")) 
summary(cluster_model0) #AIC -3

# random intercept model
cluster_model1 <- glmer(ratio ~ 1 + (1 |animal/image), data = df2_w, family = Gamma(link="inverse")) 
summary(cluster_model1) #AIC -56
cluster_model2 <- glmer(ratio ~ 1 + (1 |image), data = df2_w, family = Gamma(link="inverse")) 
summary(cluster_model2) #AIC -57 -> smallest AIC, use this intercept for predictor model
cluster_model3 <- glmer(ratio ~ 1 + (1 |animal), data = df2_w, family = Gamma(link="inverse")) 
summary(cluster_model3) #AIC -32

# compare fit of random vs fixed intercept model 
anova(cluster_model2, cluster_model0, test="Chisq") #significant, continue with predictor model

# predictor model 
cluster_model4 <- glmer(ratio ~ 1 + condition * reactivation_status + (1 |image), data = df2_w, family = Gamma(link="inverse")) 
summary(cluster_model4) # AIC -60 

#compare fit of predictor vs random intercept model 
anova(cluster_model4, cluster_model2, test="Chisq") #significant, final model can be accepted

emmeans(cluster_model4, pairwise ~ reactivation_status * condition, adjust = "fdr") # post-hoc comparisons


# check assumptions
# normality of residuals
qqPlot(residuals(cluster_model4))
# majority of values falls within prediction, no violation of normality

# homoscedasticity
df2_w$spine_res <- residuals(cluster_model4) # extracts the residuals and places them in a new column in our original data table
df2_w$spine_res <- abs(df2_w$spine_res) # creates a new column with the absolute value of the residuals
df2_w$spine_res2 <- df2_w$spine_res^2 # squares the absolute values of the residuals to provide the more robust estimate
Levene.model <- lm(spine_res2 ~  reactivation_status * condition, data = df2_w) # ANOVA of the squared residuals
anova(Levene.model) # displays the results
# non-significant, no violation of homoscedasticity

# multicollinearity
vif(cluster_model4)
# all values < 5, no violation of multicollinearity 

# linearity and homogeneity of residuals
plot(resid(cluster_model4, type = "pearson") ~ fitted(cluster_model4),
     xlab = "Fitted values", ylab = "Pearson residuals"
)
abline(h = 0, col = "red")
# small values not perfectly randomly distributed, other variable might lead to this pattern in certain values
# BUT all other assumptions are met
# --> accept final model! 

