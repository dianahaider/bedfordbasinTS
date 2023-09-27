# Install the required packages
install.packages("vegan")
install.packages("labdsv")
install.packages("MASS")
install.packages("ggplot2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")
library(BiocManager)
BiocManager::install("microbiome")

# install mvpart from package archive file
install.packages("remotes")
remotes::install_url("https://cran.r-project.org/src/contrib/Archive/mvpart/mvpart_1.6-2.tar.gz")

# Load the required packages
library(labdsv)
library(vegan)
library(MASS)
library(mvpart)
library(ggplot2)
library(microbiome)  

# Import df of samples x features
spe <- read.csv("data/shallowdfs.csv", row.names = 1)

#Import env data
env <- read.csv("data/env_data_shallow_clean.csv", row.names = 1)

dim(spe) # dataset dimensions

str(spe)  # structure of objects in dataset

# Count number of species frequencies in each abundance class
ab <- table(unlist(spe))

# Count the number of zeros in the dataset
sum(spe == 0)

# Calculate proportion of zeros in the dataset
sum(spe == 0)/(nrow(spe) * ncol(spe))

#clr transform the data here
speclr <- microbiome::transform(spe, 'clr')

#explore the env data
names(env)
dim(env)
head(env)

str(env)
summary(env)

# We can visually look for correlations between variables:
heatmap(abs(cor(env)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))


# Scale and center variables
env.z <- decostand(env, method = "standardize")

# Variables are now centered around a mean of 0
round(apply(env.z, 2, mean), 1)

# Model the effect of all environmental variables on fish
# community composition
spe.rda <- rda(speclr ~ ., data = env.z)

summary(spe.rda)

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(speclr ~ 1, data = env.z), # lower model limit (simple!)
                      scope = formula(spe.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = TRUE) # change to TRUE to see the selection process!


rda(formula = speclr ~ temperature + Nitrite + POC, data = env.z)







