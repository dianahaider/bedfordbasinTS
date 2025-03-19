# Install the required packages
#install.packages("vegan")
#install.packages("labdsv")
#install.packages("MASS")
#install.packages("ggplot2")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.17")
#BiocManager::install("microbiome")

# install mvpart from package archive file
#install.packages("remotes")
remotes::install_url("https://cran.r-project.org/src/contrib/Archive/mvpart/mvpart_1.6-2.tar.gz")

# Load the required packages
library(BiocManager)
library(labdsv)
library(vegan)
library(MASS)
library(mvpart)
library(ggplot2)
library(microbiome)  

# Import df of samples x features
spe <- read.csv("data/subyears_tables.csv", row.names = 1)

#Import env data
env <- read.csv("data/env_data_deep.csv", row.names = 1)

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



# Write our new model
spe.rda.signif <- rda(speclr ~ Temperature_y + conductivity + Chlorophyll.A +
                        oxygen + Nitrate + 
                        Phosphate , data = env.z)
# check the adjusted R2 (corrected for the number of
# explanatory variables)
RsquareAdj(spe.rda.signif)


anova.cca(spe.rda.signif, step = 1000)

anova.cca(spe.rda.signif, step = 1000, by = "term")

anova.cca(spe.rda.signif, step = 1000, by = "axis")


env$year <- as.factor(env$year)
View(env)

# Type 1 scaling
ordiplot(spe.rda.signif, scaling = 1, type='n')
points (spe.rda.signif, col = env$month, pch=as.integer(env$year) )
legend("bottomright", legend=unique(env$year), pch=unique(env$year) )
legend("topleft", legend=unique(env$month), col=unique(env$month), pch= 15 )


# Type 2 scaling
ordiplot(spe.rda.signif, scaling = 2, type='n')
points (spe.rda.signif, col = env$month, pch = as.integer(env$year))
text(spe.rda.signif,display="cn",cex=.8,col="blue")
legend("bottomright", legend=unique(env$year), pch=unique(env$year) )
legend("topright", legend=unique(env$month), col=unique(env$month), pch= 15 )



Bedford_clr <- microbiome::transform(Bedford_no_mito, "clr")   
out.pcoa.logt <- ordinate(Bedford_clr, method = "RDA", distance = "euclidean")
evals <- out.pcoa.logt$CA$eig
p3<-plot_ordination(Bedford_clr, out.pcoa.logt, type = "Sample", 
                    color = "Month" ) 





################
## extract % explained by the first 2 axes
perc <- round(100*(summary(spe.rda.signif)$cont$importance[2, 1:2]), 2)

## extract scores - these are coordinates in the RDA space
sc_si <- scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(spe.rda.signif, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(spe.rda.signif, display="bp", choices=c(1, 2), scaling=1)

## Custom triplot, step by step

# Set up a blank plot with scaling, axes, and labels
plot(spe.rda.signif,
     scaling = 1, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-1,1), 
     ylim = c(-1,1),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = "steelblue", # fill colour
       cex = 1.2) # size
# add points for species scores
points(sc_sp, 
       pch = 22, # set shape (here, square with a fill colour)
       col = "black",
       bg = "#f2bd33", 
       cex = 1.2)
# add text labels for species abbreviations
text(sc_sp + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_sp), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)
# add arrows for effects of the expanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2], # end them at the score value
       col = "red", 
       lwd = 3)
# add text labels for arrows
text(x = sc_bp[,1] -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)

