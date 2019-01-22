###
# This is the main script implementing the methodology on the Quebec Database
###


###
# Outline of the key steps:
# 0. Check for missing packages and install them
# 1. Compute the characteristics from the Quebec Database
#    (relying on functions defined in QuebecDB_glue.R)
# 2. Do the correlation analyses
# 3. Infer optimal transition matrices using the Gibbs Sampling
#    (relying on functions defined in MCMC.R)
# 4. Validation of the transition matrix
# 5. Prediction of characteristic states for 2010's, 2020's, 2030's, and long-term
###


###
# Step 0: Check for missing packages and install them
###

pkg.inst <- installed.packages()
pkgs <- c("Hmisc", "MCMCpack", "corrgram", "plotrix", "RColorBrewer")#, "wordcloud")
have.pkg <- pkgs %in% rownames(pkg.inst)
if (any(!have.pkg)) {
  cat("Some packages need to be installed\n")
  r <- readline("Install necessary packages [y/n]? ")
  if (tolower(r) == "y") {
    need <- pkgs[!have.pkg]
    message("installing packages ", paste(need, collapse = ", "))
    install.packages(need)
  }
}
for (pkg in pkgs) {
  library(pkg,character.only=T)
}


###
# Step 1: Compute the characteristics from the Quebec Database
###

source("transition_quebec_estimate_database.R")

RawInv.data <- mdb.get('BD_PEP.mdb') # loads the MDB database
Ref.species <- read.csv("REF_SPECIES_SUCCESSIONAL_QUEBEC.CSV", header=T)  ## loads Jenkins and shade tolerance species-dependend coefficients

TREES <- RawInv.data$ETUDARBR # get the individual trees info from database table ETUDARBR
TREES <- find_years(TREES, RawInv.data$INFOGEN) # parse for the years from database table INFOGEN

CHARACTERISTICS_FULL <- as.data.frame(get.summaries(TREES, RawInv.data))  # this computes all the stand characteristics


###
# Step 2: Do the Correlation Analyses
###

# First, the correlation in the whole database
latex(cor(scale(CHARACTERISTICS_FULL[,3:8]), use="pairwise.complete.obs"),
      ctable=T,landscape=T,digits=2,size="footnotesize",
      caption="Correlations in the whole database")

# Then, the correlation decade per decade
latex(cor(scale(CHARACTERISTICS_FULL[which(CHARACTERISTICS$YEAR>=2000),3:8]),use="pairwise.complete.obs"),
      ctable=T, landscape=T, digits=2, size="footnotesize",
      caption="Correlations in the 2000's")
latex(cor(scale(CHARACTERISTICS_FULL[which((CHARACTERISTICS$YEAR>=1990) & (CHARACTERISTICS$YEAR<2000)),3:8]),use="pairwise.complete.obs"),
      ctable=T, landscape=T, digits=2, size="footnotesize",
      caption="Correlations in the 1990's")
latex(cor(scale(CHARACTERISTICS_FULL[which((CHARACTERISTICS$YEAR>=1980) & (CHARACTERISTICS$YEAR<1990)),3:8]),use="pairwise.complete.obs"),
      ctable=T, landscape=T, digits=2, size="footnotesize",
      caption="Correlations in the 1980's")
latex(cor(scale(CHARACTERISTICS_FULL[ which((CHARACTERISTICS$YEAR>=1970) & (CHARACTERISTICS$YEAR<1980)),3:8]),use="pairwise.complete.obs"),
      ctable=T, landscape=T, digits=2, size="footnotesize",
      caption="Correlations in the 1970's")

# visual output of the correlation matrix
par(mfrow = c(1,1), oma = c(0,0,0,0), mar=c(0,0,0,0))
corrgram(CHARACTERISTICS[,c(3,4,5,6,7,8)], labels=c("Biomass","Basal Area","Ext. Diversity","Int. Diversity","Shade Index","Average Age"), order=TRUE, lower.panel=panel.conf, upper.panel=panel.pie, text.panel=panel.txt,cex=1.5)#,cex.cor=1.5)
dopdf('correlations_visu')


# Finally, the PCA
CHARACTERISTICS_NONA <- na.omit(CHARACTERISTICS_FULL)
pca <- prcomp(CHARACTERISTICS_NONA[,3:8],scale.=T)
screeplot(pca)
latex(pca$rotation,
      ctable=T, landscape=T, digits=2, size="footnotesize",
      caption="PCA loadings")

par(mfrow=c(2,3))
biplot(pca, choices=c(1,2),cex=0.8)
biplot(pca, choices=c(1,3),cex=0.8)
biplot(pca, choices=c(2,3),cex=0.8)
biplot(pca, choices=c(1,4),cex=0.8)
biplot(pca, choices=c(2,4),cex=0.8)
biplot(pca, choices=c(3,4),cex=0.8)
par(mfrow=c(1,1))


###
# Step 3: Infer optimal transition matrices using the Gibbs Sampling
###

source("transition_quebec_estimate_gibbs.R")

cut_high <- 50 # the threshold to keep biomass is 50x10³ kg/ha
prop.cut <- sum(CHARACTERISTICS_NONA$BIOMASS>cut_high) / length(CHARACTERISTICS_NONA$BIOMASS)
print(paste("got rid of", as.character(signif(prop.cut*100, 2)), "% of data above the threshold"))
CHARACTERISTICS_analysis <- CHARACTERISTICS_NONA[-which(CHARACTERISTICS_NONA$BIOMASS>cut_high),]

# We define the breaks defining the states ranges
biomass.breaks <- seq(from=min(CHARACTERISTICS_analysis$BIOMASS), to=max(CHARACTERISTICS_analysis$BIOMASS), length.out=26)
divint.breaks <- seq(from=min(CHARACTERISTICS_analysis$INT.DIV), to=max(CHARACTERISTICS_analysis$INT.DIV), length.out=11)
successional.breaks <- seq(from=min(CHARACTERISTICS_analysis$SHADE.TOLERANCE), to=max(CHARACTERISTICS_analysis$SHADE.TOLERANCE), length.out=11)
age.breaks <- seq(from=min(CHARACTERISTICS_analysis$MEAN.AGE), to=max(CHARACTERISTICS_analysis$MEAN.AGE), length.out=11)

# We compute the transition matrix over the full database (for prediction)
biomass.transitions <- get.transitions(CHARACTERISTICS_analysis, 3, biomass.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
divint.transitions <- get.transitions(CHARACTERISTICS_analysis, 7, divint.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
successional.transitions <- get.transitions(CHARACTERISTICS_analysis, 5, successional.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
age.transitions <- get.transitions(CHARACTERISTICS_analysis, 8, age.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)

# We compute the transition matrix over the first half of the database (for validation)
first_half <- which(CHARACTERISTICS_analysis$YEAR<=1988)
biomass.fh_transitions <- get.transitions(CHARACTERISTICS_analysis[first_half,], 3, biomass.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
divint.fh_transitions <- get.transitions(CHARACTERISTICS_analysis[first_half,], 7, divint.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
successional.fh_transitions <- get.transitions(CHARACTERISTICS_analysis[first_half,], 5, successional.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
age.fh_transitions <- get.transitions(CHARACTERISTICS_analysis[first_half,], 8, age.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)

# We compute the transition matrix over 50% of the data, taken at random over the full span 1971-2007 (for validation)
random_half <- sample.int(nrow(CHARACTERISTICS_analysis),round(nrow(CHARACTERISTICS_analysis)/2))
biomass.rh_transitions <- get.transitions(CHARACTERISTICS_analysis[random_half,], 3, biomass.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
divint.rh_transitions <- get.transitions(CHARACTERISTICS_analysis[random_half,], 7, divint.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
successional.rh_transitions <- get.transitions(CHARACTERISTICS_analysis[random_half,], 5, successional.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
age.rh_transitions <- get.transitions(CHARACTERISTICS_analysis[random_half,], 8, age.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
# And the other complimentary 50%:
random_half2 <- rep(T,nrow(CHARACTERISTICS_analysis))
random_half2[random_half] <- F
random_half2 <- which(random_half2)
biomass.rh2_transitions <- get.transitions(CHARACTERISTICS_analysis[random_half2,], 3, biomass.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
divint.rh2_transitions <- get.transitions(CHARACTERISTICS_analysis[random_half2,], 7, divint.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
successional.rh2_transitions <- get.transitions(CHARACTERISTICS_analysis[random_half2,], 5, successional.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)
age.rh2_transitions <- get.transitions(CHARACTERISTICS_analysis[random_half2,], 8, age.breaks, runs=50, burnin=100, Hmax=1000, Tmax=12, step=3, duplicate=T)


###
# Step 4: Validation of the transition matrix
###

# Step 4.1: study of the transition matrix errors

# we extract the mean of all transitions
biomass.transition <- apply(biomass.transitions,c(2,3),mean)
divint.transition <- apply(divint.transitions,c(2,3),mean)
successional.transition <- apply(successional.transitions,c(2,3),mean)
age.transition <- apply(age.transitions,c(2,3),mean)

# plots the transition matrices
pretty.matrix(biomass.transition)
pretty.matrix(divint.transition)
pretty.matrix(successional.transition)
pretty.matrix(age.transition)

# we extract the standard deviation of all transitions
biomass.error <- apply(biomass.transitions,c(2,3),sd)
divint.error <- apply(divint.transitions,c(2,3),sd)
successional.error <- apply(successional.transitions,c(2,3),sd)
age.error <- apply(age.transitions,c(2,3),sd)


# Step 4.2: prediction of the second half of the dataset with the first half

# we extract the mean of all transitions from the first half of the dataset
biomass.fh_transition <- apply(biomass.fh_transitions,c(2,3),mean)
divint.fh_transition <- apply(divint.fh_transitions,c(2,3),mean)
successional.fh_transition <- apply(successional.fh_transitions,c(2,3),mean)
age.fh_transition <- apply(age.fh_transitions,c(2,3),mean)

# we predict the states for the second half
sh_predicted_repartition.biomass <- do_prediction(years=1989:2007,CHARACTERISTICS_analysis$BIOMASS[first_half],CHARACTERISTICS_analysis$YEAR[first_half],biomass.fh_transition,biomass.breaks)
sh_predicted_repartition.divint <- do_prediction(years=1989:2007,CHARACTERISTICS_analysis$INT.DIV[first_half],CHARACTERISTICS_analysis$YEAR[first_half],divint.fh_transition,divint.breaks)
sh_predicted_repartition.successional <- do_prediction(years=1989:2007,CHARACTERISTICS_analysis$SHADE.TOLERANCE[first_half],CHARACTERISTICS_analysis$YEAR[first_half],successional.fh_transition,successional.breaks)
sh_predicted_repartition.age <- do_prediction(years=1989:2007,CHARACTERISTICS_analysis$MEAN.AGE[first_half],CHARACTERISTICS_analysis$YEAR[first_half],age.fh_transition,age.breaks)

# we get the real distribution for the states of the second half
second_half = which(CHARACTERISTICS_analysis$YEAR>1988)
sh_real_repartition.biomass <- get_repartition(values=CHARACTERISTICS_analysis$BIOMASS[second_half],biomass.breaks,freq=T)
sh_real_repartition.divint <- get_repartition(values=CHARACTERISTICS_analysis$INT.DIV[second_half],divint.breaks,freq=T)
sh_real_repartition.successional <- get_repartition(values=CHARACTERISTICS_analysis$SHADE.TOLERANCE[second_half],successional.breaks,freq=T)
sh_real_repartition.age <- get_repartition(values=CHARACTERISTICS_analysis$MEAN.AGE[second_half],age.breaks,freq=T)


# Step 4.3: prediction of two random halves of the dataset

# we extract the mean of all transitions
biomass.rh_transition <- apply(biomass.rh_transitions,c(2,3),mean)
divint.rh_transition <- apply(divint.rh_transitions,c(2,3),mean)
successional.rh_transition <- apply(successional.rh_transitions,c(2,3),mean)
age.rh_transition <- apply(age.rh_transitions,c(2,3),mean)
biomass.rh2_transition <- apply(biomass.rh2_transitions,c(2,3),mean)
divint.rh2_transition <- apply(divint.rh2_transitions,c(2,3),mean)
successional.rh2_transition <- apply(successional.rh2_transitions,c(2,3),mean)
age.rh2_transition <- apply(age.rh2_transitions,c(2,3),mean)

rh_predicted_repartition.biomass <- matrix(get_equilibrium(biomass.rh_transition),ncol=1)
rh_predicted_repartition.divint <- matrix(get_equilibrium(divint.rh_transition),ncol=1)
rh_predicted_repartition.successional <- matrix(get_equilibrium(successional.rh_transition),ncol=1)
rh_predicted_repartition.age <- matrix(get_equilibrium(age.rh_transition),ncol=1)
rh2_predicted_repartition.biomass <- matrix(get_equilibrium(biomass.rh2_transition),ncol=1)
rh2_predicted_repartition.divint <- matrix(get_equilibrium(divint.rh2_transition),ncol=1)
rh2_predicted_repartition.successional <- matrix(get_equilibrium(successional.rh2_transition),ncol=1)
rh2_predicted_repartition.age <- matrix(get_equilibrium(age.rh2_transition),ncol=1)


###
# Step 5: Prediction of characteristic states for 2010's, 2020's, 2030's, and long-term
###

# we first deduce the long-term equilibrium
predicted_repartition.biomass <- matrix(get_equilibrium(biomass.transition),ncol=1)
predicted_repartition.divint <- matrix(get_equilibrium(divint.transition),ncol=1)
predicted_repartition.successional <- matrix(get_equilibrium(successional.transition),ncol=1)
predicted_repartition.age <- matrix(get_equilibrium(age.transition),ncol=1)

#We then compute explicitely the predictions for years 2010's, 2020's, 2030's based on data from the 2000's
d2000 <- which(CHARACTERISTICS_analysis$YEAR>=2000)
for (years in list(2030:2039,2020:2029,2010:2019)) {
  predicted_repartition.biomass <- cbind(t(do_prediction(years=years,CHARACTERISTICS_analysis$BIOMASS[d2000],CHARACTERISTICS_analysis$YEAR[d2000],biomass.transition,biomass.breaks)),predicted_repartition.biomass)
  predicted_repartition.divint <- cbind(t(do_prediction(years=years,CHARACTERISTICS_analysis$INT.DIV[d2000],CHARACTERISTICS_analysis$YEAR[d2000],divint.transition,divint.breaks)),predicted_repartition.divint)
  predicted_repartition.successional <- cbind(t(do_prediction(years=years,CHARACTERISTICS_analysis$SHADE.TOLERANCE[d2000],CHARACTERISTICS_analysis$YEAR[d2000],successional.transition,successional.breaks)),predicted_repartition.successional)
  predicted_repartition.age <- cbind(t(do_prediction(years=years,CHARACTERISTICS_analysis$MEAN.AGE[d2000],CHARACTERISTICS_analysis$YEAR[d2000],age.transition,age.breaks)),predicted_repartition.age)
}
