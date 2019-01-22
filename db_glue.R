###
# The functions here compute characteristics from the Quebec Database
###


#This function concatenates trees with a new column containing the measurement year and a unique plot id (from INFOGEN table)
find_years=function(etudarbr,infogen)
{ 
  dates <- structure(numeric(dim(etudarbr)[1]), class="numeric") # initialize "dates"
  ids <- rep(NA,dim(etudarbr)[1]) # initialize "ids"
  cat("Will now parse the INFOGEN table to get the years of measurments:\n")
  pb <- txtProgressBar(style = 3)
  n <- nrow(etudarbr)
  for (j in 1:n) {
    setTxtProgressBar(pb, j / n)
    i = which(infogen$ID.PEP.MES==etudarbr$ID.PEP.MES[[j]])
    if (length(i) > 1) {
      print("duplicates entries...")
      dates[j] = NA
      ids[j] = NA
      next
    } else if (length(i) == 0) {
      print("no date available...")
      dates[j] = NA
      ids[j] = NA
      next
    }
    yy = substr(as.character(infogen$DATE.SOND[[i]]),7,8)
    if (as.numeric(yy) < 14) {
      # we are in the 20xx years
      dates[j] = as.numeric(paste0('20', yy))
      ids[j] = as.numeric(paste0('20', yy, as.character(infogen$ID.PEP[[i]])))
    } else {
      # we are in the 19xx years
      dates[j] = as.numeric(paste0('19', yy))
      ids[j] = as.numeric(paste0('19', yy, as.character(infogen$ID.PEP[[i]])))
    }
  }
  result = cbind(etudarbr,INVYR=dates, UNIQUEPLOTID=ids)
  return(result)
}


#This function sorts and finds the relevant values from the Reference document
ref.values.upd=function(species.vector, xy)
{
  result=matrix(c(species.vector,rep(NA,3*length(species.vector))),ncol=4);
  for(i in 1:length(xy[,1])) {
    # parse the reference DB
    result[as.character(species.vector)==as.character(xy$QUEBEC_NAME[i]),2]=xy$JENKINS_TOTAL_B1[i];
    result[as.character(species.vector)==as.character(xy$QUEBEC_NAME[i]),3]=xy$JENKINS_TOTAL_B2[i];
    result[as.character(species.vector)==as.character(xy$QUEBEC_NAME[i]),4]=xy$SUCCESSIONAL_ID[i];
  }
  return(result)
}


## This is the main function that calculates plot summaries needed
get.summaries = function(trees,Raw,xy=Ref.species)
{
  # we can not use a direct comparison in addressing, as it results in NA when x$ETAT==NA, and x[NA,whatever] is going to fill the table with NA values
  y <- trees[which(trees$ETAT == 10),c(which(colnames(trees)=="UNIQUEPLOTID"),
                            which(colnames(trees)=="INVYR"),
                            which(colnames(trees)=="ESSENCE"),
                            which(colnames(trees)=="DHPMM"),
                            which(colnames(trees)=="HAUTEUR"),
                            which(colnames(trees)=="AGE")
  )]
  colnames(y) <- c("UNIQUEPLOTID","INVYR","SPCD_QUEBEC","DIA","ACTUALHT","AGE")
  
  # store the Jenkins' B1 and B2, as well as the successional id
  B1B2Succ <- ref.values.upd(y$SPCD_QUEBEC,xy)
  y <- y[!is.na(B1B2Succ[,4]),]
  y <- droplevels(y)
  B1B2Succ <- B1B2Succ[!is.na(B1B2Succ[,4]),]
  # infuriating and vexating question of units in Jenkins formula: DIA should be in cm; the resulting BIOMASS is in kg
  # as we want to have the biomass in 10³kg/ha, we scale according to the plot size (1/16 ha)
  BIOMASS <- (exp(B1B2Succ[,2]+B1B2Succ[,3]*log(y$DIA/10)))*16/1000 # required for the Quebec Database, as no dry biomass is available
  y <- cbind(y,BIOMASS=BIOMASS,B1B2Succ=B1B2Succ[,4])
  results <- NULL
  years <- sort(unique(y$INVYR))
  
  cat("Will now compute stand level characteristics for each measurment:\n")
  pb <- txtProgressBar(style = 3)
  
  for(year in years) {
    setTxtProgressBar(pb, (year-years[1]) / (years[length(years)]+1-years[1]))
    z <- y[which(y$INVYR==year),]
    plot.number <- unique(z$UNIQUEPLOTID);
    for (i in plot.number) {
      relevant <- z[which(z$UNIQUEPLOTID==i),];
      relevant <- droplevels(relevant)
      nr.trees <- nrow(relevant)
      average.age <- mean(relevant$AGE, na.rm=T) # note: this results in NA when no tree has a recorded age
      stand.biomass <- sum(relevant$BIOMASS)
      stand.basal.area <- sum(pi*relevant$DIA^2) * 16 / (10^6)# this is in m²/ha
      nr.species <- length(unique(relevant$SPCD_QUEBEC))
      count.trees <- table(relevant$SPCD_QUEBEC)
      prop.count.trees <- count.trees/sum(count.trees)
      simpson.number <- 1 - sum(prop.count.trees^2)
      succesional.number <- sum(relevant$B1B2Succ)/nr.trees
      #this puts it all together
      results <- rbind(results,c(
        as.numeric(substring(as.character(i),5)),
        as.numeric(as.character(year)),
        as.numeric(stand.biomass),
        as.numeric(stand.basal.area),
        as.numeric(succesional.number),
        as.numeric(nr.species),
        as.numeric(simpson.number),
        as.numeric(as.character(average.age))
      ))
    }
  }
  
  # we normalize the number of species comparatively to the maximal number of species in the whole DB
  # so we manipulate an indice between 0 and 1
  results[,6] <- results[,6]/max(results[,6])
  
  colnames(results)=c("ID","YEAR","BIOMASS","BASAL.AREA","SHADE.TOLERANCE","EXT.DIV","INT.DIV","MEAN.AGE")
  return(as.data.frame(results))
}
