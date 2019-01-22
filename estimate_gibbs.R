###
# The functions here implement the MCMC algorithm
###


# initialization of sequences: complete missing data
initialize_missing <- function(S,n,method="determ",overall_distrib=NA)
{
  y <- S
  Tmax = dim(S)[2]
  
  if (method == "rnd") {
    if (is.na(overall_distrib)) {
      ## ## ## completely random
      y[missing_data] <- sample(1:n, sum(missing_data), replace=T)
    } else {
      ## ## ## random but so that the distribution of new data matches the distribution of existent data
      for (t in 1:Tmax) {
        for (missing in which(missing_data[,t]==T)) {
          y[missing,t] <- which(overall_distrib >= runif(1))[1]
        }
      }
    }
  } else if (method == "determ") {
    ## ## ## determinist: new samples are the same that the last known sample
    last = y[,1]
    for (t in 2:Tmax) {
      missing = is.na(y[,t])
      y[missing,t] = last[missing]
      last[!missing] = y[!missing,t]
    }
  } else {
    stop("method not recognized")
  }
  
  return(y)
}


# complete sequences using the Gibbs sampling
gibbs_sampling <- function(S,n,Hmax=200,burnin=100,gamma=1)
{ 
  # initialization of various variables
  m = dim(S)[1]
  Tmax = dim(S)[2]
  missing_data <- is.na(S)
  transition.h = matrix(NA,nrow=n,ncol=n)
  omega.h = matrix(NA,nrow=n,ncol=n)
  transitions = array(NA,dim=c(Hmax-burnin,n,n))
  
  # compute the overall distribution
  overall_distrib = rep(0,n)
  previous = 0
  for (i in 1:n) {
    overall_distrib[i] = sum(S==i,na.rm=T) + previous
    previous = overall_distrib[i]
  }
  overall_distrib = overall_distrib/overall_distrib[n]
  
  # initialize the sequences
  y <- initialize_missing(S,n)
  
  pb <- txtProgressBar(style = 3)
  for (h in 1:Hmax) {
    ########################
    # parameter estimation #
    ########################
    omega.h <- get_statistics(y,n)
    
    for (i in 1:n) {
      transition.h[i,] <- rdirichlet(1,1+omega.h[i,])
    }
    
    #####################
    # data augmentation #
    #####################
    for (t in 2:(Tmax-1)) {
      for (missing in which(missing_data[,t]==T)) {
        proba <- transition.h[y[missing,t-1],] * transition.h[,y[missing,t+1]]
        if (sum(proba,na.rm=T)==0) {
          # this check is for a rare but nasty bug with very low values for:
          # transition.h[y[missing,t-1],] and transition.h[,y[missing,t+1]]
          # resulting in sum(proba)==0|NA => divide by 0 throws the error
          y[missing,t] <- which(overall_distrib >= runif(1))[1]
        } else {
          y[missing,t] <- which(cumsum(proba/sum(proba)) >= runif(1))[1]
        }
        if (is.na(y[missing,t])) {browser()}
      }
    }    
    if (h>burnin) {
      transitions[h-burnin,,] = transition.h
    }
    setTxtProgressBar(pb, h / Hmax)
  }
  
  return(transitions)
}

# compute the statistics indicator of the sequence
get_statistics  <- function(S,n)
{
  m = dim(S)[1]
  Tmax = dim(S)[2]
  # estimate from the sequences a transition matrix
  omega <- matrix(0,ncol=n,nrow=n)
  for (i in 1:n) {
    for (j in 1:n) {
      for (t in 2:Tmax) {
        omega[i,j] <- omega[i,j] + sum((S[,t-1] == i & S[,t] == j),na.rm=T)
      } 
    }  
  }
  return(omega)
}

# when duplicate==True, each year is considered both as a first measurement and as a re-measurement
isolate_repeated_measurements = function(y_id_indic, duplicate=T)
{
  total_length = dim(y_id_indic)[1]
  data.repeated = cbind(y_id_indic[,2],10*y_id_indic[,1],y_id_indic[,3],rep(1,total_length))
  i = 1
  pb <- txtProgressBar(style = 3)
  while (i < dim(data.repeated)[1]) {
    if (data.repeated[i,4]==0) {i=i+1;next}
    j = which((data.repeated[,2]==data.repeated[i,2]) & (data.repeated[,4]==1))
    if (length(j) == 0) {
      # measures already parsed before: do nothing
      i = i + 1
    } else if (length(j) == 1) {
      data.repeated = data.repeated[-j,] # get rid of the unique measure
      #i = i + 1
    } else {
      data.repeated[j,4] = 0 # we won't parse them anymore
      if (duplicate) {
        # case duplicate: will duplicate as much as possible data
        stored_j = j # will use to delete at the end
        subid = 1
        while (length(j) > 1) {
          m = which(data.repeated[j,1] == min(data.repeated[j,1]))
          data.repeated[j,1] = data.repeated[j,1] - data.repeated[j[m],1]
          for (jj in 1:length(j)) {
            data.repeated = rbind(data.repeated,c(data.repeated[j[jj],1]+1,data.repeated[j[jj],2]+subid,data.repeated[j[jj],3],0))
          }
          j = j[-m]
          subid = subid + 1
        }
        data.repeated = data.repeated[-stored_j,]
      } else {
        # case no-duplicate: won't complete the data with duplicates
        data.repeated[j,1] = data.repeated[j,1] - min(data.repeated[j,1]) + 1
        i = i + 1
      }
    }
    setTxtProgressBar(pb, i / dim(data.repeated)[1])
  }
  
  return(as.data.frame(cbind(val=data.repeated[,3],year=data.repeated[,1],id=data.repeated[,2])))
}

# builds the sequences from yearly data and given breaks
construct_sequences = function(data, states.nb=15, Tmax=NA, breaks=NA)
{
  if (is.na(Tmax)) {
    Tmax = max(data$year)
  }
  r = data[order(data$id,data$year),]
  r = r[which(r$year <= Tmax),]
  len = dim(r)[1]
  if (any(is.na(breaks))) {
    b = seq(min(r$val), max(r$val), by=(max(r$val)-min(r$val))/states.nb) # this defines the states
  } else {
    b = breaks
  }
  b[length(b)] = b[length(b)]+1 # to make sure we include data in the last bin
  S = matrix(NA,ncol=Tmax,nrow=1)
  pb <- txtProgressBar(style = 3)
  i = 1
  while(i <= len) {
    s = rep(NA,Tmax)
    ii = which(r$id == r$id[i])
    for(y in 1:Tmax) {
      datai = which(r$year[ii] == y)
      if (length(datai) > 1) {
        browser("debug breakpoint: this should never happen")
      } else if (length(datai) == 1) {
        # there is a data for this year
        s[y] = which(b > r$val[ii[datai]])[1] -1 # keep the state, not the value
      }
    }
    if (sum(!is.na(s[2:Tmax]))>=1 & !is.na(s[1])) {
      S = rbind(S,s)      
    }
    i = max(ii)+1
    len = dim(r)[1]
    setTxtProgressBar(pb, i / len)
  }
  return(S[2:dim(S)[1],])
}


compress_sequence = function(S,step=2)
{
  Tmax = dim(S)[2]
  m = dim(S)[1]
  cS = matrix(NA,ncol=round(Tmax/step),nrow=m)
  toremove = NULL
  for (i in 1:m) {
    t = 1
    while (t<Tmax) {
      t2 = min(c(Tmax,t+step-1))
      if (sum(!is.na(S[i,t:t2])) == 1) {
        cS[i,(t+step-1)/step] = na.omit(S[i,t:t2])        
      } else if (sum(!is.na(S[i,t:t2])) > 1) {
        print(paste("Note: skipping NA entry at i=",i))
        toremove = c(toremove,i)
      }
      t = t + step
    }
  }
  if (length(toremove)) {
    cS = cS[-toremove,]
  }
  return(cS)
}

# main function - workflow of the approach
get.transitions <- function(CHARACTERISTICS, c_indic, breaks, runs=5, burnin=100, Hmax=500, Tmax=12, step=3, duplicate=T)
{
  dup <- isolate_repeated_measurements(CHARACTERISTICS[,c(1,2,c_indic)], duplicate=duplicate)
  cat("\nConstructing sequences:\n")
  S_dup <- construct_sequences(dup, states.nb=n, Tmax=Tmax, breaks=breaks)
  cat(paste("\nDiscretizing into ranges of",step,"years:\n"))
  cS_dup <- compress_sequence(S_dup, step=step)
  l = 1
  pooled.transitions = array(NA,dim=c((Hmax-burnin)*runs,length(breaks)-1,length(breaks)-1))
  for (i in 1:runs) {
    cat(paste("\nFull MCMC search (", i, "/", runs,") with", burnin, "ignored iterations followed by", Hmax-burnin, "used iterations:\n"))
    transitions <- gibbs_sampling(cS_dup, length(breaks)-1, Hmax=Hmax, burnin=burnin)
    for (j in 1:nrow(transitions)){
      pooled.transitions[l,,] <- transitions[j,,]
      l = l + 1
    }
  }
  return (pooled.transitions)
}

# empirical deduction of equilibrium state
get_equilibrium = function(transition,n=1000,min_n=10,epsilon=1e-10,plot=F)
{
  for (i in 1:n) {
    transition = transition %*% transition
    sumdif = 0
    for (j in 1:nrow(transition)) {
      sumdif = sumdif + (max(transition[,j]) - min(transition[,j]))
    }
    if (sumdif < epsilon & i >= min_n) {
      cat(paste('stopped after',i,'iterations'))
      break
    }
  }
  if (plot==T) {
    pretty.matrix(transition)
  }
  return (transition[1,])
}

# multiply a given state with the transition matrix, `step` times`
get_prediction = function(step=1,transition,repartition)
{
  for (s in 1:step) {
    repartition = repartition %*% transition
  }
  
  return(repartition)
}

# get the distributions of values within breaks (cut)
get_repartition = function(values,breaks,freq=F)
{
  repartition = NULL
  for (i in 2:length(breaks)) {
    repartition[i-1] = sum( (values>=breaks[i-1]) & (values<breaks[i]))
  }
  # the following accounts for the values exactly equal to the last
  # breaks (so the state intervals are [a,b[, [b,c[, ..., [y,z] )
  i=length(breaks)
  repartition[i-1] = repartition[i-1] + sum(values==breaks[i])
  if (freq)
    return(repartition/sum(repartition))
  else
    return(repartition)
}

# predict using a given transition matrix, timeframe, and initial state
do_prediction = function(years=2010:2019,CHARACTERISTICS.VARIABLE,CHARACTERISTICS.YEAR,transition,breaks)
{
  repartition = rep(0,length.out=length(breaks)-1)
  db.years = unique(CHARACTERISTICS.YEAR)
  for (y1 in db.years) {
    relevant_plots = which(CHARACTERISTICS.YEAR==y1)
    for (y2 in years) {
      if ((y2-y1)/3==round((y2-y1)/3)) {
        # we do the prediction
        #         cat(paste('doin the prediction with year',y1,'to year',y2,'\n'))
        repartition = repartition + get_prediction(step=(y2-y1)/3,transition,get_repartition(CHARACTERISTICS.VARIABLE[relevant_plots],breaks))
      }
    }
  }
  return(repartition/sum(repartition))
}


# function to output a transition matrix with probabilities color-coded
pretty.matrix <- function(m,round.nb=2, col=NA, doplot=T, labels=NA, normalize_color=NA, topleft=NA,
                          #                           par.mar=c(0.5,0.5,0.5,0.5),
                          par.mar=c(2,2,2,2),
                          colorbar=NA,
                          top=NA, left=NA, mtext.cex=2, colormap=heat.colors(100)[100:1], ...)
{
  write.table(format(round(m,digits=round.nb), justify="right"), row.names=F, col.names=F, quote=F)
  
  if (doplot) {
    if (!is.na(colorbar)) {
      layout(matrix(c(1,2),ncol=2),widths=c(100-colorbar,colorbar))
    } else {
      par(mfrow = c(1,1))
    }
    par(oma = c(1,1,1,1), mar=par.mar)
    if (any(is.na(col))) {
      m.image = m
      zlim=c(min(m),max(m))
    } else {
      m.image = col
      zlim=c(min(col),max(col))
    } 
    if (!any(is.na(normalize_color))) {
      m.image[which(m.image>=normalize_color)] = normalize_color
      zlim=c(0,normalize_color)
    }
    if (length(labels)>1) {
      if (any(is.na(normalize_color))) {
        image(3+(1:ncol(m.image)), 1+(1:nrow(m.image)), t(m.image[nrow(m.image):1,]),xlim=c(1.5,(3.5+ncol(m.image))), ylim=c(1.5,(3.5+nrow(m.image))),col=colormap, ylab="",xlab="",xaxt="n", yaxt="n")        
      } else {
        image(3+(1:ncol(m.image)), 1+(1:nrow(m.image)), t(m.image[nrow(m.image):1,]),xlim=c(1.5,(3.5+ncol(m.image))), ylim=c(1.5,(3.5+nrow(m.image))),zlim=zlim, col=colormap, ylab="",xlab="",xaxt="n", yaxt="n")
      }
      m2 = matrix(NA,nrow=nrow(m.image)+2,ncol=ncol(m.image)+2)
      m2[1:2,] = 1
      m2[,1:2] = 1
      #m2[3:nrow(m2),3:ncol(m2)] = m
      image(1+(1:ncol(m2)), 1+(1:nrow(m2)), t(m2[nrow(m2):1,]),add=T,xlim=c(1.5,(1.5+ncol(m))), ylim=c(1.5,(1.5+nrow(m))),zlim=c(0,1), col=rgb(0.8,0.8,0.8), ylab="",xlab="",xaxt="n", yaxt="n")
    } else {
      image((1:ncol(m.image)), (1:nrow(m.image)), t(m[nrow(m.image):1,]), zlim=zlim, col=colormap, ylab="",xlab="",xaxt="n", yaxt="n")
    }
    axis(side=1, labels=F, at=NULL, tick=F)
    axis(side=2, labels=F, at=NULL, tick=F)
    for (x in 1:ncol(m)) {
      if (length(labels)>1) {
        text(2.5, nrow(m)-x+2, labels[x],...)
        text(x+3, nrow(m)+2.5, labels[x],srt=60,...)
        for (y in 1:nrow(m)) {
          text(x+3, nrow(m)-y+1+1, round(m[y,x],digits=round.nb),...)
        }
      } else {
        for (y in 1:nrow(m)) {
          text(x, nrow(m)-y+1, round(m[y,x],digits=round.nb),...)
        }
      }
    }
  }
  if (!is.na(topleft)) {
    text(2.5,nrow(m)-1+3.5,topleft,...)
  }
  box()
  if (!is.na(top)) {
    mtext(top,font=2,cex=mtext.cex,padj=-0.5)
  }
  if (!is.na(left)) {
    mtext(left,side=2,font=2,cex=mtext.cex,padj=-0.5)
  }
  
  if (!is.na(colorbar)) {
    mar = par("mar")
    mar[2]=1
    mar[4]=3
    par(mar=mar)
    labelbreaks=round(seq(zlim[1],zlim[length(zlim)],length.out=3+1),digits=round.nb)
    image.scale(z=zlim,labelbreaks=labelbreaks,labelbreaks.str=paste(labelbreaks,'%'),horiz=F,col=colormap,mtext.cex=mtext.cex)
    box()
  }
  print(paste("max value is",as.character(max(m))))
}
