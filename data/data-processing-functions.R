# functions to process data for open robust design model

# from Kery and Schaub 2012
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}


# create vector with occasion of marking
get.first <- function(x) min(which(x!=0))

# generate all 3 m-arrays for ORD model

# marray_ord takes a data frame with 5 columns:
# flag - unique individual ID
# year - year of observation
# sample - sampling occasion of observation
# type - resight or capture
# firstyear - logical, T/F year of first marking

marray_ord = function(dat){
  
  # across years
  ch_yr = as.matrix(table(dat$flag, dat$year))
  ch_yr[ch_yr > 1] = 1
  ch_yr = ch_yr[-1,]
  marr.p = marray(ch_yr)
  rel = rowSums(marr.p)
  
  # within years
  # note: "dummy" flags may be added during pre-processing to ensure that 
  # there is a row in the data for every year-sample combination 
  # for year-sample occs without real observations, detection probability should be
  # constrained accordingly in the model
  dat_sec <- dat %>% 
    filter(!first_year | flag == "dummy") %>% 
    filter(type == "resight" | flag == "dummy")
  
  n.years = length(unique(dat_sec$year))
  n.sec = length(unique(dat_sec$sample))
  years = unique(dat_sec$year)
  
  all.marr.r = all.marr.t = array(NA, dim = c(n.years, n.sec-1, n.sec))
  rel.sec.t = rel.sec.r = matrix(0, nrow = n.years, ncol = n.sec-1)
  new.s = matrix(0, nrow = n.years, ncol = n.sec)
  
  for(j in 1:length(years)){
    
    yrdat = filter(dat_sec, year == years[j])
    ch = as.matrix(table(yrdat$flag, yrdat$sample))
    ch[ch > 1] = 1
    ch = ch[-1,]
    
    age = rep(1, nrow(ch))
    mAge = 2
    maxAge <- max(c(max(age), mAge))
    nind <- nrow(ch)
    n.occasions <- ncol(ch)
    first <- apply(ch, 1, get.first)
    age.matrix <- matrix(0, ncol = n.occasions, nrow = nind)
    for (i in 1:nind){
      age.matrix[i,first[i]:n.occasions] <- 1:(n.occasions-first[i]+1)+(age[i]-1)
    }
    age.matrix[age.matrix > maxAge] <- maxAge
    
    # Recode capture history
    ch.rec <- ch
    for (i in 1:nind){
      h <- which(ch.rec[i,]==1)
      for (k in 1:length(h)){
        ch.rec[i,h[k]] <- k
      } # j
    } # i
    ch.rec[ch.rec > maxAge] <- maxAge
    ch.split <- array(0, dim = c(nrow(ch), ncol(ch), maxAge))
    
    for (a in 1:maxAge){
      for (i in 1:nind){
        k <- which(ch.rec[i,]==a | ch.rec[i,]==(a+1))
        if (length(k)==0) next
        ch.split[i,k[1:2],age.matrix[i,k[1]]] <- 1
        if (length(k)>1){
          ch.split[i,k[2:length(k)],age.matrix[i,k[2]]] <- 1
        }
      } # i
    } # a
    
    all.marr.t[j,,] <- marray(ch.split[,,1])
    all.marr.r[j,,] <- marray(ch.split[,,2])
    
    rel.sec.t[j,] <- rowSums(all.marr.t[j,,])
    rel.sec.r[j,] <- rowSums(all.marr.r[j,,])
    
    if(length(table(first)) < n.sec){
      occ = as.numeric(names(table(first)))
      to.add = which(!c(1:n.sec) %in% occ)
      new = numeric(n.sec)
      for(i in occ){
        new[occ] <- table(first)[as.character(occ)]
      }
      new[to.add] <- 0
    } else new <- as.numeric(table(first))
    
    new.s[j,] <- new
  }
  
  dat_out = list(
    marr.p = marr.p,
    new.s = new.s,
    marr.s.r = all.marr.r,
    marr.s.t = all.marr.t,
    rel = rel,
    rel.sec.r = rel.sec.r,
    rel.sec.t = rel.sec.t,
    tot = rowSums(new.s)
  )
  
  return(dat_out)
}


