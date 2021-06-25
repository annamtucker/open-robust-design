# helper functions 


# perform GoF testing with R2ucare
gof_which_mod = function(ch){
  
  caphx = group_data(ch, rep(1, dim(ch)[1]))
  hx = caphx[,1:(ncol(caphx)-1)]
  freq = caphx[,ncol(caphx)]
  
  overall = overall_CJS(hx, freq)
  
  test_3sr = test3sr(hx, freq)
  test_3sm = test3sm(hx, freq)
  test_2ct = test2ct(hx, freq)
  test_2cl = test2cl(hx, freq)
  
  test = list(test_2ct, test_3sr, test_3sm, test_2cl)
  tests = c("test_2ct", "test_3sr", "test_3sm", "test_2cl")
  sig = c()
  for(j in 1:4){
    sigtest = ifelse(test[j][[1]][[1]][3] < 0.05, tests[j], "")
    sig = c(sig, " ", sigtest)
  }
  
  # overall test
  model = "standard CJS model"
  if(overall$p_value < 0.05){
    
    # test 2CT = yes
    if(test[1][[1]][[1]][3] < 0.05){
      
      # test 3SR = yes
      if(test[2][[1]][[1]][3] < 0.05){
        
        # test 3SM or 2CL = yes
        if(test[3][[1]][[1]][3] < 0.05 | test[4][[1]][[1]][3] < 0.05){
          model = "trap-dependence, transience, overdispersion"
        }else model = "trap-dependence and transience"
      }
      
      # test 3SR = no
      if(test[2][[1]][[1]][3] > 0.05){
        
        # test 3SM or 2CL = yes
        if(test[3][[1]][[1]][3] < 0.05 | test[4][[1]][[1]][3] < 0.05){
          model = "trap-dependence and overdispersion"
        } else model = "trap-dependence"
      }
    }
    
    # test 2CT = no
    if(test[1][[1]][[1]][3] > 0.05){
      
      # test 3SR = yes
      if(test[2][[1]][[1]][3] < 0.05){
        
        # test 3SM or 2CL = yes
        if(test[3][[1]][[1]][3] < 0.05 | test[4][[1]][[1]][3] < 0.05){
          model = "transience and overdispersion"
        }else model = "transience"
      }
      
      # test 3SR = no
      if(test[2][[1]][[1]][3] > 0.05){
        
        # test 3SM or 2CL = yes
        if(test[3][[1]][[1]][3] < 0.05 | test[4][[1]][[1]][3] < 0.05){
          model = "overdispersion"
        }
      }
    }
  }
  
  return(list(overall$p_value, paste(sig, collapse = ""), model))  
  
}

Rhat = function(x){
  
  M = length(x)
  N = dim(x[[1]])[1]
  
  theta_m = sapply(x, function(x) apply(x, 2, mean))
  var_m = sapply(x, function(x) apply(x, 2, var))
  theta = apply(theta_m, 1, mean)
  theta_dif_sq = (theta_m - theta)^2
  
  B = apply(theta_dif_sq, 1, function(x) (N/(M-1))*sum(x))
  W = apply(var_m, 1, mean)
  
  V = (((N-1)/N)*W) + ((M+1)/(M*N)*B)
  
  d = N - 1
  
  Rhat = sqrt(((d+3)/(d+1))*(V/W))
  
  return(Rhat)
  
}

safe_Rhat <- safely(Rhat, otherwise = NA_real_)


extract_safe_rhat = function(x){
  y = x$rhat
  y$result
}

extract_safe_parms = function(x){
  names(x$rhat$result)
}

extract_rhat = function(x){
  x$rhat
}

extract_parms = function(x){
  names(x$rhat)
}


edit_colnames =function(z){
  x = z
  cnames = colnames(x$chain1)
  cnames = gsub(" ", "", cnames)
  
  colnames(x$chain1) <- cnames
  colnames(x$chain2) <- cnames
  colnames(x$chain3) <- cnames
  return(x)
}

spread_results = function(x){
  x %>% 
    spread_draws(mean.gammaII, mean.gammaOI,
                 mean.tau, mean.phi, mean.p,
                 sig.gammaII, sig.gammaOI, 
                 sig.phi, sig.tau, sig.p,
                 var.gammaII, var.gammaOI, 
                 var.tau, var.phi,
                 pstar[year], phi[year], p[year,occ], gammaII[year], gammaOI[year],
                 psi[year,occ], tau[year], delta[year,occ])
}


est_to_df = function(x){
  parmnames = gsub(" ", "", rownames(x))
  parm = ifelse(str_detect(parmnames, "\\["),
                substr(parmnames, 1, str_locate(parmnames, "\\[")-1),
                parmnames)
  year = as.numeric(ifelse(parm %in% na.year, NA,
                           ifelse(parm %in% na.occ, substr(parmnames, 
                                                           str_locate(parmnames, "\\[")+1,
                                                           str_locate(parmnames, "\\]")-1),
                                  substr(parmnames, str_locate(parmnames, ",")+1,
                                         str_locate(parmnames, "\\]")-1))))
  occ = as.numeric(ifelse(parm %in% na.occ, NA,
                          substr(parmnames, str_locate(parmnames, "\\[")+1,
                                 str_locate(parmnames, ",")-1)))
  tibble(parm, year, occ) %>% 
    bind_cols(as.tibble(x))
}

mcmc_to_df = function(x, parms){
  
  parmnames = colnames(x$chain1)
  
  parm = ifelse(str_detect(parmnames, "\\["),
                substr(parmnames, 1, str_locate(parmnames, "\\[")-1),
                parmnames)
  
  keep = which(parm %in% parms)
  parm = parm[keep]
  parmnames = parmnames[keep]
  
  na.year = parm[!str_detect(parmnames, "\\[")]
  na.occ = parm[!str_detect(parmnames, ",")]
  reps = dim(x$chain1)[1]
  
  year = as.numeric(ifelse(parm %in% na.year, NA,
                           ifelse(parm %in% na.occ, substr(parmnames, 
                                                           str_locate(parmnames, "\\[")+1,
                                                           str_locate(parmnames, "\\]")-1),
                                  substr(parmnames, str_locate(parmnames, "\\[")+1,
                                         str_locate(parmnames, ",")-1))))
  
  occ = as.numeric(ifelse(parm %in% na.occ, NA,
                          substr(parmnames, str_locate(parmnames, ",")+1,
                                 str_locate(parmnames, "\\]")-1)))
  
  c1 <- tibble(parm = rep(parm, each = reps),
         year = rep(year, each = reps),
         occ = rep(occ, each = reps),
         val = c(x$chain1[,keep]),
         draw = rep(c(1:reps),length(parmnames)))
  
  c2 <- tibble(parm = rep(parm, each = reps),
               year = rep(year, each = reps),
               occ = rep(occ, each = reps),
               val = c(x$chain2[,keep]),
               draw = rep(c((reps+1):(reps*2)),length(parmnames)))
               
  c3 <- tibble(parm = rep(parm, each = reps),
               year = rep(year, each = reps),
               occ = rep(occ, each = reps),
               val = c(x$chain3[,keep]),
               draw = rep(c(((reps*2)+1):(reps*3)), length(parmnames)))
  
  c1 %>% 
    rbind(c2) %>% 
    rbind(c3)
}




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


# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))



marray_age <- function(ch, age, mAge = 1){
  # Input variables
  # ch: matrix with capture histories.
  # Note: the capture history file is a single file that includes the individuals of all age classes
  # age: vector with the age class at first capture for each individual
  # mAge: maximal number of age classes for which m-arrays are constructed. 
  # Input is optional and only required if the age matrix has fewer age classes 
  # as we want to separate (e.g. CH contains only individuals marked as juveniles, and we want 2 age classes)
  #
  # Output
  # marr: 3-d array with the m-array. The third dimension is the age class. The last column of each m-array is the number of released individuals that were never recaptured. Thus, the total number of released individuals per occasion is the row sum of the m-array.
  # 1. Define subfunctions
  # 1.1. Function to create a m-array from capture-histories (ch)
  marray <- function(ch, unobs = 0){
    ns <- length(table(ch)) - 1 + unobs
    no <- ncol(ch)
    out <- matrix(0, ncol = ns*(no-1)+1, nrow = ns*(no-1))
    # Remove capture histories of ind. that are marked at last occasion
    get.first <- function(x) min(which(x!=0))
    first <- apply(ch, 1, get.first)
    last <- which(first==no)
    if (length(last) > 0) ch <- ch[-last,]
    # Compute m-array
    for (i in 1:nrow(ch)){
      cap.occ <- which(ch[i,]!=0)
      state <- ch[i,cap.occ]
      if (length(state) == 1) {
        out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] <- out[state[1]+ns*(cap.occ[1]-1), ns*(no-1)+1] + 1
      }
      if (length(state) > 1) {
        for (t in 2:length(cap.occ)){
          out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*ns+state[t]] <- out[(cap.occ[t-1]-1)*ns+state[t-1], (cap.occ[t]-2)*ns+state[t]] + 1
        } #
        if (max(cap.occ) < no){
          out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] <- out[(cap.occ[t]-1)*ns+state[t], ns*(no-1)+1] + 1
        } # if
      } # if
    } # t
    return(out)
  }
  # 1.2. Function to remove capture histories without any capture
  clean.ch <- function(ch){
    incl <- which(rowSums(ch)>=1)
    ch <- ch[incl,]
    return(ch)
  }
  # 1.3. Function to remove all first captures from the capture-histories
  rm.first <- function(ch) {
    get.first <- function(x) min(which(x==1))
    first <- apply(ch, 1, get.first)
    for (i in 1:nrow(ch)){
      ch[i,first[i]] <- 0
    }
    return(ch)
  }
  # 1.4. Function to calculate the occasion of first capture
  get.first <- function(x) min(which(x==1))
  # 2. Calculations
  if (is.matrix(ch)==FALSE) ch <- matrix(ch, nrow = 1)
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
    for (j in 1:length(h)){
      ch.rec[i,h[j]] <- j
    } # j
  } # i
  ch.rec[ch.rec > maxAge] <- maxAge
  ch.split <- array(0, dim = c(nrow(ch), ncol(ch), maxAge))
  for (a in 1:maxAge){
    for (i in 1:nind){
      j <- which(ch.rec[i,]==a | ch.rec[i,]==(a+1))
      if (length(j)==0) next
      ch.split[i,j[1:2],age.matrix[i,j[1]]] <- 1
      if (length(j)>1){
        ch.split[i,j[2:length(j)],age.matrix[i,j[2]]] <- 1
      }
    } # i
  } # a
  marr <- array(0, dim = c(n.occasions-1, n.occasions, maxAge))
  for (a in 1:(maxAge-1)){
    for (i in 1:nind){
      u <- which(ch.split[i,,a]==1)
      if (length(u)==0) next
      if (u[1]==n.occasions) next
      if (length(u)==1) marr[u,n.occasions,a] <- marr[u,n.occasions,a] + 1
      if (length(u)==2) marr[u[1],u[2]-1,a] <- marr[u[1],u[2]-1,a] + 1
    } # i
  } # a
  a <- maxAge
  if (is.matrix(ch.split[,,a])==FALSE){
    ch.split1 <- matrix(ch.split[,,a], nrow = 1)
    marr[,,a] <- marray(ch.split1)
  } # if
  else marr[,,a] <- marray(ch.split[,,a])
  return(marr)
}



# function to process shorebird data for ORD model

marray_ord_2 = function(dat){
  
  # across years
  ch_yr = as.matrix(table(dat$flag, dat$year))
  ch_yr[ch_yr > 1] = 1
  ch_yr = ch_yr[-1,]
  
  marr.p = marray_age(ch_yr, age = rep(1, nrow(ch_yr)), mAge = 2)
  
  marr.p.t = marr.p[,,1]
  marr.p.r = marr.p[,,2]
  
  rel.t = rowSums(marr.p.t)
  rel.r = rowSums(marr.p.r)
  
  # within years
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
    
    # marr <- array(0, dim = c(n.sec-1, n.sec, maxAge))
    # for (i in 1:nind){
    #   u <- which(ch.split[i,,1]==1)
    #   if (length(u)==0) next
    #   if (u[1]==n.occasions) next
    #   if (length(u)==1) marr[u,n.occasions,1] <- marr[u,n.occasions,1] + 1
    #   if (length(u)==2) marr[u[1],u[2]-1,1] <- marr[u[1],u[2]-1,1] + 1
    # }
    # 
    # marr[,,2] <- marray(ch.split[,,2])
    
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
    marr.p.t = marr.p.t,
    marr.p.r = marr.p.r,
    rel.t = rel.t,
    rel.r = rel.r,
    new.s = new.s,
    marr.s.r = all.marr.r,
    marr.s.t = all.marr.t,
    rel.sec.r = rel.sec.r,
    rel.sec.t = rel.sec.t,
    tot = rowSums(new.s)
  )
  
  return(dat_out)
}



marray_ord = function(dat){
  
  # across years
  ch_yr = as.matrix(table(dat$flag, dat$year))
  ch_yr[ch_yr > 1] = 1
  ch_yr = ch_yr[-1,]
  marr.p = marray(ch_yr)
  rel = rowSums(marr.p)
  
  # within years
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
    
    # marr <- array(0, dim = c(n.sec-1, n.sec, maxAge))
    # for (i in 1:nind){
    #   u <- which(ch.split[i,,1]==1)
    #   if (length(u)==0) next
    #   if (u[1]==n.occasions) next
    #   if (length(u)==1) marr[u,n.occasions,1] <- marr[u,n.occasions,1] + 1
    #   if (length(u)==2) marr[u[1],u[2]-1,1] <- marr[u[1],u[2]-1,1] + 1
    # }
    # 
    # marr[,,2] <- marray(ch.split[,,2])
    
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



marray_ord_day = function(dat){
  
  # across years
  ch_yr = as.matrix(table(dat$flag, dat$year))
  ch_yr[ch_yr > 1] = 1
  ch_yr = ch_yr[-1,]
  marr.p = marray(ch_yr)
  rel = rowSums(marr.p)
  
  # within years
  dat_sec <- dat %>% 
    filter(!first_year | flag == "dummy") %>% 
    filter(type == "resight" | flag == "dummy")
  
  n.years = length(unique(dat_sec$year))
  n.sec = length(unique(dat_sec$day))
  years = unique(dat_sec$year)
  
  all.marr.r = all.marr.t = array(NA, dim = c(n.years, n.sec-1, n.sec))
  rel.sec.t = rel.sec.r = matrix(0, nrow = n.years, ncol = n.sec-1)
  new.s = matrix(0, nrow = n.years, ncol = n.sec)
  
  for(j in 1:length(years)){
    
    yrdat = filter(dat_sec, year == years[j])
    ch = as.matrix(table(yrdat$flag, yrdat$day))
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
    
    # marr <- array(0, dim = c(n.sec-1, n.sec, maxAge))
    # for (i in 1:nind){
    #   u <- which(ch.split[i,,1]==1)
    #   if (length(u)==0) next
    #   if (u[1]==n.occasions) next
    #   if (length(u)==1) marr[u,n.occasions,1] <- marr[u,n.occasions,1] + 1
    #   if (length(u)==2) marr[u[1],u[2]-1,1] <- marr[u[1],u[2]-1,1] + 1
    # }
    # 
    # marr[,,2] <- marray(ch.split[,,2])
    
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





