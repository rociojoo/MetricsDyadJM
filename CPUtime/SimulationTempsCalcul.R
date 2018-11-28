#
library('pdist')
pdist.pdist <- pdist::pdist
library('polyclip', 'colorspace')
library(geoR)
library(parallel)
library(pbmcapply)

# library(microbenchmark)
#
#
pre.metrics <- function(A, B) {
  distances <- sqrt(rowSums(A - B) ^ 2)
  longueur <- length(A[, 1])
  distance1 <- sqrt((diff(A[, 1])) ^ 2 + (diff(A[, 2])) ^ 2)
  distance2 <- sqrt((diff(B[, 1])) ^ 2 + (diff(B[, 2])) ^ 2)
  return(
    list(
      distances = distances,
      longueur = longueur,
      distance1 = distance1,
      distance2 = distance2
    )
  )
}

Prox <- function(A, B, delta) {
  pre <- pre.metrics(A, B)
  return(sum(as.numeric(pre$distances < delta)) / pre$longueur)
}

# Prox(A,B,1)

Cs.function <- function(A, B) {
  pre <- pre.metrics(A, B)
  do <- mean(pre$distances)
  distances.de <- as.matrix(pdist.pdist(X = A, Y = B))
  de <- mean(distances.de)
  Cs <- (de - do) / (de + do)
  return(Cs)
}

# Cs.function(A,B)

pre.lixn.HAI <- function(A, B, puntos.area = NULL) {
  A.sAB <- pnt.in.poly(A, as.data.frame(puntos.area))
  B.sAB <- pnt.in.poly(B, as.data.frame(puntos.area))
  longueur <- dim(A)[1]
  nA0 = sum(A.sAB$pip == 1 & B.sAB$pip == 0)
  n0B = sum(A.sAB$pip == 0 & B.sAB$pip == 1)
  nAB = sum(A.sAB$pip == 1 & B.sAB$pip == 1)
  n00 = sum(A.sAB$pip == 0 & B.sAB$pip == 0)
  pAB = ((sum(A.sAB$pip == 1)) / longueur) * ((sum(B.sAB$pip == 1)) / longueur)
  pA0 = ((sum(A.sAB$pip == 1)) / longueur) * (1 - (sum(B.sAB$pip == 1)) /
                                                longueur)
  p0B = (1 - (sum(A.sAB$pip == 1)) / longueur) * ((sum(B.sAB$pip == 1)) /
                                                    longueur)
  p00 = (1 - (sum(A.sAB$pip == 1)) / longueur) * (1 - (sum(B.sAB$pip == 1)) /
                                                    longueur)
  ind.in = which(A.sAB$pip == 1 & B.sAB$pip == 1)
  return(
    list(
      nA0 = nA0,
      n0B = n0B,
      nAB = nAB,
      n00 = n00,
      pA0 = pA0,
      p0B = p0B,
      pAB = pAB,
      p00 = p00,
      points.in = ind.in
    )
  )
}

LixnT <- function(A, B, puntos.area = NULL) {
  pre <- pre.lixn.HAI(A, B, puntos.area)
  if (pre$pA0 == 0 && pre$p0B == 0 && pre$p00 == 0) {
    if (pre$pAB == 0) {
      lixnT <- 0
    } else{
      lixnT <- 1
    }
  } else if (pre$pA0 == 0 && pre$p0B == 0) {
    lixnT <- NA
  } else{
    lixn <-
      log(((pre$nAB / pre$pAB) + (pre$n00 / pre$p00)) / ((pre$nA0 / pre$pA0) +
                                                           (pre$n0B / pre$p0B)))
    lixnT <- 1 / (1 + exp(-lixn))
  }
  return(lixnT)
}

# LixnT(A,B,puntos.area=Z)

HAI <- function(A, B, puntos.area = NULL, delta) {
  pre <- pre.lixn.HAI(A, B, puntos.area)
  pre2 <- pre.metrics(A[pre$points.in, ], B[pre$points.in, ])
  K <- sum(as.numeric(pre2$distances < delta))
  return(K / (K + (pre$nA0 + pre$n0B) / 2))
}

# HAI(A,B,puntos.area=Z,1)

getEllipse <-
  function (X_, Y_, D_, npts = 300) {
    ## computing ellipses for jPPA
    x0 = mean(X_)
    y0 = mean(Y_)
    a <- D_ / 2
    c2 <- (diff(X_) ^ 2 + diff(Y_) ^ 2) / 2
    if (a ^ 2 - c2 < 0) {
      print('Error: max dist is too small regarding the distance between points. Change it. ')
    } else{
      b <- sqrt(a ^ 2 - c2)
      phi <- atan(diff(Y_) / diff(X_))
      
      theta <- seq(0, 2 * pi, length = npts)
      xh <-  a * cos(theta)
      yh <-  b * sin(theta)
      co <- cos(phi)
      si <- sin(phi)
      x  <- x0 + co * xh - si * yh
      y  <- y0 + si * xh + co * yh
      # x  <- x0 + xh
      # y  <- y0 + yh
      
      return(list(x = x, y = y))
    }
  }

jPPA <- function(A,
                 B,
                 Dmax,
                 factor.color = 1,
                 tam.cel = 3) {
  # in contrast to the jppa function in the wildlifeTG package, there is no need to have ltraj objects, nor spatial polygons
  # by gridding the space to compute the areas, the function efficiently computes the intersection and union areas without encountering problems related to disjoint sets
  XA <- A[, 1]
  YA <- A[, 2]
  XB <- B[, 1]
  YB <- B[, 2]
  testA <-
    lapply(1:(length(XA) - 1), function(i) {
      getEllipse(X_ = XA[c(i, i + 1)], Y_ = YA[c(i, i + 1)], D_ = Dmax)
    })
  testB <-
    lapply(1:(length(XB) - 1), function(i) {
      getEllipse(X_ = XB[c(i, i + 1)], Y_ = YB[c(i, i + 1)], D_ = Dmax)
    })
  range.x <-
    range(do.call(rbind, testA)[, 1], do.call(rbind, testB)[, 1], na.rm = TRUE)
  range.y <-
    range(do.call(rbind, testA)[, 2], do.call(rbind, testB)[, 2], na.rm = TRUE)
  xgrid = seq(from = range.x[1], to = range.x[2], by = tam.cel)
  ygrid = seq(from = range.y[1], to = range.y[2], by = tam.cel)
  grilla.inter <- grilla.union <- matrix(0,
                                         ncol = length(xgrid),
                                         nrow = length(ygrid),
                                         byrow = FALSE)
  test.2 <- sapply(1:length(testA), function(i) {
    UnionAB <-
      polyclip(
        testA[[i]],
        testB[[i]],
        op = "union",
        fillA = 'positive',
        fillB = 'positive',
        x0 = mean(c(XA[i], XB[i])),
        y0 = mean(c(YA[i], YB[i])),
        eps = 1e-5
      )
    InterAB <-
      polyclip(
        testA[[i]],
        testB[[i]],
        op = "intersection",
        fillA = 'positive',
        fillB = 'positive',
        x0 = mean(c(XA[i], XB[i])),
        y0 = mean(c(YA[i], YB[i])),
        eps = 1e-5
      )
    if (length(InterAB) > 0) {
      poly.inter.sub <-
        polygrid(
          xgrid = seq(
            from = range.x[1],
            to = range.x[2],
            by = tam.cel
          ),
          ygrid = seq(
            from = range.y[1],
            to = range.y[2],
            by = tam.cel
          ),
          borders = InterAB[[1]],
          vec = TRUE
        )
      toto <-
        matrix(
          as.numeric(poly.inter.sub$vec.inout),
          ncol = length(xgrid),
          nrow = length(ygrid),
          byrow = FALSE
        )
      grilla.inter <<- toto + grilla.inter
    }
    test <- sapply(1:length(UnionAB), function(j) {
      poly.union.sub <-
        polygrid(
          xgrid = seq(
            from = range.x[1],
            to = range.x[2],
            by = tam.cel
          ),
          ygrid = seq(
            from = range.y[1],
            to = range.y[2],
            by = tam.cel
          ),
          borders = UnionAB[[j]],
          vec = TRUE
        )
      toto <-
        matrix(
          as.numeric(poly.union.sub$vec.inout),
          ncol = length(xgrid),
          nrow = length(ygrid),
          byrow = FALSE
        )
      grilla.union <<- toto + grilla.union
    })
    
  })
  grilla.den <- sum(grilla.union > 0)
  grilla.num <- sum(grilla.inter > 0)
  
  jPPA <- grilla.num / grilla.den
  
  return(jPPA)
}

# jPPA(A = A,B=B,Dmax=7)

CSEM.function <- function(A, B, delta) {
  pre <- pre.metrics(A, B)
  m <- 1
  PM <- 1
  P <- NULL
  while (PM[m] != 0 & m < pre$longueur) {
    seq1 <- 1:(pre$longueur - m)
    seq2 <- seq1 + m
    pairs <- as.matrix(cbind(seq1, seq2))
    P[[m]] <-
      as.numeric(sapply(1:length(seq1), function(x) {
        max(pre$distances[pairs[x, 1]:pairs[x, 2]])
      }) < delta)
    PM <- c(PM, sum(P[[m]]))
    (m <- m + 1)
  }
  PM <- PM[-1]
  
  if (PM[length(P)] == 0) {
    CSEM <- (length(P) - 1) / (pre$longueur - 1)
  } else{
    CSEM <- (length(P)) / (pre$longueur - 1)
  }
  return(CSEM)
}

# CSEM.function(A,B,delta)

DI.function <-
  function(A, B, delta.DI) {
    # adapted from the wildlifeDI package without requiring ltraj objects
    pre <- pre.metrics(A, B)
    DId0 <-
      1 - (abs(pre$distance1 - pre$distance2) / (pre$distance1 + pre$distance2)) ^
      delta.DI
    DId0[which(pre$distance1 + pre$distance2 == 0)] <- 1
    DI0.d <- mean(DId0, na.rm = TRUE)
    x1 <- A[-1,]
    x2 <- A[-pre$longueur,]
    dx <- c(x1[, 1] - x2[, 1])
    dy <- c(x1[, 2] - x2[, 2])
    abs.angle.A <- atan2(dy, dx)
    x1 <- B[-1,]
    x2 <- B[-pre$longueur,]
    dx <- c(x1[, 1] - x2[, 1])
    dy <- c(x1[, 2] - x2[, 2])
    abs.angle.B <- atan2(dy, dx)
    DItheta <- cos(abs.angle.A - abs.angle.B)
    DItheta[which(is.na(abs.angle.A) == TRUE &
                    is.na(abs.angle.B) == TRUE)] <- 1
    DItheta[which(is.na(DItheta))] <- 0
    DI.theta <- mean(DItheta, na.rm = TRUE)
    
    # DI
    DI0 <- mean(DItheta * DId0, na.rm = TRUE)
    return(list(
      DI = DI0
      # DI.d = DI0.d,
      # DI.theta = DI.theta
    ))
  }

# DId
DId.function <-
  function(A, B, delta.DI) {
    # adapted from the wildlifeDI package without requiring ltraj objects
    pre <- pre.metrics(A, B)
    DId0 <-
      1 - (abs(pre$distance1 - pre$distance2) / (pre$distance1 + pre$distance2)) ^
      delta.DI
    DId0[which(pre$distance1 + pre$distance2 == 0)] <- 1
    DI0.d <- mean(DId0, na.rm = TRUE)
    # x1 <- A[-1,]
    # x2 <- A[-pre$longueur,]
    # dx <- c(x1[, 1] - x2[, 1])
    # dy <- c(x1[, 2] - x2[, 2])
    # abs.angle.A <- atan2(dy, dx)
    # x1 <- B[-1,]
    # x2 <- B[-pre$longueur,]
    # dx <- c(x1[, 1] - x2[, 1])
    # dy <- c(x1[, 2] - x2[, 2])
    # abs.angle.B <- atan2(dy, dx)
    # DItheta <- cos(abs.angle.A - abs.angle.B)
    # DItheta[which(is.na(abs.angle.A) == TRUE &
    #                 is.na(abs.angle.B) == TRUE)] <- 1
    # DItheta[which(is.na(DItheta))] <- 0
    # DI.theta <- mean(DItheta, na.rm = TRUE)
    # 
    # # DI
    # DI0 <- mean(DItheta * DId0, na.rm = TRUE)
    return(list(
      # DI = DI0,
      DI.d = DI0.d
      # DI.theta = DI.theta
    ))
  }
# DI.function(A=A,B=B,delta.DI)


DItheta.function <-
  function(A, B, delta.DI) {
    # adapted from the wildlifeDI package without requiring ltraj objects
    pre <- pre.metrics(A, B)
    # DId0 <-
    #   1 - (abs(pre$distance1 - pre$distance2) / (pre$distance1 + pre$distance2)) ^
    #   delta.DI
    # DId0[which(pre$distance1 + pre$distance2 == 0)] <- 1
    # DI0.d <- mean(DId0, na.rm = TRUE)
    x1 <- A[-1,]
    x2 <- A[-pre$longueur,]
    dx <- c(x1[, 1] - x2[, 1])
    dy <- c(x1[, 2] - x2[, 2])
    abs.angle.A <- atan2(dy, dx)
    x1 <- B[-1,]
    x2 <- B[-pre$longueur,]
    dx <- c(x1[, 1] - x2[, 1])
    dy <- c(x1[, 2] - x2[, 2])
    abs.angle.B <- atan2(dy, dx)
    DItheta <- cos(abs.angle.A - abs.angle.B)
    DItheta[which(is.na(abs.angle.A) == TRUE &
                    is.na(abs.angle.B) == TRUE)] <- 1
    DItheta[which(is.na(DItheta))] <- 0
    DI.theta <- mean(DItheta, na.rm = TRUE)
    
    # # DI
    # DI0 <- mean(DItheta * DId0, na.rm = TRUE)
    return(list(
      # DI = DI0
      # DI.d = DI0.d,
      DI.theta = DI.theta
    ))
  }



# Measuring CPU time: # compare with Simul06092017.R




DI.function.CPU <-
  function(A,
           B,
           delta.DI,
           longueur,
           distance1,
           distance2) {
    start.time <- proc.time()
    DId0 <-
      1 - (abs(distance1 - distance2) / (distance1 + distance2)) ^ delta.DI
    DId0[which(distance1 + distance2 == 0)] <- 1
    end.time <- proc.time()
    time.taken <- difftime(end.time, start.time, units = "secs")
    DId.time1 <- as.numeric(time.taken)
    
    start.time <- proc.time()
    DI0.d <- mean(DId0, na.rm = TRUE)
    end.time <- proc.time()
    time.taken <- difftime(end.time, start.time, units = "secs")
    DId.time <- as.numeric(time.taken) + DId.time1
    
    
    # DItheta
    start.time <- proc.time()
    x1 <- A[-1,]
    x2 <- A[-longueur,]
    dx <- c(x1[, 1] - x2[, 1])
    dy <- c(x1[, 2] - x2[, 2])
    abs.angle.A <- atan2(dy, dx)
    x1 <- B[-1,]
    x2 <- B[-longueur,]
    dx <- c(x1[, 1] - x2[, 1])
    dy <- c(x1[, 2] - x2[, 2])
    abs.angle.B <- atan2(dy, dx)
    DItheta <- cos(abs.angle.A - abs.angle.B)
    DItheta[which(is.na(abs.angle.A) == TRUE &
                    is.na(abs.angle.B) == TRUE)] <- 1
    DItheta[which(is.na(DItheta))] <- 0
    end.time <- proc.time()
    time.taken <- difftime(end.time, start.time, units = "secs")
    DItheta.time1 <- as.numeric(time.taken)
    
    start.time <- proc.time()
    DI.theta <- mean(DItheta, na.rm = TRUE)
    end.time <- proc.time()
    time.taken <- difftime(end.time, start.time, units = "secs")
    DItheta.time <- as.numeric(time.taken) + DItheta.time1
    
    
    # DI
    start.time <- proc.time()
    DI0 <- mean(DItheta * DId0, na.rm = TRUE)
    end.time <- proc.time()
    time.taken <- difftime(end.time, start.time, units = "secs")
    DI.time <- as.numeric(time.taken) + DItheta.time1 + DId.time1
    
    return(
      list(
        DI = DI0,
        DI.d = DI0.d,
        DI.theta = DI.theta,
        DId.time = DId.time,
        DItheta.time = DItheta.time,
        DI.time = DI.time
      )
    )
  }


sim_dyad <- function(length = 100) {
  A  <-
    cbind(cumsum(rnorm(Nstep, mean = 0, sd = 1)), cumsum(rnorm(Nstep, mean =
                                                                 0, sd = 1)))
  probas <- rbinom(n = Nstep, size = 100, prob = 0.5) / 100
  B  <-
    probas * jitter(A, amount = 0.05) + (1 - probas) * cbind(cumsum((rnorm(
      Nstep, mean = 0, sd = 1
    ))),
    cumsum(rnorm(Nstep, mean =
                   0, sd = 1)))
  return(list(A = A, B = B))
}


# Tuning parameter
Dmax <- 10
delta.DI <- 1
n.cores  <- 4



eval_time <- function(Nstep,type) {
  # start.time <- proc.time()[2]
  sim_data  <- mclapply(1:NSim, function(i_) {
    set.seed(i_)
    sim_dyad(Nstep)
  }, mc.cores = n.cores)
  # end.time <- proc.time()[2]
  time.taken <- system.time(  
    sim_data  <- mclapply(1:NSim, function(i_) {
    set.seed(i_)
    sim_dyad(Nstep)
  }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  # time.taken <- end.time - start.time
  simu.time <- as.numeric(time.taken)
  
  
  # A eviter
  # deb <- proc.time()
  # distances_v1 <- mclapply(sim_data, function(d_)  {diag(as.matrix(pdist.pdist(X = d_$A,Y = d_$B)))}, mc.cores = n.cores)
  # fin <- proc.time()
  # fin-deb
  
  
  
  # start.time <- proc.time()[2]
  # distances <-
  #   mclapply(sim_data, function(d_)  {
  #     sqrt(rowSums(d_$A - d_$B) ^ 2)
  #   }, mc.cores = n.cores)
  # end.time <- proc.time()[2]
  # time.taken <- difftime(end.time, start.time, units = "secs")
  # distances.time <- as.numeric(time.taken)
  
  #### prox1  CPU time
  # start.time <- proc.time()[2]
  time.taken <- system.time(
    distances <-
    mclapply(sim_data, function(d_)  {
      (rowSums(d_$A - d_$B) ^ 2)
    },  mc.cores = n.cores)
  )[type]
  time.taken <- time.taken + system.time(
  prox1 <- mclapply(distances, function(d_) {
    sum(d_ < 1) / length(d_)
  }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()[2]
  # time.taken <- difftime(end.time, start.time, units = "secs")
  prox1.time <- as.numeric(time.taken)
  
  
  
  #### CS  CPU time
  time.taken <- system.time(
  distances <-
    mclapply(sim_data, function(d_)  {
      sqrt(rowSums(d_$A - d_$B) ^ 2)
    },  mc.cores = n.cores),gcFirst = TRUE
  )[type]
  time.taken <- time.taken + system.time(
  do <- mclapply(distances, mean, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  time.taken <- time.taken + system.time(
  distances.de <-
    mclapply(sim_data, function(d_) {
      as.matrix(pdist.pdist(X = d_$A, Y = d_$B))
    }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  time.taken <- time.taken + system.time(
  de <- mclapply(distances.de, mean, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  time.taken <- time.taken + system.time(
  Cs <-
    lapply(1:NSim, function(sim_) {
      (de[[sim_]] - do[[sim_]]) / (de[[sim_]] + do[[sim_]])
    }),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()
  # time.taken <- difftime(end.time, start.time, units = "secs")
  Cs.time <- as.numeric(time.taken)
  
  # 
  # 
  # #### jPPA  CPU time
  # start.time <- proc.time()
  time.taken <- system.time(
  jPPA.index <-
    mclapply(sim_data, function(d_)  {
      jPPA(A = d_$A,
           B = d_$B,
           Dmax = Dmax)
    }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()
  # time.taken <- difftime(end.time, start.time, units = "secs")
  jPPA.time <- as.numeric(time.taken)
  # 
  
  # start.time <- proc.time()
  time.taken <- system.time(
  Lon  <-
    mclapply(sim_data, function(d_) {
      cor(d_$A[, 1], d_$B[, 1])
    }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()
  # time.taken <- difftime(end.time, start.time, units = "secs")
  rlon.time <- as.numeric(time.taken)
  
  
  # start.time <- proc.time()
  time.taken <- system.time(
  Lat  <-
    mclapply(sim_data, function(d_) {
      cor(d_$A[, 2], d_$B[, 2])
    }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()
  # time.taken <- difftime(end.time, start.time, units = "secs")
  rlat.time <- as.numeric(time.taken)
  
  # start.time <- proc.time()
  time.taken <- system.time(
  Lonlat  <-
    mclapply(sim_data, function(d_) {
      cor(rowMeans(d_$A[, ]), rowMeans(d_$B[, ]))
    }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()
  # time.taken <- difftime(end.time, start.time, units = "secs")
  rLongLat.time <- as.numeric(time.taken)
  
  
  # start.time <- proc.time()
  time.taken <- system.time(
  Speed  <- mclapply(sim_data, function(d_) {
    distance1 <- sqrt((diff(d_$A[, 1])) ^ 2 + (diff(d_$A[, 2])) ^ 2)
    distance2 <- sqrt((diff(d_$B[, 1])) ^ 2 + (diff(d_$B[, 2])) ^ 2)
    cor(distance1, distance2, use = 'pairwise.complete.obs')
  }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()
  # time.taken <- difftime(end.time, start.time, units = "secs")
  Speed.time <- as.numeric(time.taken)
  
  # start.time <- proc.time()
  time.taken <- system.time(
  CSEM <- mclapply(sim_data, function(d_) {
    CSEM.function(A = d_$A,
                  B = d_$B,
                  delta = 1)
  }),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()
  # time.taken <- difftime(end.time, start.time, units = "secs")
  CSEM.time <- as.numeric(time.taken)
  
  
  # start.time <- proc.time()
  time.taken <- system.time(
  DI <- mclapply(sim_data, function(d_) {
    DI.function(A = d_$A,
                B = d_$B,
                delta.DI = delta.DI)
  }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()
  # time.taken <- difftime(end.time, start.time, units = "secs")
  DI.time <- as.numeric(time.taken)
  
  time.taken <- system.time(
    DId <- mclapply(sim_data, function(d_) {
      DId.function(A = d_$A,
                  B = d_$B,
                  delta.DI = delta.DI)
    }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()
  # time.taken <- difftime(end.time, start.time, units = "secs")
  DId.time <- as.numeric(time.taken)
  
  time.taken <- system.time(
    DItheta <- mclapply(sim_data, function(d_) {
      DItheta.function(A = d_$A,
                   B = d_$B,
                   delta.DI = delta.DI)
    }, mc.cores = n.cores),gcFirst = TRUE
  )[type]
  # end.time <- proc.time()
  # time.taken <- difftime(end.time, start.time, units = "secs")
  DItheta.time <- as.numeric(time.taken)
  
  res <-
    data.frame(
      simu.time = simu.time,
      # distances.time = distances.time,
      prox1.time = prox1.time,
      Cs.time = Cs.time,
      jPPA.time = jPPA.time,
      rlon.time = rlon.time,
      rlat.time = rlat.time,
      rLongLat.time = rLongLat.time,
      Speed.time = Speed.time,
      CSEM.time = CSEM.time,
      DI.time = DI.time,
      DId.time = DId.time,
      DItheta.time = DItheta.time,
      Nstep = Nstep
    )
  print(res)
  return(res)
}


# first, simulating  Nsim Dyades of size Nstep
NSim  <- 5000
Nstep <- sort(rep(c(100, 500, 1000, 5000, 1e4, 1e5), 10))

time_fin <- lapply(Nstep, function(s_) {
  eval_time(s_,type=3)
})
time_fin_unlist <- do.call('rbind', time_fin)
save(list = c('time_fin'), file = 'timing.RData')

# 
# boxplot(time_fin_unlist$prox1.time~time_fin_unlist$Nstep)
# # boxplot(time_fin_unlist$jPPA.time~time_fin_unlist$Nstep)
# boxplot(time_fin_unlist$Cs.time~time_fin_unlist$Nstep)
# boxplot(time_fin_unlist$CSEM.time~time_fin_unlist$Nstep)
# boxplot(time_fin_unlist$rLongLat.time~time_fin_unlist$Nstep)
# boxplot(time_fin_unlist$DI.time~time_fin_unlist$Nstep)
# 
# # with jPPA
# NSim  <- 5000
# Nstep <- sort(rep(c(100, 500),10))
# 
# time_fin <- lapply(Nstep, function(s_) {
#   eval_time(s_)
# })
# time_fin_unlist <- do.call('rbind', time_fin)
# save(list = c('time_fin'), file = 'timing_with_jppa.RData')
# 
# # with jPPA
# NSim  <- 5000
# Nstep <- sort(rep(c(100, 500),5))
# 
# time_fin <- lapply(Nstep, function(s_) {
#   eval_time(s_,type=3) # type of time to estimate: CPU (2), user (1) and total (3) 
# })
# time_fin_unlist <- do.call('rbind', time_fin)
# save(list = c('time_fin'), file = 'timing_with_jppa_total.RData')
# 
# 
# # first, simulating  Nsim Dyades of size Nstep
# NSim  <- 5000
# Nstep <- sort(rep(c(1000, 5000, 1e4, 1e5), 5))
# 
# time_fin_2 <- lapply(Nstep, function(s_) {
#   eval_time(s_,type=3)
# })
# time_fin_unlist_2 <- do.call('rbind', time_fin_2)
# save(list = c('time_fin_2'), file = 'timing_with_jppa_total_2.RData')
# 
# time_fin_unlist_3 <- rbind(time_fin_unlist,time_fin_unlist_2)
# 
# library(tidyverse)
# time_fin_unlist_3 %>% 
#   group_by(Nstep) %>% 
#   summarise(mean=mean(DItheta.time),min=min(DItheta.time),max=max(DItheta.time))
# 
# boxplot(time_fin_unlist_3$prox1.time~time_fin_unlist_3$Nstep)
# boxplot(time_fin_unlist_3$jPPA.time~time_fin_unlist_3$Nstep)
# boxplot(time_fin_unlist_3$Cs.time~time_fin_unlist_3$Nstep)
# boxplot(time_fin_unlist_3$CSEM.time~time_fin_unlist_3$Nstep)
# boxplot(time_fin_unlist_3$rLongLat.time~time_fin_unlist_3$Nstep)
# boxplot(time_fin_unlist_3$DI.time~time_fin_unlist_3$Nstep)
