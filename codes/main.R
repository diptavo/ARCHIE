
estimate_uv <- function(W,pen.u,pen.v,init_u, init_v,prec = 1e-03,max.iter = 200, verbose = TRUE){
  
  dim.u <- dim(W)[1]; dim.v <- dim(W)[2]
  if(is.null(init_u))
    init_u <- normalize(rnorm(dim.u))
  if(is.null(init_v))
    init_v <- normalize(rnorm(dim.v))
  u0 <- init_u; v0 <- init_v; diffs <- 2;iter <- 1
  while((diffs > prec) && (iter < max.iter)){
    if(verbose)
      print(paste("Iteration ",iter))
    u1 <- as.numeric(normalize(soft(W%*%v0,para = pen.u)));
    v1 <- as.numeric(normalize(soft(t(W)%*%u0,para = pen.v)));
    q <- as.numeric(t(u1)%*%W%*%v1)
    if(sum(is.nan(u1)>0) ||sum(is.nan(v1) > 0)){
      return(list("u" = rep(0,dim.u), "v"= rep(0,dim.v),"q" = 0));
      exit(0)
    }
    diffs <- sum((u0-u1)^2)/dim.u + sum((v0-v1)^2)/dim.v
    if(verbose)
      print(paste0("the difference is: ",diffs))
    u0 <- u1; v0 <- v1;
    iter <- iter+1
  }
  return(list("u" = u0, "v"= v0,"q" = q))
  
}


estimate_uv_sp <- function(W,pen.u,pen.v,pos.u = NULL, pos.v = NULL,init_u, init_v,prec = 1e-03,max.iter = 200, verbose = TRUE){
  
  dim.u <- dim(W)[1]; dim.v <- dim(W)[2]
  if(is.null(init_u))
    init_u <- normalize(rnorm(dim.u))
  if(is.null(init_v))
    init_v <- normalize(rnorm(dim.v))
  u0 <- init_u; v0 <- init_v; diffs <- 2;iter <- 1
  while((diffs > prec) && (iter < max.iter)){
    if(verbose)
      print(paste("Iteration ",iter))
    u1 <- as.numeric((soft(W%*%v0,para = pen.u)));u1[pos.u] <- 0; u1 <- normalize(u1);
    v1 <- as.numeric(soft(t(W)%*%u0,para = pen.v)); v1[pos.v] <- 0; v1 <- normalize(v1);
    q <- as.numeric(t(u1)%*%W%*%v1)
    diffs <- sum((u0-u1)^2)/dim.u + sum((v0-v1)^2)/dim.v
    if(verbose)
      print(paste0("the difference is: ",diffs))
    u0 <- u1; v0 <- v1;
    iter <- iter+1
  }
  return(list("u" = u0, "v"= v0,"q" = q))
  
}



intersect_min_12 <- function(W,u10,v10,u20,v20,dec.u = 0, dec.v = 0,prec.grid.u =1000, prec.grid.v =1000, max.iter = 200,method = c("both"),verbose = TRUE){
  
  dec.u = 0; dec.v = 0; prec.grid.u =100; prec.grid.v =100; max.iter = 2000; method = c("both"); verbose = TRUE;
  
  grid.u <- seq(0,1,by = 1/prec.grid.u); 
  grid.v <- seq(0,1,by = 1/prec.grid.v);
  max.iter <- min(prec.grid.v,prec.grid.u) - 1
  iter <- 1; conv_u <- length(u10); conv_v <- length(v10);move.u <- 1;move.v <- 1
  if(method == "both"){
    while((conv_u > dec.u || conv_v > dec.v) && (iter < max.iter )){
      if(iter > max.iter){
        return(list("comp1" = A1, "comp2" = A2))
        if(verbose)
          print(paste0("Iterations reached maximum limit before convergence!!"))
        exit(0)
      }
      q <- 0
      if(verbose)
        print(paste0("Iteration ",iter))
      if(move.u == 1)
        p.u <- grid.u[iter];
      if(move.v == 1)
        p.v <- grid.v[iter];
      if(verbose)
        print(paste0("penalty.u = ",p.u,"; penalty.v = ",p.v))
      W_n <- W - q*(u10%*%t(v10))
      A1 <- estimate_uv_sp(W_n,pen.u = p.u,pen.v = p.v,init_u = u10,init_v = v10,verbose = F)
      u1 <- A1$u; v1 <- A1$v; q1 <- A1$q; 
      W_n <- W - q1*(u1%*%t(v1));
      A2 <- estimate_uv_sp(W_n,pen.u = p.u,pen.v = p.v, init_u = u20,init_v = v20,verbose = F)
      u2 <- A2$u; v2 <- A2$v; q2 <- A2$q;
      conv_u <- length(intersect(which(abs(u1) > 0), which(abs(u2) > 0)))
      conv_v <- length(intersect(which(abs(v1) > 0), which(abs(v2) > 0)))
      move.u <- as.numeric(conv_u > dec.u); move.v <- as.numeric(conv_v > dec.v)
      iter <- iter+1
    }
    print("Iterations ended since both u and v converged")
  }
  
  
  if(method == "first"){
    while(conv_u > dec.u && conv_v > dec.v){
      if(verbose)
        print(paste0("Iteration ",iter))
      if(move.u == 1)
        p.u <- grid.u[iter];
      if(move.v == 1)
        p.v <- grid.v[iter];
      if(verbose)
        print(paste0("penalty.u = ",p.u,"; penalty.v = ",p.v))
      W_n <- W - q*(u10%*%t(v10))
      A1 <- estimate_uv(W_n,pen.u = p.u,pen.v = p.v,init_u = u10,init_v = v10,verbose = F)
      u1 <- A1$u; v1 <- A1$v; q1 <- A1$q;
      W_n <- W - q1*(u1%*%t(v1))
      A2 <- estimate_uv(W_n,pen.u = p.u,pen.v = p.v,init_u = u20,init_v = v20,verbose = F)
      u2 <- A2$u; v2 <- A2$v; q2 <- A2$q;
      conv_u <- length(intersect(which(abs(u1) > 0), which(abs(u2) > 0)))
      conv_v <- length(intersect(which(abs(v1) > 0), which(abs(v2) > 0)))
      move.u <- as.numeric(conv_u > dec.u); move.v <- as.numeric(conv_v > dec.v)
      iter <- iter+1
    }
    print("Iterations ended since at least on of u or v has converged")
  }
  return(list("comp1" = A1, "comp2" = A2))
  
}


intersect_min_3on <- function(W,u,v,u30,v30,dec.u = 0, dec.v = 0,prec.grid.u =1000, prec.grid.v =1000,method = c("first","both"),verbose = TRUE,max.iter = 200){
  
  
  #  dec.u = 0; dec.v = 0; prec.grid.u =100; prec.grid.v =100; max.iter = 2000; method = c("both"); verbose = TRUE;
  
  grid.u <- seq(0,1,by = 1/prec.grid.u); 
  grid.v <- seq(0,1,by = 1/prec.grid.v); 
  iter <- 1; conv_u <- length(u30); conv_v <- length(v30);move.u <- 1;move.v <- 1
  max.iter <- min(prec.grid.v,prec.grid.u) - 1
  W_n <- W
  if(method == "both"){
    while((conv_u > dec.u || conv_v > dec.v) && (iter < max.iter )){
      q <- 0
      if(verbose)
        print(paste0("Iteration ",iter))
      if(move.u == 1)
        p.u <- grid.u[iter];
      if(move.v == 1)
        p.v <- grid.v[iter];
      if(verbose)
        print(paste0("penalty.u = ",p.u,"; penalty.v = ",p.v))
      
      u1 <- u; v1 <- v; 
      
      A2 <- estimate_uv(W_n,pen.u = p.u,pen.v = p.v,init_u = u30,init_v = v30,verbose = F)
      u2 <- A2$u; v2 <- A2$v; q2 <- A2$q;
      conv_u <- length(intersect(which(abs(u1) > 0), which(abs(u2) > 0)))
      conv_v <- length(intersect(which(abs(v1) > 0), which(abs(v2) > 0)))
      move.u <- as.numeric(conv_u > dec.u); move.v <- as.numeric(conv_v > dec.v)
      iter <- iter+1
      W_n <- W - q*(u1%*%t(v1))
    }
    print("Iterations ended since both u and v converged")
  }
  
  if(method == "first"){
    while((conv_u > dec.u && conv_v > dec.v) && (iter < max.iter )){
      if(verbose)
        print(paste0("Iteration ",iter))
      if(move.u == 1)
        p.u <- grid.u[iter];
      if(move.v == 1)
        p.v <- grid.v[iter];
      if(verbose)
        print(paste0("penalty.u = ",p.u,"; penalty.v = ",p.v))
      
      u1 <- u; v1 <- v; 
      W_n <- W - q*(u1%*%t(v1))
      A2 <- estimate_uv(W_n,pen.u = p.u,pen.v = p.v,init_u = u20,init_v = v20,verbose = F)
      u2 <- A2$u; v2 <- A2$v; q2 <- A2$q;
      conv_u <- length(intersect(which(abs(u1) > 0), which(abs(u2) > 0)))
      conv_v <- length(intersect(which(abs(v1) > 0), which(abs(v2) > 0)))
      move.u <- as.numeric(conv_u > dec.u); move.v <- as.numeric(conv_v > dec.v)
      iter <- iter+1
    }
    print("Iterations ended since at least on of u or v has converged")
  }
  
  return(list("comp" = A2))
}

archie_work <- function(Sigma_GE,Sigma_GG,Sigma_EE,K = NULL, dec.u = 0, dec.v = 0,prec.grid.u =1000, prec.grid.v =1000,method = c("both"),verbose = TRUE){
  
  W <- Sigma_GG%*%Sigma_GE%*%Sigma_EE;
  dim.u <- dim(W)[1]; dim.v <- dim(W)[2]
  if(verbose)
    print(paste0("summary statistics across ",dim.u," variants and ",dim.v," genes are provided..."))

    svd.g <- eigen(W%*%t(W));
  if(is.null(K)){
    d1 <- -diff(svd.g$values); K <- which.max(d1);
    print(paste0("Extracting top ",K," components"));
  }
  
  u10 <- normalize(svd.g$vectors[,1]); v10 <- normalize(t(W)%*%svd.g$vectors[,1])
  u20 <- normalize(svd.g$vectors[,2]); v20 <- normalize(t(W)%*%svd.g$vectors[,2])
  
  obj <- intersect_min_12(W,u10,v10,u20,v20,method = "both",prec.grid.u = 100,prec.grid.v = 100,max.iter = 1000,verbose = verbose)
  df.u <- cbind(obj$comp1$u,obj$comp2$u); df.v <- cbind(obj$comp1$v,obj$comp2$v);
  qs <- c(obj$comp1$q,obj$comp2$q);
  es <- eigen(W%*%t(W));ev <- es$values; uvecs <- es$vectors; uvecs <- apply(uvecs,2,normalize); vvecs <-t(W)%*%uvecs; vvecs <- apply(vvecs,2,normalize)
  W_n <- W - ev[1]*(uvecs[,1]%*%t(vvecs[,1])) - ev[2]*(uvecs[,2]%*%t(vvecs[,2]))
  qs[1] <- qs[1]^2/sum(svd.g$values^2); W_2 <-  W - ev[1]*(uvecs[,1]%*%t(vvecs[,1])); svd.2 <- eigen(W_2%*%t(W_2));qs[2] <- qs[2]^2/sum(svd.2$values^2)
  if(K>2){
    for(i in 3:K){
      
      init.u <- normalize(uvecs[,i]); init.v <- normalize(vvecs[,i])
      obj.1 <- intersect_min_3on(W_n,u30 = init.u,v30 = init.v,u = df.u[,(i-1)],v = df.v[,(i-1)],dec.u = 0,dec.v = 0,method = "both",prec.grid.u = 100,prec.grid.v = 100,verbose = verbose)
      df.u <- cbind(df.u,obj.1$comp$u); df.u <- force_sp(df.u); df.v <- cbind(df.v,obj.1$comp$v);df.v <- force_sp(df.v)
      W_n <- W_n - ev[i]*df.u[,i]%*%t(df.v[,i]); svd.n <- eigen(W_n%*%t(W_n));
      qs <- c(qs,obj.1$comp$q/sum(svd.n$values^2))  
    }
  }
  df.u.sel <- df.u; df.u.sel[df.u.sel !=0 ] <- 1; 
  df.v.sel <- df.v; df.v.sel[df.v.sel !=0 ] <- 1; 
  # for(i in 1:K){
  #   if(qs[i] < 0.00001)
  #     qs[i:K] <- 0; df.u.sel[,i:K] <- 0; df.v.sel[,i:K] <- 0
  # }
  return(list("us" = df.u.sel,"vs" = df.v.sel, "qs" = qs, "K" = K))
  
  
  
}





