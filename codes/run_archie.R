library(optparse)
library(data.table)

option_list <- list(
  make_option(c("-i","--sigmage"),type="character", default = "Sigma_GE.txt",
              help = "file containing sigma_ge (cross-covariance) values"),
  make_option(c("-g","--sigmagg"),type="character", default = "Sigma_GG.txt",
              help = "file containing sigma_gg (snp-covariance/LD) values"),
  make_option(c("-e","--sigmaee"),type="character", default = "Sigma_EE.txt",
              help = "file containing sigma_ee (gene-covariance/co-expression) values"),
  make_option(c("-k","--K"),type="character", default = "opt",
              help = "number of components to be extracted"),
  make_option(c("-d","--dir"),type="character", default ="~/archie/",
              help = "directory fore all the base code files"),
  make_option(c("-u","--precu"),type="integer", default =100,
              help = "grid search parameter for snp-component"),
  make_option(c("-v","--precv"),type="integer", default =100,
              help = "grid search parameter for gene-component"),
  make_option(c("-r","--decu"),type="integer", default =0,
              help = "maximum intersect for snp-component"),
  make_option(c("-s","--decv"),type="integer", default =0,
              help = "maximum intersect for gene-component"),
  make_option(c("-q", "--verbose"), type="character", default=FALSE,
              help="Print extra output [default]"),
  make_option(c("-o","--out"),type="character",default = "out",
              help = "prefix of output rds file")
  
)

opt <- parse_args(OptionParser(option_list=option_list))
f.i <- normalizePath(as.character(opt$sigmage));
f.g <- normalizePath(as.character(opt$sigmagg));f.e <- normalizePath(as.character(opt$sigmaee)); out <- as.character(opt$out)
K <- as.character(opt$K); dir <- as.character(opt$dir); dec.u <- as.numeric(opt$decu); dec.v <- as.numeric(opt$dec.v);
prec.grid.u <- as.numeric(opt$precu);prec.grid.v <- as.numeric(opt$precv); verbose = as.character(opt$verbose) 

if(K == "opt"){
  K = NULL}else{K = as.integer(K)}

source(paste0(normalizePath(dir),"/helper.R"));
source(paste0(normalizePath(dir),"/main.R"));
print("reading files....")
Sigma_GG <- as.matrix(read.table(f.g));
Sigma_EE <- as.matrix(read.table(f.e));
Sigma_GE <- as.matrix(read.table(f.i));
outf <- paste0(out,".rds");

A1 <- archie_work(Sigma_GE,Sigma_GG,Sigma_EE,K = K,dec.u = dec.u,dec.v = dec.v,prec.grid.u = prec.grid.u,prec.grid.v = prec.grid.v,method = "both",verbose = verbose)
saveRDS(object = A1,file = outf);

