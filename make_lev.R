# Script to make levfile from .atm
args = commandArgs(trailingOnly=TRUE)
source("./make_atm.R")

atmpath = args[1]
nlev    = as.numeric(scan(atmpath,skip=2,nmax=1)[1])
nlines  = ceiling(nlev/5)
z       = read_atm(atmpath,skip=4+0*(1+nlines),nlev=nlev)  # km

levfile = paste("./lev/z",nlev,".lev",sep="")
write(z,levfile,sep=" ")
