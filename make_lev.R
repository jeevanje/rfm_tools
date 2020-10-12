# Script to make levfile from .atm
# input atm name
args = commandArgs(trailingOnly=TRUE)
source("./make_atm.R")

atm     = args[1]
atmpath = paste("atm/",atm,".atm",sep="")
nlev    = as.numeric(scan(atmpath,skip=2,nmax=1)[1])
nlines  = ceiling(nlev/5)
z       = read_atm(atmpath,skip=4+0*(1+nlines),nlev=nlev)  # km

levfile = paste("./lev/",atm,".lev",sep="")
write(z,levfile,sep=" ")
