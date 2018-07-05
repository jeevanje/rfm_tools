# Script to convert ncdf file into rfm atm file
args = commandArgs(trailingOnly=TRUE)

library(ncdf4)
Rtools_dir  = "~/Rtools/"
project_dir = "~/18h2o_feedback/"
source(paste(Rtools_dir,"nc_tools.R",sep=""))
source(paste(Rtools_dir,"thermo_tools.R",sep=""))
source(paste(project_dir,"rfm/make_atm.R",sep=""))

#======#
# Data #
#======#
varlist = c("p","tabs","qv")
period  = 20  # days
ncpath  = args[1]
nc      = nc_open(ncpath)
z	= ncvar_get(nc,"z")
for (var in varlist){
    profile = get_profile(nc,period,var)
    assign(var,profile)
}
nz    = length(z)
n_h2o = qv/m_h2o*m_air
n_co2 = rep(280e-6,times=nz)

#=======#
# Write #
#=======#
sst  = substr(ncpath,41,43)
file = paste(project_dir,"rfm/atm/",sst,"K_atm.txt",sep="")
description = paste("RFM atm file from ",ncpath,sep="")
if (file.exists(file)) {
   file.remove(file)
}
make_atm(file,description,z,p,tabs,n_h2o,n_co2)

levfile = paste(project_dir,"rfm/lev/z",nz,".lev",sep="")
write(1e-3*z,levfile,sep=" ")
