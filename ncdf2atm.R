# Script to convert ncdf file into rfm atm file
args = commandArgs(trailingOnly=TRUE)

library(ncdf4)
Rtools_dir  = "~/Rtools/"
project_dir = "~/18co2/"
source(paste(Rtools_dir,"nc_tools.R",sep=""))
source(paste(Rtools_dir,"thermo_tools.R",sep=""))
source(paste(project_dir,"rfm/make_atm.R",sep=""))

#======#
# Data #
#======#
ncpath    = args[1]
ncname    = args[2]
model     = args[3]
nco2_ppmv = as.numeric(args[4])
nc        = nc_open(ncpath)

if (model=="dam"){
   varlist = c("p","tabs","qv")
   period  = 20  # days
   z	= ncvar_get(nc,"z")
   for (var in varlist){
        profile = get_profile(nc,period,var)
    	assign(var,profile)
   }
   np    = length(z)
   n_co2 = rep(nco2_ppmv*1e-6,times=nz)
} else if (model=="gfdl"){
   p     = ncvar_get(nc,"pressm")
   np    = length(p)
   dz    = ncvar_get(nc,"deltaz")
   z     = numeric(np)
   z[np] = dz[np]/2
   for (k in (np-1):1){
       z[k] = z[k+1] + (dz[k]+dz[k+1]/2)
   }
   tabs  = ncvar_get(nc,"temp")
   r     = ncvar_get(nc,"rh2o")
   qv    = r/(1+r)   
   n_co2 = rep(nco2_ppmv*1e-6,times=np)
   n_h2o = qv/m_h2o*m_air
   varlist = c("z","p","tabs","qv","n_h2o","n_co2")
   for (var in varlist){
        vals  = rev(eval(as.name(var)))
    	assign(var,vals)
   }
}



#=======#
# Write #
#=======#
file = paste(project_dir,"rfm/atm/",ncname,"_",nco2_ppmv,".atm",sep="")
description = paste("RFM atm file from ",ncpath,sep="")
if (file.exists(file)) {
   file.remove(file)
}
make_atm(file,description,z,p,tabs,n_h2o,n_co2)

levfile = paste(project_dir,"rfm/lev/z",np,".lev",sep="")
write(1e-3*z,levfile,sep=" ")
