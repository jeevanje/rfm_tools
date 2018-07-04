args = commandArgs(trailingOnly=TRUE)

library(ncdf4)

#=======================#
# Function read_atm     #
#                       #
# Reads .atm RFM input  #
# and outputs single    #
# atmospheric field     #
#=======================#

read_atm = function(atmpath,skip,nlev){
	var = scan(atmpath,skip=skip,sep=",",nlines=ceiling(nlev/5),
			strip.white=TRUE,skipNul=TRUE,quiet=TRUE) 
        var = var[!is.na(var)]
        dim(var) <- nlev
        return(var)
}

#===========#
# Get data  #   
#===========#

case       = args[1]
rfmdir     = "~/18h2o_feedback/rfm"
casedir    = paste(rfmdir,case,sep="/")
drvfile    = paste(casedir,"/rfm.drv",sep="")
atmfile    = scan(drvfile,skip=9,sep="/",nmax=3,what="raw")[3]
atmpath    = paste(rfmdir,"/atm/",atmfile,sep="")

# z,p,tabs
nlev       = as.numeric(scan(atmpath,skip=2,nmax=1)[1])
nlines     = ceiling(nlev/5)
z	   = 1e3*read_atm(atmpath,skip=4,nlev=nlev)   #m
p	   = 1e2*read_atm(atmpath,skip=4+(1+nlines),nlev=nlev)  #Pa
tabs       = read_atm(atmpath,skip=4+2*(1+nlines),nlev=nlev)    

# k
kdata      = scan(drvfile,skip=5,nmax=3,what="raw")
k1	   = 1e2*as.numeric(kdata[1])  # m^-1
k_nk       = 1e2*as.numeric(kdata[2])  # m^-1
dk	   = 1e2*as.numeric(kdata[3])  # m^-1
k	   = seq(from = k1, to = k_nk,by=dk)
nk	   = length(k)

# 2d fields
for (var2d in c("coo","opt","flx")){
    field      = array(dim=c(nk,nlev))
    datadir    = paste(rfmdir,case,var2d,sep="/")
    for (m in 1:nlev){
       zval      = z[m]    # m
       if (zval < 1e5){
	  zstring = formatC(round(zval,digits=0),format="d",width=5,flag="0")
       } else {
	  zstring = formatC(zval/1000,format="d",width=5,flag="0")
       }
       file   = paste(datadir,"/",var2d,"_",zstring,".asc",sep="")
       field[ ,m] = scan(file,skip=4,quiet=TRUE)
    }
    assign(var2d,field)
}

#===========#
# make ncdf #   
#===========#

kdim = ncdim_def("k","m^-1",k,longname="wavenumber")
zdim = ncdim_def("z","m",z,longname="height")
    
vars1d = list()

vars1d[['p']]	          = list()
vars1d[['p']]$longname    = "Pressure"
vars1d[['p']]$units       = "Pa"
vars1d[['p']]$data        = p

vars1d[['tabs']]          = list()
vars1d[['tabs']]$longname = "Absolute Temperature"
vars1d[['tabs']]$units    = "K"
vars1d[['tabs']]$data     = tabs

vars2d = list()

vars2d[['coo']]           = list()
vars2d[['coo']]$longname  = "Radiative Cooling"
vars2d[['coo']]$units     = "K/day/cm^-1"
vars2d[['coo']]$data      = coo

vars2d[['opt']]           = list()
vars2d[['opt']]$longname  = "Optical Depth"
vars2d[['opt']]$units     = "Unitless"
vars2d[['opt']]$data      = opt

vars2d[['flx']]           = list()
vars2d[['flx']]$longname  = "Net Upward Flux"
vars2d[['flx']]$units     = "W/m^2/cm^-1"
vars2d[['flx']]$data      = flx

# Make vardefs
vardef1d <- list()
for (name in names(vars1d)) {
  vardef1d[[name]] <- ncvar_def(name,vars1d[[name]]$units,zdim,
     missval=-999,longname=vars1d[[name]]$longname)
}

vardef2d <- list()
for (name in names(vars2d)) {
  vardef2d[[name]] <- ncvar_def(name,vars2d[[name]]$units,list(kdim,zdim),
     missval=-999,longname=vars2d[[name]]$longname)
}

# Create the NetCDF file with vardef
filename = paste(casedir,"/",case,".nc",sep="")
nc 	 = nc_create(filename,c(vardef1d,vardef2d))

# Fill the NetCDF with data
vars = c(vars1d,vars2d)
for (name in names(vars)){
   ncvar_put(nc,name,vars[[name]]$data)
}

nc_close(nc)

   

