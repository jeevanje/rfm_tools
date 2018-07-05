args = commandArgs(trailingOnly=TRUE)
source("../Rtools/thermo_tools.R")
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
rfmdir     = "~/17rad_cooling2/rfm"
casedir    = paste(rfmdir,case,sep="/")
drvfile    = paste(casedir,"/rfm.drv",sep="")
atmfile    = scan(drvfile,skip=9,sep="/",nmax=3,what="raw")[3]
atmpath    = paste(rfmdir,"/atm/",atmfile,sep="")

# z,p,tabs,qh2o,qco2
print("Reading thermodynamic profiles")
nlev       = as.numeric(scan(atmpath,skip=2,nmax=1)[1])
nlines     = ceiling(nlev/5)
nlines_tab = 1+ceiling((nlev-3)/5)
z	   = 1e3*read_atm(atmpath,skip=4+0*(1+nlines),nlev=nlev)  # m
p	   = 1e2*read_atm(atmpath,skip=4+1*(1+nlines),nlev=nlev)  # Pa
tabs       = 1e0*read_atm(atmpath,skip=4+2*(1+nlines),nlev=nlev)  # K  
q_h2o      = m_h2o/m_air*1e-6*read_atm(atmpath,skip=4+3*(1+nlines),nlev=nlev) 
q_co2      = m_co2/m_air*1e-6*read_atm(atmpath,skip=4+4*(1+nlines),nlev=nlev)   
# q in kg/kg

# k
kdata      = scan(drvfile,skip=5,nmax=3,what="raw")
k1	   = 1e2*as.numeric(kdata[1])  # m^-1
k_nk       = 1e2*as.numeric(kdata[2])  # m^-1
dk	   = 1e2*as.numeric(kdata[3])  # m^-1
k	   = seq(from = k1, to = k_nk,by=dk)
nk	   = length(k)

# tab
print("Reading tab data")
tabfile    = paste(casedir,"/tab/tab.asc",sep="")
#tab_data   = read.table(tabfile,skip=10)
#sigma      = exp(tab_data[[2]])/1000/N_avo 
sigma      = array(dim=c(nk,nlev))
for (i in 1:nk){
    #print(paste("k=",k[i],sep=""))
    skip      = 5+3*nlines+2+(i-1)*nlines_tab
    lnsigma_i = scan(tabfile,skip=skip,nmax=nlev+1,quiet=TRUE)[2:(nlev+1)]
    sigma[i,] = exp(lnsigma_i)/1000/N_avo   # m^2/molecule
}

# 2d fields
for (var2d in c("coo","opt","flx")){
    print(paste("Reading ",var2d," data",sep=""))
    field      = array(dim=c(nk,nlev))
    datadir    = paste(rfmdir,case,var2d,sep="/")
    for (m in 1:nlev){
       zval      = z[m]    # m
       if (zval < 1e5){
	  zstring = formatC(zval,format="d",width=5,flag="0")
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
pdim = ncdim_def("p","Pa",p,longname="pressure")

vars_p = list()

vars_p[['z']]	          = list()
vars_p[['z']]$longname    = "Height"
vars_p[['z']]$units       = "m"
vars_p[['z']]$data        = z

vars_p[['tabs']]          = list()
vars_p[['tabs']]$longname = "Absolute Temperature"
vars_p[['tabs']]$units    = "K"
vars_p[['tabs']]$data     = tabs

vars_p[['q_h2o']]          = list()
vars_p[['q_h2o']]$longname = "H2O specific concentration"
vars_p[['q_h2o']]$units    = "(kg h2o)/(kg moist air)"
vars_p[['q_h2o']]$data     = q_h2o

vars_p[['q_co2']]          = list()
vars_p[['q_co2']]$longname = "CO2 specific concentration"
vars_p[['q_co2']]$units    = "(kg co2)/(kg air)"
vars_p[['q_co2']]$data     = q_co2

vars_kp = list()

vars_kp[['sigma']]          = list()
vars_kp[['sigma']]$longname = "Absorption Coefficient"
vars_kp[['sigma']]$units    = "m^2/molec."
vars_kp[['sigma']]$data     = sigma
    
vars_kp[['coo']]           = list()
vars_kp[['coo']]$longname  = "Radiative Cooling"
vars_kp[['coo']]$units     = "K/day/cm^-1"
vars_kp[['coo']]$data      = coo

vars_kp[['opt']]           = list()
vars_kp[['opt']]$longname  = "Optical Depth"
vars_kp[['opt']]$units     = "Unitless"
vars_kp[['opt']]$data      = opt

vars_kp[['flx']]           = list()
vars_kp[['flx']]$longname  = "Net Upward Flux"
vars_kp[['flx']]$units     = "W/m^2/cm^-1"
vars_kp[['flx']]$data      = flx

# Make vardefs
vardef_p <- list()
for (name in names(vars_p)) {
  vardef_p[[name]] <- ncvar_def(name,vars_p[[name]]$units,pdim,
     missval=-999,longname=vars_p[[name]]$longname)
}

vardef_kp <- list()
for (name in names(vars_kp)) {
  vardef_kp[[name]] <- ncvar_def(name,vars_kp[[name]]$units,list(kdim,pdim),
     missval=-999,longname=vars_kp[[name]]$longname)
}

# Create the NetCDF file with vardef
filename = paste(casedir,"/",case,".nc",sep="")
nc 	 = nc_create(filename,c(vardef_p,vardef_kp))

# Fill the NetCDF with data
vars = c(vars_p,vars_kp)
for (name in names(vars)){
   ncvar_put(nc,name,vars[[name]]$data)
}

nc_close(nc)

   

