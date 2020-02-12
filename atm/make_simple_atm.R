#args = as.numeric(commandArgs(trailingOnly=TRUE))

source("~/Rtools/thermo_tools.R")
source("~/Rtools/rfm_tools.R")
source("../make_rfm_grd.R")
source("../make_atm.R")

# Set these params only!
Ts        = args[1]  # K
rh        = args[2]
gamma     = args[3]  # K/km or "m" 
gamma_st  = args[4]  # K/km, nonnegative!
co2_ppmv  = args[5]
Ttp	  = 200   # K

rfmdir	    = "~/19feedbacks/rfm/"
atmfile     = paste(rfmdir,"atm/Ts",Ts,"_rh",rh,"_gamma_st",gamma_st,"_q",co2_ppmv,".atm",sep="")
description = paste("! Simple atmosphere with co2_ppmv = ",co2_ppmv,sep="")

if (file.exists(atmfile)){
   print(paste("File ",atmfile," exists, removing ...",sep=""))
   file.remove(atmfile)
}


ptemp = 1e2*seq(1000,10,by=-10)
Htemp = 8e3   # crude scale height, km
z     = -Htemp*log(ptemp/ps)
z     = round(z,digits=-1)
nz    = length(z)

# tabs, qv
if (gamma != "m"){
   gamma = 1e-3*gamma  # K/m 
   tabs  = Ts - gamma*z
   p     = ps*(tabs/Ts)^(g/gamma/Rd)
} else if (gamma == "m"){
   tabs    = numeric(nz)
   p       = numeric(nz)
   tabs[1] = Ts
   p[1]	   = ps
   for (k in 2:nz){
      tabs[k] = tabs[k-1] - gamma_m(tabs[k-1],p[k-1])*(diff(z)[k-1])
      p[k]    = p[k-1]    - g*p[k-1]/Rd/tabs[k-1]*(diff(z)[k-1])
   }
}   
qv = rh*qsat(tabs,p)

# add strat
k_tp  = max(which(tabs>Ttp))+1
z_tp  = z[k_tp]
p_tp  = p[k_tp]
qv_tp = qv[k_tp]
tabs[(k_tp):nz] = Ttp + 1e-3*gamma_st*(z[k_tp:nz]-z_tp)
qv[(k_tp):nz]   = qv_tp
p[(k_tp):nz]    = p_tp*exp(-g/Rd/Ttp*(z[(k_tp):nz]-z_tp))
n_h20 = qv/m_h2o*m_air
n_co2 = rep(co2_ppmv*1e-6,times=nz)

make_atm(atmfile,description,z,p,tabs,n_h20,n_co2)

levfile     = paste(rfmdir,"lev/z",nz,".lev",sep="")
write(1e-3*z,levfile,sep=" ")
