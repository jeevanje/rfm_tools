args = as.numeric(commandArgs(trailingOnly=TRUE))

source("~/Rtools/thermo_tools.R")
source("~/Rtools/rfm_tools.R")
source("../make_rfm_grd.R")
source("../make_atm.R")

Ts	  = args[1]  # K
rh        = args[2]
gamma     = args[3]  # K/km !
co2_ppmv  = args[4]

rfmdir	    = "~/18co2/rfm/"
atmfile     = paste(rfmdir,"atm/Ts",Ts,"_rh",rh,"_q",co2_ppmv,".atm",sep="")
if (file.exists(atmfile)){
   print(paste("File ",atmfile," exists, removing ...",sep=""))
   file.remove(atmfile)
}

Ttp   = 200   # K
ps    = 1e5   # Pa
gamma = 1e-3*gamma  # K/m 

z = c( 0.0,  0.5,  1.0,  1.5,  2.0,  
       2.5,  3.0,  3.5,  4.0,  4.5,
       5.0,  5.5,  6.0,  6.5,  7.0,  
       7.5,  8.0,  8.5,  9.0,  9.5,
      10.0, 10.5, 11.0, 11.5, 12.0, 
      12.5, 13.0, 13.5, 14.0, 14.5, 
      15.0, 16.0, 17.0,
      18.0, 19.0, 20.0, 21.0, 22.0,
      23.0, 24.0, 25.0, 27.5, 30.0,
      32.5, 35.0, 37.5, 40.0, 42.5,
      45.0, 47.5, 50.0)*1e3   # m
z     = make_rfm_grd(c(50,500,2500),5e4,c(1e3,18e3),8)
z     = seq(0,3.5e4,by=100)  # m
ptemp = 1e2*seq(1000,10,by=-10)
Htemp = 8e3   # crude scale height, km
z     = -Htemp*log(ptemp/ps)
z     = round(z,digits=-1)
tabs  = Ts - gamma*z
p     = ps*(tabs/Ts)^(g/gamma/Rd)
qv    = rh*qsat(tabs,p)
nz    = length(z)

# add strat
k_tp  = max(which(tabs>Ttp))+1
z_tp  = z[k_tp]
p_tp  = p[k_tp]
qv_tp = qv[k_tp]
tabs[(k_tp):nz] = Ttp
qv[(k_tp):nz]   = qv_tp
p[(k_tp):nz]    = p_tp*exp(-g/Rd/Ttp*(z[(k_tp):nz]-z_tp))
n_h20 = qv/m_h2o*m_air
n_co2 = rep(co2_ppmv*1e-6,times=nz)

description = paste("! Simple atmosphere with co2_ppmv = ",co2_ppmv,sep="")
make_atm(atmfile,description,z,p,tabs,n_h20,n_co2)

levfile     = paste(rfmdir,"lev/z",nz,".lev",sep="")
write(1e-3*z,levfile,sep=" ")
