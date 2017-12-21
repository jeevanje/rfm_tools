case       = "test1" 
rfmdir     = "~/rad_cooling2/rfm/"
casedir    = paste(rfmdir,case,"/",sep="")
zfile      = paste(rfmdir,"lev/z.lev",sep="")
z	   = scan(zfile)   #km
nz         = length(z)

radfiles   = list.files(casedir,pattern="rad_")
nk         = length(scan(paste(casedir,radfiles[1],sep=""),skip=4))

radiance = array(dim=c(nk,nz))

for (k in 1:nz){
   zval      = z[k]    # km
   if (zval < 100){
      zstring = formatC(zval*1000,format="d",width=5,flag="0")
   } else {
      zstring = formatC(zval,format="d",width=5,flag="0")
   }
   file   = paste(casedir,"rad_",zstring,".asc",sep="")
   radiance[ ,k] = scan(file,skip=4)
}
