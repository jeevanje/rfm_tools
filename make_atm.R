#=====================#
# Function make_atm   #
#                     #
# Make atm file from  #
# given profiles      #
#=====================#

make_atm = function(file,description,z,p,tabs,h2o,co2){
	 # Preliminaries
	 labs  = c("*HGT [km]","*PRE [mb]","*TEM [K]",
	 	    "*H2O [ppmv]","*CO2 [ppmv]")
	 vars  = c("z","p","tabs","h2o","co2")  # SI units
	 facs  = c(1e-3,1e-2,1,1e6,1e6)		# For conversion
	 nvars = length(vars)
	 write_file = function(record){
	 	   write(record,file=file,append=TRUE,sep=", ")
         }
	 nlev = length(z)

	 # Write
	 write_file(description)
	 write_file("! Produced by make_atm in rfm_tools.R")
	 write_file(nlev)
	 for (i in 1:nvars){ 
	    lab  = labs[i]
	    var  = vars[i]
	    fac  = facs[i]
	    write_file(lab)
	    write_file(fac*eval(as.name(var)))
	 }
	 write_file("*END")
}