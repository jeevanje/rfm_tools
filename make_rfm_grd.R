# Function make_grd
make_rfm_grd <- function(dz_vals,max_height,transition_heights,transition_steps){

   B <- transition_steps

   k <- 1
   grd <- dz_vals[1]/2
   C <- dz_vals[1]
   D <- grd - C*k

   for (i in 1:length(transition_heights)) {

      A <- (pi/2/B)*(dz_vals[i+1]-dz_vals[i])
      k0 <- (transition_heights[i] - D + A*B^2*(1/pi^2 - 1/(2*pi))) / C - B/2

      if ((k+1) <= ceiling(k0)-1) {
         k_range <- (k+1):(ceiling(k0)-1)
         grd <- c(grd,C*k_range + D)
      }
      k <- ceiling(k0)-1
   
      k_range <- (k+1):(ceiling(k0)+B-1)
      grd <- c(grd,C*k_range + D - A*B^2/(pi^2)*sin(pi*(k_range-k0)/B) + A*B/pi*(k_range-k0))
      k <- ceiling(k0)+B-1
   
      C <- dz_vals[i+1]
      D <- max(grd) - C*k
 
   }

   k0 <- (max_height - D) / C

   if ((k+1) <= ceiling(k0)) {
      k_range <- (k+1):(ceiling(k0))
      grd <- c(grd,C*k_range + D)
   }
   k <- ceiling(k0)

   grd<-round(grd-grd[1], digits=1)
   return(grd)

}
