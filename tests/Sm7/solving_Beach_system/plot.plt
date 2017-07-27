# The plot command for visualizing the spectra for different Least Square routines applied for the Beach matrix.
# The settings are documented in the pade.par file. 
# For the gelsd simulation using fortran, the code was recompiled so it would use the ZGELSD routine instead of the ZGELS routine and with rcond=-1. 
# Using instead rcond=10**(-40d0) in ZGELSD did not change the spectra much.
p "out_Beach.dat" u 1:(-1/3.14*$3) w l t 'python linalg.lstsq: gelsd?', "pade_A_64_gelsd" w l t '64 gelsd', "pade_A_64" w l t '64 gels', "pade_A_128" w l t '128 gels', "../exact.dat" u 1:(-1/pi*$3) w l lc 0 t 'exact'
