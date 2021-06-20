to compile this program in cygwins gfortran:

gfortran -o ../compiled/vaccine_agestratified_101.exe vaccine_age_stratified_101.F95 -fcheck=bounds -ffree-line-length-0

gfortran -o ../compiled/vaccine2019_sia_singlematrix.exe vaccine2019_sia_singlematrix.F95 -fcheck=bounds -ffree-line-length-0

gfortran -o name_of_compiled name_of_fortran options:

 -fcheck=bounds / fcheck=all
 -ffree-line-length-0 (doesn't cut the line-length after a certain number of characters)
 -Wall displays some warnings
 -Wtabs or Wno-tabs -depending on version of gfortran: doesn't display unnecessary tabs warning
