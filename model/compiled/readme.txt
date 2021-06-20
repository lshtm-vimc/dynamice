Call the executable as follows:

./executable years tstep folder country psa OS debug_debug debug_compartments debug_age debug_timepoints debug_relative

 ./executable		-	name of executable (vaccine.exe)
 years				-	years that the model is ran (should correspond to input-data)
 tstep				-	tstep that is used in the model (should correspond to beta-matrices used)
 folder				-	foldername that is parent directory of input and output folders
 country			-	ISO3 name of country that is used in model (should correspond to input-data)
 PSA				-	ran of psa (0 if psa is disabled; should correspond to folderstructure of input-data)
 OS					-	operating system ("WIN" for windows, others don't matter)
 
 debug_debug		-	should debug-data be generated
						0	-	no
						1	-	only spin-up
						2	-	only model after spin-up
						3	-	spin-up + model after spin_up
					
 debug_compartments	-	type of data
						0	-	number of cases
						1	-	number of people in compartments
						
 debug_age			-	how should age-groups be reported
						0	-	all annual age-bands
						1	-	0-2yo in weekly age-bands, 3-100yo in annual age bands
						2	-	sum all age-bands
	
 debug_timepoints	-	how should time be reported
						0	-	per year
						1	-	per timestep, but only report first 25% of total timesteps
						2	-	per timestep, report all timesteps
						3 	- 	per week, report every 19-th timestep (19.18 timesteps in a week)XX
						
 debug_relative		-	how should cases be reported
						0	-	as absolute numbers
						1	-	as proportion
						
example
./vaccine.exe 121 1000 20171031_v0_s0_deter COD 0 WIN 2 0 2 1 0
./vaccine 121 1000 20171031_v0_s0_deter COD 0 LIN 0 0 0 0 0
