@echo off
REM to run >> run_test.cmd INPUTFILE.obj
@echo Running test script for ScatInterpol
REM first argument is the input obj file name 
REM set input_file=%1
REM set patch_size=%2
@setlocal EnableDelayedExpansion

cd build/Release


REM Experimenting with num_data_points 
REM for %%i in (100 1000 10000 100000 1000000 10000000)  do @(
REM 	call ScatInterpol.exe -output %%i_spheres_H_L_RE -scat 7 -slice 0.5 -res_x 1024 -res_y 1024 -proj 1 -data_fun 1 -num_data %%i -k 10 -r 0.1)


REM Experimenting with R^2 factor for Hardy method
REM for %%i in (0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 .8 0.9 0.99)  do @(
REM 	call ScatInterpol.exe -output %%i_spheres_H_G_MQ -scat 4 -slice 0.5 -res_x 1024 -res_y 1024 -proj 1 -data_fun 1 -num_data 1000 -k 10 -r %%i
REM 	)
 	
REM for %%i in (0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 .8 0.9 0.99)  do @(
REM  	call ScatInterpol.exe -output %%i_spheres_H_L_MQ -scat 5 -slice 0.5 -res_x 1024 -res_y 1024 -proj 1 -data_fun 1 -num_data 10000 -k 10 -r %%i
REM 	)
 
REM for %%i in (0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 .8 0.9 0.99)  do @(
REM  	call ScatInterpol.exe -output %%i_spheres_H_G_RE -scat 6 -slice 0.5 -res_x 1024 -res_y 1024 -proj 1 -data_fun 1 -num_data 1000 -k 10 -r %%i
REM 	)
	
REM for %%i in (0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 .8 0.9 0.99)  do @(
REM  	call ScatInterpol.exe -output %%i_spheres_H_L_RE -scat 7 -slice 0.5 -res_x 1024 -res_y 1024 -proj 1 -data_fun 1 -num_data 10000 -k 10 -r %%i
REM 	)


REM Experimenting with K nearest neighbours for localized methods 
for %%i in (3 6 9 12 15 18 21 24 27 30 33 36 39)  do @(
  	call ScatInterpol.exe -output %%i_spheres_S2_L -scat 3 -slice 0.5 -res_x 1024 -res_y 1024 -proj 1 -data_fun 1 -num_data 10000 -k %%i -r 0.1
 	)
	
for %%i in (3 6 9 12 15 18 21 24 27 30 33 36 39)  do @(
  	call ScatInterpol.exe -output %%i_spheres_H_L_MQ -scat 5 -slice 0.5 -res_x 1024 -res_y 1024 -proj 1 -data_fun 1 -num_data 10000 -k %%i -r 0.1
 	)
	
for %%i in (3 6 9 12 15 18 21 24 27 30 33 36 39)  do @(
  	call ScatInterpol.exe -output %%i_spheres_H_L_RE -scat 7 -slice 0.5 -res_x 1024 -res_y 1024 -proj 1 -data_fun 1 -num_data 10000 -k %%i -r 0.1
 	)
	
	