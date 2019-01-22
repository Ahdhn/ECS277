#include <string>
#include <iostream>
#include <stdint.h>

#include "common.h"
#include "grid.h"
#include "rc.h"
#include "image.h"

#define DIM 3

typedef uint32_t index_t;
typedef int32_t int_t;
typedef double data_t;

template<typename T, typename T_d>
void get_coordinates(T i, T j, T_d){

}

template<typename T, typename T_d>
void init_grid(std::string filename, Grid<T, T_d, DIM>*&grid){

	//read input file and fill in the data 

	//testing on five cocentric spheres 
	//sphere0, f=0.1, radius = 0.2
	//sphere1, f=0.2, radius = 0.5
	//sphere2, f=0.5, radius = 0.9

	T_d r0_sq(0.2*0.2), r1_sq(0.5*0.5), r2_sq(0.9*0.9);

	T n[3];
	n[0] = n[1] = n[2] = 2;
	T_d x[3];
	x[0] = x[1] = x[2] = 0;
	T_d s[3]; 
	s[0] = s[1] = s[2] = T_d(n[0]-1) / T_d(1);

	color_t background_color;
	background_color.r = 1;
	background_color.g = 1;
	background_color.b = 1;
	background_color.a = 1;

	static Grid<T, T_d, 3> my_grid(n, x, s, background_color);
		
	for (uint32_t i = 0; i < n[0]; i++){		
		for (uint32_t j = 0; j < n[1]; j++){
			for (uint32_t k = 0; k < n[2]; k++){

				T_d val = 0.0;

				T_d dist = Dist(0.5, 0.5, 0.5, i*s[i], j*s[j], k*s[k]);

				if (dist < r0_sq){ val = 0.1; }
				else if (dist < r1_sq){ val = 0.2; }
				else if (dist < r2_sq){ val = 0.5; }
				
				my_grid.fill_grid(i, j, k, val);
			}			
		}
	}

	grid = &my_grid;

	//return &my_grid;

}

int main(int argc, char**argv){
	
	
	
	std::string inputfilename = STRINGIFY(INPUT_DIR)"torus.obj";

	//Command line args 
	if (argc > 1){
		if (cmdOptionExists(argv, argc + argv, "-h")){
			std::cout << "\n\nUsage: RayCasting.exe <-option X>" << std::endl;
			std::cout << " -h:       Display this massage and exits\n" << std::endl;

			std::cout << " -input:      Input file. Input file should under the data/ subdirectory" << std::endl;
			std::cout << "              The default is " << inputfilename << std::endl;			
		}

		if (cmdOptionExists(argv, argc + argv, "-input")){
			inputfilename = STRINGIFY(INPUT_DIR) +
				std::string(getCmdOption(argv, argv + argc, "-input"));
			std::cout << "	input= " << inputfilename << std::endl;
		}
	}

	index_t image_res[2];
	int_t image_x0[3], image_xn[3];
	image_res[0] = image_res[1] = 512;
	image_x0[0] = image_x0[1] = image_x0[2] = -1;
	image_xn[0] = image_xn[1] = 2;
	image_xn[2] = -1;


	//Init grid 
	Grid<index_t, data_t, DIM> *my_grid = NULL;
	init_grid<index_t, data_t>(inputfilename, my_grid);

	

	Image<index_t, int_t, data_t> my_image(image_res, image_x0, image_xn);

	RC<index_t, int_t, data_t> rc(my_grid, &my_image, 100);
	
	rc.run_raycasting();
	

	//Do ray casting 
	//pointCSV("points.csv", num_points, points);	


	//Write output
	
	

	return 0;
}
