#include <string>
#include <iostream>
#include <stdint.h>
#include <chrono>

#include "common.h"
#include "image.h"
#include "colormap.h"
#include "compare_images.h"

#include "scat_data_2d.h"
#include "vorodel.h"


#define DIM 3

typedef int32_t index_t;
typedef double data_t;

void points_gen(ScatData2D<index_t, data_t>*my_scat_data, const index_t 
	num_points){

	//generate four corner points and then fill the interior randomly 
	data_t x, y, f;

	//(0,0)
	x = 0;
	y = 0;
	f = sqrt(x*x + y*y);
	my_scat_data->add_data(x, y, f);

	//(1,0)
	x = 1;
	y = 0;
	f = sqrt(x*x + y*y);
	my_scat_data->add_data(x, y, f);

	//(1,1)
	x = 1;
	y = 1;
	f = sqrt(x*x + y*y);
	my_scat_data->add_data(x, y, f);

	//(0,1)
	x = 0;
	y = 1;
	f = sqrt(x*x + y*y);
	my_scat_data->add_data(x, y, f);

	//(0.5,0)
	x = 0.5;
	y = 0;
	f = sqrt(x*x + y*y);
	my_scat_data->add_data(x, y, f);
	
	//(0.0,0.5)
	x = 0;
	y = 0.5;
	f = sqrt(x*x + y*y);
	my_scat_data->add_data(x, y, f);

	//(0.5,1.0)
	x = 0.5;
	y = 1.0;
	f = sqrt(x*x + y*y);
	my_scat_data->add_data(x, y, f);

	//(1.0,0.5)
	x = 1.0;
	y = 0.5;
	f = sqrt(x*x + y*y);
	my_scat_data->add_data(x, y, f);



	for (index_t i = 0; i < num_points; i++){
		data_t x = data_t(rand()) / data_t(RAND_MAX);
		data_t y = data_t(rand()) / data_t(RAND_MAX);

		x *= 0.9;
		y *= 0.9;
		x += 0.05;
		y += 0.05;

		if (i == 0){
			x = 0.5;
			y = 0.5;
		}

		data_t f = sqrt(x*x + y*y);
		my_scat_data->add_data(x, y, f);
	}
}

int main(int argc, char**argv){

	
	index_t num_points = 500;	

	index_t knn = 60;


	//Command line args 
	if (argc > 1){
		if (cmdOptionExists(argv, argc + argv, "-h")){
			std::cout << "\n\nUsage:  RayCasting.exe <-option X>" << std::endl;
			std::cout << " -h:        Display this massage and exits\n" << std::endl;
			std::cout << " -points:   Input number of points. Default is " << num_points << "\n" << std::endl;
			std::cout << " -K:        K nearest neighbour used for Delaunay construction. Default is " << knn << "\n" << std::endl;
			exit(EXIT_SUCCESS);
		}	
		
		if (cmdOptionExists(argv, argc + argv, "-points")){
			num_points = atoi(getCmdOption(argv, argv + argc, "-input"));
			std::cout << "	num_points= " << num_points << std::endl;
		}	
	}
		
	
	//Init data 	
	ScatData2D<index_t, data_t> my_scat_data(knn);	
	points_gen(&my_scat_data, num_points);	
	my_scat_data.build_kdtree();

	my_scat_data.run_multi_res();



	//Voronoi shading lambda function 
	/*auto voro_shading = [&my_scat_data](index_t p){
		data_t f = my_scat_data.get_data_value(p);
		color_t color;
		COLOR_MAP cm = COLOR_MAP::RAINBOW;
		colormap(cm,
			(f - my_scat_data.get_min_val()) / 
			(my_scat_data.get_max_val() - my_scat_data.get_min_val()) 
			, color.r, color.g, color.b);
		color.a = 1.0;
		return color;

	};
	//construct delaunay over the data 
	VoroDel<index_t, data_t> my_vorodel(knn, my_scat_data.get_num_data());

	std::vector<index_t>points(my_scat_data.get_num_data());
	fill_sequential_numbers(points.data(), points.size());
	rand_permute_subset(points.data(), points.size());

	my_vorodel.construct_del(&my_scat_data,points);
	my_vorodel.plot("test.ps", &my_scat_data, false, true, voro_shading);*/

	//my_vorodel.sibson_interpol(&my_scat_data,1);


	

	return 0;
}
