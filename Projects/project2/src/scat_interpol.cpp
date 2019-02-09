#include <string>
#include <iostream>
#include <stdint.h>

#include "common.h"
#include "grid.h"
#include "renderer.h"
#include "image.h"

#include "colormap.h"

#include "compare_images.h"

#include "scat_data.h"



#define DIM 3

typedef uint32_t index_t;
typedef int32_t int_t;
typedef double data_t;


template<typename T>
void read_file_raw(std::string filename, uint32_t x0, uint32_t x1, uint32_t x2,
	T*&data){

	//read a raw file where x0, x1, x2 are the dimension of the data in x,y,z
	//directions and data is a unallocated buffer to store the file content 

	//the type T depends on the size of the info in the raw file. 
	//for "head256.raw", it is 8-bit (unsigned char)
	
	if (data != NULL){
		PRINT_ERROR("read_file_raw() data should be uninitialized ");
		exit(EXIT_FAILURE);
	}
		
	uint32_t total_size = x0*x1*x2;
	data = (T*)malloc(total_size *sizeof(T));

	FILE *fp;

	if (!(fp = fopen(filename.c_str(), "rb"))){
		PRINT_ERROR("read_file_raw() Can not open " + filename);
		exit(EXIT_FAILURE);
	}

	if (fread(data, sizeof(char), total_size, fp) != total_size){
		PRINT_ERROR("read_file_raw() Failed to read " + filename);
		exit(EXIT_FAILURE);
	}

	fclose(fp);
}

template<typename T, typename T_d, typename bits_type>
void init_grid(std::string filename, Grid<T, T_d, DIM>*&grid, T_d&f_value_min,
	T_d&f_value_max, T n_grid[3], T_d grid_lower[3], T_d grid_upper[3],
	T_d bk_color[4], T_d l_color[4], ScatData<T, T_d>*data = nullptr, 
	INTERPOL_METHOD scat_data_interpol_method = INTERPOL_METHOD::S2_G, 
	const T K = 5, const T_d R = 0.1){
	
	//read input file and fill in the data 
	//unsigned short*raw_data = NULL;
	//unsigned char*raw_data = NULL;
	bits_type* raw_data = nullptr;

	if (filename.size() > 1){
		read_file_raw(filename, n_grid[0], n_grid[1], n_grid[2], raw_data);
	}

	//testing on five cocentric spheres 
	//sphere0, f=0.5, radius = 0.2 
	//sphere1, f=0.2, radius = 0.3
	//else f = 0.1
	//T_d r0_sq(0.2*0.2), r1_sq(0.4*0.4);
	//T_d f0(0.5), f1(0.4), f2(0.1);
	

	T n[3];	
	n[0] = n_grid[0];
	n[1] = n_grid[1];
	n[2] = n_grid[2];

	
	T_d x_lower[3], x_upper[3], s[3];

	x_lower[0] = grid_lower[0];
	x_lower[1] = grid_lower[1];
	x_lower[2] = grid_lower[2];

	x_upper[0] = grid_upper[0];
	x_upper[1] = grid_upper[1];
	x_upper[2] = grid_upper[2];

	s[0] = T_d(x_upper[0] - x_lower[0]) / T_d(n[0] - 1);
	s[1] = T_d(x_upper[1] - x_lower[1]) / T_d(n[1] - 1);
	s[2] = T_d(x_upper[2] - x_lower[2]) / T_d(n[2] - 1);
	
	color_t background_color;
	background_color.r = bk_color[0];
	background_color.g = bk_color[1];
	background_color.b = bk_color[2];
	background_color.a = bk_color[3];

	color_t light_color;
	light_color.r = l_color[0];
	light_color.g = l_color[1];
	light_color.b = l_color[2];
	light_color.a = l_color[3];

	T_d light_pos[3];
	light_pos[0] = -2.0;
	light_pos[1] = -2.0;
	light_pos[2] = 2.0;

	T_d light_dir[3];
	light_dir[0] = -1.0;
	light_dir[1] = -1.0;
	light_dir[2] = -1.0;

	color_t ambient_light;
	ambient_light.r = 25;
	ambient_light.g = 25;
	ambient_light.b = 25;
	ambient_light.a = 1.0;

	T_d ambient_const = 0.5;  
	T_d diffuse_const = 0.5;
	T_d specular_const = 1.0;
	
	T_d phong_power = 100.0;
	T_d attenuation_const = 2.0;
	T_d light_intensity = 1.0;
	

	static Grid<T, T_d, 3> my_grid(n, x_lower, s, background_color, light_color,
		light_pos, light_dir, phong_power, attenuation_const, light_intensity,
		ambient_light, ambient_const, diffuse_const, specular_const);

	f_value_min = std::numeric_limits<T_d>::max();
	f_value_max = -std::numeric_limits<T_d>::max();
		
	uint32_t total_size = n[0] * n[1] * n[2];


	uint32_t num_processed = 0;
	std::cout << " init_grid::filling the grid 0% ";
	std::vector<int> percentage(1, 0);

	
	for (uint32_t i = 0; i < total_size; i++){

		num_processed++;
		int percetnage_processed = int(100 * double(num_processed) /
			double(total_size));
		if (percetnage_processed % 10 == 0
			&& percetnage_processed != percentage.back()){
			std::cout << percetnage_processed << "%  ";
			percentage.push_back(percetnage_processed);
		}


		T_d my_data;
		if (raw_data != nullptr){
			my_data = T_d(raw_data[i]);
		}
		else if (data != nullptr){
			T_d xx, yy, zz;
			my_grid.get_location(i, xx, yy, zz);
			my_data = data->interpolate(xx, yy, zz,
				scat_data_interpol_method, K, R);
		}
		else{
			PRINT_ERROR("init_grid():: no vaild source to fill in the grid");
		}
		my_grid.fill_grid(i, my_data);
		f_value_min = std::min(f_value_min, my_data);
		f_value_max = std::max(f_value_max, my_data);
	}
	
	
	if (f_value_min != 0){
		//shift values so that f_value_min = 0		
		for (uint32_t i = 0; i < total_size; i++){
			my_grid.fill_grid(i, my_grid.get_value(i) - f_value_min);
		}
		f_value_max -= f_value_min;
		f_value_min = 0;
	}

	double avg(0), stddev(0);

	if (raw_data != nullptr){
		compute_avg_stddev(raw_data, total_size, avg, stddev);
		std::cout << "\n        Function max value=  " << f_value_max << std::endl;
		std::cout << "        Function min value=  " << f_value_min << std::endl;
		std::cout << "        Function avg value=  " << avg << std::endl;
		std::cout << "        Function std dev=    " << stddev << std::endl;
	}


	//my_grid.export_grid("grid.csv");
	grid = &my_grid;
}

template<typename T, typename T_d, typename bits_type>
void init_data(std::string filename, ScatData<T, T_d>*&data, T_d&f_value_min,
	T_d&f_value_max, const T num_data =100){
	
	f_value_min = std::numeric_limits<T_d>::max();
	f_value_max = -std::numeric_limits<T_d>::max();

	static ScatData<T, T_d> my_data(true, num_data);
	
	if (filename.size() <= 1){
		//analytical function
		
		for (uint32_t i = 0; i < num_data;i++){
			T_d x = T_d(rand()) / T_d(RAND_MAX);
			T_d y = T_d(rand()) / T_d(RAND_MAX);
			T_d z = T_d(rand()) / T_d(RAND_MAX);
			T_d f = sqrt(x*x + y*y + z*z);
			f_value_min = std::min(f_value_min, f);
			f_value_max= std::max(f_value_max, f);
			my_data.fill_data(i, x, y, z, f);
		}
	}
	else{

	}

	data = &my_data;
}


int main(int argc, char**argv){

	std::string inputfilename = STRINGIFY(INPUT_DIR)"stent/stent16.raw";

	std::string model_name = "lol";

	index_t samples_per_cell = 1;

	index_t n_grid[3]; n_grid[0] = n_grid[1] = n_grid[2] = 50;
		
	index_t bits = 16;

	data_t grid_lower[3]; grid_lower[0] = grid_lower[1] = grid_lower[2] = 0.0;

	data_t grid_upper[3]; grid_upper[0] = grid_upper[1] = grid_upper[2] = 1.0;	
		
	data_t bg_color[4]; bg_color[0] = bg_color[1] = bg_color[2] = 51.0; bg_color[3] = 1.0;

	data_t light_color[4]; light_color[0] = light_color[1] = light_color[2] = 255.0; light_color[3] = 1.0;

	index_t image_res[2]; image_res[0] = image_res[1] = 512;

	int projection = 1;

	if (model_name == "stent"){
		n_grid[0] = n_grid[1] = 512;
		n_grid[2] = 174;
		grid_upper[2] = 1.0;
		grid_upper[0] = grid_upper[1] = (512 * 0.8398) / (174 * 3.2);
		projection = 2;
	}

	if (model_name == "tooth"){
		n_grid[0] = n_grid[2] = 256;
		n_grid[1] = 161;

		grid_upper[0] = grid_upper[1] = 1.0;
		grid_upper[2] = 161.0 / 256.0;
		
		projection = -3;

	}

	if (model_name == "comunix"){
		n_grid[0] = n_grid[1] = 128;
		n_grid[2] = 83;

		grid_upper[0] = grid_upper[1] = (128.0 * 2.0593637) / (83 * 3.3750001);
		grid_upper[2] = 1.0;

		projection = -3;
	}

	if (model_name == "wbscan"){
		n_grid[0] = n_grid[1] = 128;
		n_grid[2] = 180;

		grid_upper[0] = grid_upper[1] = 128.0/180.0;
		grid_upper[2] = 1.0;

		projection = -3;
	}

	if (model_name == "vessels"){
		n_grid[0] = n_grid[1] = 384;
		n_grid[2] = 72;

		grid_upper[0] = grid_upper[1] = 1.0;
		grid_upper[2] = (72.0*1.5)/(384*1.0026);

		projection = -3;
	}

	if (model_name == "set1"){
		n_grid[0] = n_grid[1] = 256;
		n_grid[2] = 247;

		grid_upper[0] = grid_upper[1] = 1.0;
		grid_upper[2] = 247.0/256.0;

		projection = -3;
	}

	if (model_name == "liver16"){
		n_grid[0] = n_grid[1] = 256;
		n_grid[2] = 16;

		grid_upper[0] = grid_upper[1] = 1.0;
		grid_upper[2] = (256 * 1.6769) / (12 * 16);

		projection = -3;
	}

	if (model_name == "liver22"){
		n_grid[0] = n_grid[1] = 256;
		n_grid[2] = 22;

		grid_upper[0] = grid_upper[1] = 1.0;
		grid_upper[2] = (256 * 1.6769) / (10 * 22);

		projection = -3;
	}


	//Command line args 
	if (argc > 1){
		if (cmdOptionExists(argv, argc + argv, "-h")){
			std::cout << "\n\nUsage:  RayCasting.exe <-option X>" << std::endl;
			std::cout << " -h:        Display this massage and exits\n" << std::endl;

			std::cout << " -input:    Input file. Input file should under the data/ subdirectory" << std::endl;
			std::cout << "            The default is " << inputfilename << std::endl << std::endl;

			std::cout << " -model:    Model name used to name output images" << std::endl;
			std::cout << "            The default is " << model_name << std::endl << std::endl;

			std::cout << " -xn:       Number of grid points in the x-direction" << std::endl;
			std::cout << "            The default is " << n_grid[0] << std::endl << std::endl;

			std::cout << " -yn:       Number of grid points in the y-direction" << std::endl;
			std::cout << "            The default is " << n_grid[1] << std::endl << std::endl;

			std::cout << " -zn:       Number of grid points in the z-direction" << std::endl;
			std::cout << "            The default is " << n_grid[2] << std::endl << std::endl;

			std::cout << " -x_lower:  X coordinates of the gird lower corner" << std::endl;
			std::cout << "            The default is " << grid_lower[0] << std::endl << std::endl;

			std::cout << " -y_lower:  Y coordinates of the gird lower corner" << std::endl;
			std::cout << "            The default is " << grid_lower[1] << std::endl << std::endl;

			std::cout << " -z_lower:  Z coordinates of the gird lower corner" << std::endl;
			std::cout << "            The default is " << grid_lower[2] << std::endl << std::endl;

			std::cout << " -x_upper:  X coordinates of the gird upper corner" << std::endl;
			std::cout << "            The default is " << grid_upper[0] << std::endl << std::endl;

			std::cout << " -y_upper:  Y coordinates of the gird upper corner" << std::endl;
			std::cout << "            The default is " << grid_upper[1] << std::endl << std::endl;

			std::cout << " -z_upper:  Z coordinates of the gird upper corner" << std::endl;
			std::cout << "            The default is " << grid_upper[2] << std::endl << std::endl;

			std::cout << " -bits:     bit size of the input .raw file" << std::endl;
			std::cout << "            The default is " << bits << std::endl << std::endl;

			std::cout << " -bg_r:     Red channel for the background color bteween 0 to 255" << std::endl;
			std::cout << "            The default is " << bg_color[0] << std::endl << std::endl;

			std::cout << " -bg_g:     Green channel for the background color" << std::endl;
			std::cout << "            The default is " << bg_color[1] << std::endl << std::endl;
						
			std::cout << " -bg_b:     Blue channel for the background color bteween 0 to 255" << std::endl;
			std::cout << "            The default is " << bg_color[2] << std::endl << std::endl;

			std::cout << " -bg_a:     Alpha channel for the background color bteween 0 to 1" << std::endl;
			std::cout << "            The default is " << bg_color[3] << std::endl << std::endl;
			
			std::cout << " -light_r:  Red channel for the light color bteween 0 to 255" << std::endl;
			std::cout << "            The default is " << light_color[0] << std::endl << std::endl;

			std::cout << " -light_g:  Green channel for the light color bteween 0 to 255" << std::endl;
			std::cout << "            The default is " << light_color[1] << std::endl << std::endl;

			std::cout << " -light_b:  Blue channel for the light color bteween 0 to 255" << std::endl;
			std::cout << "            The default is " << light_color[2] << std::endl << std::endl;

			std::cout << " -light_a:  Alpha channel for the light color bteween 0 to 1" << std::endl;
			std::cout << "            The default is " << light_color[3] << std::endl << std::endl;

			std::cout << " -res_x:    Image resolution in the x-direction" << std::endl;
			std::cout << "            The default is " << image_res[0] << std::endl << std::endl;

			std::cout << " -res_y:    Image resolution in the y-direction" << std::endl;
			std::cout << "            The default is " << image_res[1] << std::endl << std::endl;

			std::cout << " -proj:     The projection direction where the value represents the axis of the normal ";
			std::cout << "            and the sign represents the direction e.g., -1 is projection along x-axis looking ";
			std::cout << "            at the grid back face. 3 is projection along z-axis pointing to the grid top face";
			std::cout << "            The default is " << projection << std::endl << std::endl;			

			std::cout << " -samples:  Number of samples per cell" << std::endl;
			std::cout << "            The default is " << samples_per_cell << std::endl << std::endl;



			exit(EXIT_SUCCESS);

		}
		if (cmdOptionExists(argv, argc + argv, "-input")){
			inputfilename = STRINGIFY(INPUT_DIR) +
				std::string(getCmdOption(argv, argv + argc, "-input"));
		}

		if (cmdOptionExists(argv, argc + argv, "-model")){
			model_name = std::string(getCmdOption(argv, argv + argc, "-model"));
		}

		if (cmdOptionExists(argv, argc + argv, "-samples")){
			samples_per_cell = atoi(getCmdOption(argv, argv + argc, "-samples"));
		}
		if (cmdOptionExists(argv, argc + argv, "-xn")){
			n_grid[0] = atoi(getCmdOption(argv, argv + argc, "-xn"));
		}
		if (cmdOptionExists(argv, argc + argv, "-yn")){
			n_grid[1] = atoi(getCmdOption(argv, argv + argc, "-yn"));
		}
		if (cmdOptionExists(argv, argc + argv, "-zn")){
			n_grid[2] = atoi(getCmdOption(argv, argv + argc, "-zn"));
		}

		if (cmdOptionExists(argv, argc + argv, "-x_lower")){
			grid_lower[0] = atoi(getCmdOption(argv, argv + argc, "-x_lower"));
		}

		if (cmdOptionExists(argv, argc + argv, "-y_lower")){
			grid_lower[1] = atoi(getCmdOption(argv, argv + argc, "-y_lower"));
		}

		if (cmdOptionExists(argv, argc + argv, "-z_lower")){
			grid_lower[2] = atoi(getCmdOption(argv, argv + argc, "-z_lower"));
		}

		if (cmdOptionExists(argv, argc + argv, "-x_upper")){
			grid_upper[0] = atoi(getCmdOption(argv, argv + argc, "-x_upper"));
		}

		if (cmdOptionExists(argv, argc + argv, "-y_upper")){
			grid_upper[1] = atoi(getCmdOption(argv, argv + argc, "-y_upper"));
		}

		if (cmdOptionExists(argv, argc + argv, "-z_upper")){
			grid_upper[2] = atoi(getCmdOption(argv, argv + argc, "-z_upper"));
		}

		if (cmdOptionExists(argv, argc + argv, "-bk_r")){
			bg_color[0] = atoi(getCmdOption(argv, argv + argc, "-bk_r"));
		}

		if (cmdOptionExists(argv, argc + argv, "-bk_g")){
			bg_color[1] = atoi(getCmdOption(argv, argv + argc, "-bk_g"));
		}

		if (cmdOptionExists(argv, argc + argv, "-bk_b")){
			bg_color[2] = atoi(getCmdOption(argv, argv + argc, "-bk_b"));
		}

		if (cmdOptionExists(argv, argc + argv, "-bk_a")){
			bg_color[3] = atoi(getCmdOption(argv, argv + argc, "-bk_a"));
		}

		if (cmdOptionExists(argv, argc + argv, "-light_r")){
			light_color[0] = atoi(getCmdOption(argv, argv + argc, "-light_r"));
		}

		if (cmdOptionExists(argv, argc + argv, "-light_g")){
			light_color[1] = atoi(getCmdOption(argv, argv + argc, "-light_g"));
		}

		if (cmdOptionExists(argv, argc + argv, "-light_b")){
			light_color[2] = atoi(getCmdOption(argv, argv + argc, "-light_b"));
		}

		if (cmdOptionExists(argv, argc + argv, "-light_a")){
			light_color[3] = atoi(getCmdOption(argv, argv + argc, "-light_a"));
		}

		if (cmdOptionExists(argv, argc + argv, "-res_x")){
			image_res[0] = atoi(getCmdOption(argv, argv + argc, "-res_x"));
		}

		if (cmdOptionExists(argv, argc + argv, "-res_y")){
			image_res[1] = atoi(getCmdOption(argv, argv + argc, "-res_y"));
		}

		if (cmdOptionExists(argv, argc + argv, "-bits")){
			bits = atoi(getCmdOption(argv, argv + argc, "-bits"));
		}

		if (cmdOptionExists(argv, argc + argv, "-proj")){
			projection = atoi(getCmdOption(argv, argv + argc, "-proj"));
		}
	}

	std::cout << "	input= " << inputfilename << std::endl;
	std::cout << "	model= " << model_name << std::endl;
	std::cout << "	xn= " << n_grid[0] << std::endl;
	std::cout << "	yn= " << n_grid[1] << std::endl;
	std::cout << "	zn= " << n_grid[2] << std::endl;
	std::cout << "	lower corner= (" << grid_lower[0] << ", " << grid_lower[1] << ", "
		<< grid_lower[2] << ") " << std::endl;
	std::cout << "	upper corner= (" << grid_upper[0] << ", " << grid_upper[1] << ", "
		<< grid_upper[2] << ") " << std::endl;
	std::cout << "	bits= " << bits << std::endl;
	std::cout << "	background color= (" << bg_color[0] << ", " << bg_color[1] << ", "
		<< bg_color[2] << ", " << bg_color[3] << ") " << std::endl;
	std::cout << "	light color= (" << light_color[0] << ", " << light_color[1] << ", "
		<< light_color[2] << ", " << light_color[3] << ") " << std::endl;
	std::cout << "	projection= " << projection << std::endl;
	std::cout << "	samples= " << samples_per_cell << std::endl;
	std::cout << "	image_resolution= " << image_res[0] << " X " << image_res[1] << std::endl;
		

	//setup the image using the right projection
	data_t image_x0[3], image_xn[3], image_normal[3];		
	data_t frac = 0.1 *(grid_upper[0] - grid_lower[0]);
	std::string image_prefix;
	if (abs(projection) == 3){
		//Image point to +ve z-axis 
		data_t dist_z = grid_upper[2] - grid_lower[2];

		image_x0[0] = grid_lower[0] - frac;
		image_x0[1] = grid_lower[1] - frac;
		image_x0[2] = (projection > 0) ? grid_lower[2] - dist_z : grid_upper[2] + dist_z;

		image_xn[0] = grid_upper[0] + frac;
		image_xn[1] = grid_upper[1] + frac;
		image_xn[2] = (projection > 0) ? grid_lower[2] - dist_z : grid_upper[2] + dist_z;

		image_normal[0] = 0;
		image_normal[1] = 0;
		image_normal[2] = (projection < 0) ? -1 : 1;

		image_prefix = model_name + ((projection < 0) ? "-z" : "+z");
	}
	else if (abs(projection) == 2){
		//Image point to -ve y-axis 
		data_t dist_y = grid_upper[1] - grid_lower[1];

		image_x0[0] = grid_lower[0] - frac;
		image_x0[1] = (projection > 0) ? grid_lower[1] - dist_y : grid_upper[1] + dist_y;
		image_x0[2] = grid_lower[2] - frac;

		image_xn[0] = grid_upper[0] + frac;
		image_xn[1] = (projection > 0) ? grid_lower[1] - dist_y : grid_upper[1] + dist_y;
		image_xn[2] = grid_upper[2] + frac;

		image_normal[0] = 0;
		image_normal[1] = (projection < 0) ? -1 : 1;
		image_normal[2] = 0;

		image_prefix = model_name + ((projection < 0) ? "-y" : "+y");
	}
	else if (abs(projection) == 1){
		//Image point to -ve x-axis 
		data_t dist_x = grid_upper[0] - grid_lower[0];

		image_x0[0] = (projection > 0) ? grid_lower[0] - dist_x : grid_upper[0] + dist_x;
		image_x0[1] = grid_lower[1] - frac;
		image_x0[2] = grid_lower[2] - frac;

		image_xn[0] = (projection > 0) ? grid_lower[0] - dist_x : grid_upper[0] + dist_x;
		image_xn[1] = grid_upper[1] + frac;
		image_xn[2] = grid_upper[2] + frac;

		image_normal[0] = (projection < 0) ? -1 : 1;
		image_normal[1] = 0;
		image_normal[2] = 0;

		image_prefix = model_name + ((projection < 0) ? "-x" : "+x");
	}
	else {
		PRINT_ERROR("Invalid projection direction");
	}

	bool flip_vertical(false), flip_horizontal(false);


	//Init data 
	data_t data_f_value_min(0), data_f_value_max(0);
	ScatData<index_t, data_t> *my_data= NULL;
	init_data<index_t, data_t, unsigned char>("", my_data, data_f_value_min, 
		data_f_value_max, 1000);

	my_data->export_data("data.csv");
		
	my_data->precompute(INTERPOL_METHOD::H_G_RE);
		
	//to indicate that the grid data comes from scatter data
	bits = 0;

	//Init grid 
	Grid<index_t, data_t, DIM> *my_grid = NULL;
	data_t f_value_min(0), f_value_max(0);

	if (bits == 0){
		//for testing with analytical function 
		init_grid<index_t, data_t, unsigned char>("", my_grid,
			f_value_min, f_value_max, n_grid, grid_lower, grid_upper, bg_color,
			light_color, my_data, INTERPOL_METHOD::H_G_RE, 10, 0.1);
		
	}
	else if (bits == 8){
		init_grid<index_t, data_t, unsigned char>(inputfilename, my_grid, 
			f_value_min, f_value_max, n_grid, grid_lower, grid_upper, bg_color,
			light_color);
	}
	else if (bits == 16) {
		init_grid<index_t, data_t, unsigned short>(inputfilename, my_grid, 
			f_value_min, f_value_max, n_grid, grid_lower, grid_upper, bg_color,
			light_color);
	}
	else{
		PRINT_ERROR("Invalid bits size" + std::to_string(bits));
	}
	
	//image for the trilinear case
	TGAImage tga_image_linear(image_res[0], image_res[1], TGAImage::RGBA);
	Image<index_t, int_t, data_t> my_image_linear(image_res, image_x0, image_xn, 
		image_normal, &tga_image_linear, image_prefix + "_.tga",
		flip_vertical, flip_horizontal);	

	//init the ray-caster   
	Renderer<index_t, int_t, data_t> my_renderer(my_grid, samples_per_cell);
	
	//skip threshold 
	data_t min_threshold(f_value_min), max_threshold(f_value_max);

	if (model_name == "stent"){
		min_threshold = 1000;
		max_threshold = 2000;
	}

	if (model_name == "tooth"){
		min_threshold = 500;
		max_threshold = 950;
	}

	if (model_name == "comunix"){
		min_threshold = 930;
		max_threshold = 10000;
	}

	if (model_name == "vessels"){
		min_threshold = 1000;
		max_threshold = 3000;
	}


	if (model_name == "set1"){
		min_threshold = 10;
		max_threshold = 90;
	}

	if (model_name == "liver22"){
		min_threshold = 100;
		max_threshold = 450;
	}

	//Alpha Transfer Function 
	auto alpha_transfer_func = [f_value_min, f_value_max, min_threshold, max_threshold,
	model_name](data_t f_value){

		/*if (f_value - EPSILON < f_value_min){
		return 0.0;
		}
		else if (f_value + EPSILON > f_value_max){
		return 1.0;
		}

		data_t val = (f_value - f_value_min) / (f_value_max - f_value_min);*/

		if (f_value - EPSILON < min_threshold){
			return 0.0;
		}
		else if (f_value + EPSILON > max_threshold){
			return 1.0;
		}

		data_t val = (f_value - min_threshold) / (max_threshold - min_threshold);

		//with stent 
		if (model_name == "stent"){
			if (f_value > min_threshold && f_value < max_threshold){
				val *= 1.4;
				val = std::min(val, 1.0);
			}
		}		

		return val;
	};
	
	//Color Transfer Function 
	auto color_transfer_func = [f_value_min, f_value_max, min_threshold, max_threshold,
		model_name](data_t f_value){
		color_t color;
		COLOR_MAP cm = COLOR_MAP::RAINBOW;
		

		if (model_name == "tooth"){
			cm = COLOR_MAP::RAINBOW;
		}

		colormap(cm,
			(f_value - min_threshold) / (max_threshold - min_threshold) //f_value / f_value_max
			,color.r, color.g, color.b);

		return color;
	};

	//Skip Transfer Function 
	auto skip_transfer_func = [f_value_min, f_value_max, min_threshold, max_threshold](data_t f_value){		
		if (f_value < min_threshold || f_value > max_threshold){
			return true;
		}
		return false;
	};

	//do slicing 
	my_renderer.slice(&my_image_linear, color_transfer_func,0.5);

	//do the ray casting 
	//my_renderer.run_raycasting(&my_image_linear, color_transfer_func, 
	//	alpha_transfer_func, skip_transfer_func, INTERPOL_TYPE::TRILINEAR);	

	
		
	
	system("pause");

	return 0;
}
