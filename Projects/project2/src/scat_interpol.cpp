#include <string>
#include <iostream>
#include <stdint.h>
#include <chrono>

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


template<class T>
inline void unflat_id(T flat_id, T&i, T&j, T&k, const T n[3]){
	k = flat_id / (n[0] * n[1]);
	flat_id -= (k * n[0] * n[1]);
	j = flat_id / n[0];
	i = flat_id % n[0];
}
template<class T, class T_d>
void get_location(const T i, const T j, const T k, T_d&x, T_d&y, T_d&z, 
	const T_d s[3], const T_d x0[]){
	x = i*s[0] + x0[0];
	y = j*s[1] + x0[1];
	z = k*s[2] + x0[2];
}


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
	T_d bk_color[4], T_d l_color[4], ScatData<T, T_d>*data, 
	INTERPOL_METHOD scat_data_interpol_method = INTERPOL_METHOD::S2_G, 
	const T K = 5, const T_d R = 0.1){
		

	bits_type* raw_data = NULL;
	if (filename.size() > 1){
		read_file_raw(filename, n_grid[0], n_grid[1], n_grid[2], raw_data);
	}
	
	T_d s[3];

	s[0] = T_d(grid_upper[0] - grid_lower[0]) / T_d(n_grid[0] - 1);
	s[1] = T_d(grid_upper[1] - grid_lower[1]) / T_d(n_grid[1] - 1);
	s[2] = T_d(grid_upper[2] - grid_lower[2]) / T_d(n_grid[2] - 1);
	
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
	

	static Grid<T, T_d, 3> my_grid(n_grid, grid_lower, s, background_color, light_color,
		light_pos, light_dir, phong_power, attenuation_const, light_intensity,
		ambient_light, ambient_const, diffuse_const, specular_const);

	f_value_min = std::numeric_limits<T_d>::max();
	f_value_max = -std::numeric_limits<T_d>::max();
		
	uint32_t total_size = n_grid[0] * n_grid[1] * n_grid[2];
	


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

		if (raw_data != NULL){
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

	if (!raw_data){
		free(raw_data);
	}
	grid = &my_grid;
}

template<typename T, typename T_d, typename bits_type, typename AnalyticalFunc>
void init_data(std::string filename, ScatData<T, T_d>*&data, T n_grid[3], 
	T_d grid_lower[3], T_d grid_upper[3],
	T_d&f_value_min, T_d&f_value_max, const T num_data,
	AnalyticalFunc analytical_function){
	
	f_value_min = std::numeric_limits<T_d>::max();
	f_value_max = -std::numeric_limits<T_d>::max();

	static ScatData<T, T_d> my_data(true, num_data);
	
	if (filename.size() <= 1){
		//analytical function defined in unit box 
		for (uint32_t i = 0; i < num_data;i++){
			T_d x = T_d(rand()) / T_d(RAND_MAX);
			T_d y = T_d(rand()) / T_d(RAND_MAX);
			T_d z = T_d(rand()) / T_d(RAND_MAX);
			
			//T_d f = sqrt(x*x + y*y + z*z);
			T_d f = analytical_function(x, y, z);

			f_value_min = std::min(f_value_min, f);
			f_value_max= std::max(f_value_max, f);
			my_data.fill_data(i, x, y, z, f);
		}
	}
	else{
		//function from input raw file 
		bits_type* raw_data = nullptr;
		read_file_raw(filename, n_grid[0], n_grid[1], n_grid[2], raw_data);

		T total_size = n_grid[0] * n_grid[1] * n_grid[2];

		T_d s[3];
		s[0] = T_d(grid_upper[0] - grid_lower[0]) / T_d(n_grid[0] - 1);
		s[1] = T_d(grid_upper[1] - grid_lower[1]) / T_d(n_grid[1] - 1);
		s[2] = T_d(grid_upper[2] - grid_lower[2]) / T_d(n_grid[2] - 1);

		std::vector<uint32_t> rand_num(total_size);
		fill_sequential_numbers(rand_num.data(), rand_num.size());
		rand_permute_subset(rand_num.data(), rand_num.size());
		for (uint32_t q = 0; q < num_data; q++){
			T id = rand_num[q];

			//map the flat id into 3d 
			T i, j, k;
			unflat_id(id, i, j, k, n_grid);
			T_d x, y, z;
			get_location(i, j, k, x, y, z, s, grid_lower);

			T_d f = T_d(raw_data[id]);

			f_value_min = std::min(f_value_min, f);
			f_value_max = std::max(f_value_max, f);
			my_data.fill_data(q, x, y, z, f);
		}

		if (f_value_min != 0){
			//shift values so that f_value_min = 0		
			for (uint32_t i = 0; i < total_size; i++){
				my_data.fill_data(i, my_data.get_data_value(i) - f_value_min);
			}
			f_value_max -= f_value_min;
			f_value_min = 0;
		}

		if (raw_data != nullptr){
			double avg(0), stddev(0);
			compute_avg_stddev(raw_data, total_size, avg, stddev);
			std::cout << "\n      Function max value=  " << f_value_max << std::endl;
			std::cout << "        Function min value=  " << f_value_min << std::endl;
			std::cout << "        Function avg value=  " << avg << std::endl;
			std::cout << "        Function std dev=    " << stddev << std::endl;
		}

		free(raw_data);
	}
	data = &my_data;
}


int main(int argc, char**argv){

	std::string inputfilename = STRINGIFY(INPUT_DIR)"stent/stent16.raw";

	std::string model_name = "stent";

	std::string output_image_name = "";

	index_t samples_per_cell = 1;

	index_t n_grid[3]; n_grid[0] = n_grid[1] = n_grid[2] = 512;
		
	index_t bits = 16;

	data_t grid_lower[3]; grid_lower[0] = grid_lower[1] = grid_lower[2] = 0.0;

	data_t grid_upper[3]; grid_upper[0] = grid_upper[1] = grid_upper[2] = 1.0;	
		
	data_t bg_color[4]; bg_color[0] = bg_color[1] = bg_color[2] = 51.0; bg_color[3] = 1.0;

	data_t light_color[4]; light_color[0] = light_color[1] = light_color[2] = 255.0; light_color[3] = 1.0;

	index_t image_res[2]; image_res[0] = image_res[1] = 1024;

	data_t slice_depth = 0.5;

	int projection = 1;

	int scat_data_interpol_method = 3;

	index_t data_fun = 0;

	INTERPOL_METHOD scat_data_interpol;

	uint32_t num_data = 1000000;
	
	uint32_t K = 10;

	double R = 0.1;

	if (model_name == "tooth"){
		n_grid[0] = n_grid[2] = 256;
		n_grid[1] = 161;
		grid_upper[0] = grid_upper[1] = 1.0;
		grid_upper[2] = 161.0 / 256.0;
		projection = -3;
		slice_depth = 0.25;
	}

	if (model_name == "stent"){
		n_grid[0] = n_grid[1] = 512;
		n_grid[2] = 174;
		grid_upper[2] = 1.0;
		grid_upper[0] = grid_upper[1] = (512 * 0.8398) / (174 * 3.2);
		projection = 2;	
		
		slice_depth = 0.3;
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

			std::cout << " -output:   Output image name. If left empty, the name will be combination " << std::endl;
			std::cout << "            of the model name and view direction which is the default" << std::endl << std::endl;

			std::cout << " -slice:    The slice depth perpedicular to the view direction normalized to the grid depth (i.e.,[0,1])" << std::endl;
			std::cout << "            The default is " << slice_depth << std::endl << std::endl;

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

			//std::cout << " -bg_r:     Red channel for the background color bteween 0 to 255" << std::endl;
			//std::cout << "            The default is " << bg_color[0] << std::endl << std::endl;
			//
			//std::cout << " -bg_g:     Green channel for the background color" << std::endl;
			//std::cout << "            The default is " << bg_color[1] << std::endl << std::endl;
			//			
			//std::cout << " -bg_b:     Blue channel for the background color bteween 0 to 255" << std::endl;
			//std::cout << "            The default is " << bg_color[2] << std::endl << std::endl;
			//
			//std::cout << " -bg_a:     Alpha channel for the background color bteween 0 to 1" << std::endl;
			//std::cout << "            The default is " << bg_color[3] << std::endl << std::endl;
			
			//std::cout << " -light_r:  Red channel for the light color bteween 0 to 255" << std::endl;
			//std::cout << "            The default is " << light_color[0] << std::endl << std::endl;
			//
			//std::cout << " -light_g:  Green channel for the light color bteween 0 to 255" << std::endl;
			//std::cout << "            The default is " << light_color[1] << std::endl << std::endl;
			//
			//std::cout << " -light_b:  Blue channel for the light color bteween 0 to 255" << std::endl;
			//std::cout << "            The default is " << light_color[2] << std::endl << std::endl;
			//
			//std::cout << " -light_a:  Alpha channel for the light color bteween 0 to 1" << std::endl;
			//std::cout << "            The default is " << light_color[3] << std::endl << std::endl;

			std::cout << " -res_x:    Image resolution in the x-direction" << std::endl;
			std::cout << "            The default is " << image_res[0] << std::endl << std::endl;

			std::cout << " -res_y:    Image resolution in the y-direction" << std::endl;
			std::cout << "            The default is " << image_res[1] << std::endl << std::endl;

			std::cout << " -proj:     The projection direction where the value represents the axis of the normal " << std::endl;
			std::cout << "            and the sign represents the direction e.g., -1 is projection along x-axis looking " << std::endl;
			std::cout << "            at the grid back face. 3 is projection along z-axis pointing to the grid top face" << std::endl;
			std::cout << "            The default is " << projection << std::endl << std::endl;			

			//std::cout << " -samples:  Number of samples per cell" << std::endl;
			//std::cout << "            The default is " << samples_per_cell << std::endl << std::endl;

			std::cout << " -scat:     Indicates scatter data interpolation method where" << std::endl;
			std::cout << "            0 = Global Shepard 1" << std::endl;
			std::cout << "            1 = Localized Shepard 1" << std::endl;
			std::cout << "            2 = Global Shepard 2" << std::endl;
			std::cout << "            3 = Localized Shepard 2" << std::endl;
			std::cout << "            4 = Global Multiquadric Hardy's" << std::endl;
			std::cout << "            5 = Localized Multiquadric Hardy's" << std::endl;
			std::cout << "            6 = Global Reciprocal Hardy's" << std::endl;
			std::cout << "            7 = Localized Reciprocal Hardy's" << std::endl;
			std::cout << "            The default is " << scat_data_interpol_method << std::endl;

			std::cout << " -data_fun: Indicates the type of the scatter data where " << std::endl;
			std::cout << "            0 = Data from input file (.raw file)" << std::endl;
			std::cout << "            1 = sqrt(x^2 + y^2 + z^2)" << std::endl;

			std::cout << " -num_data: Number of scatter data points." << std::endl;
			std::cout << "            The default is " << num_data << std::endl << std::endl;

			std::cout << " -k:        Number of K nearest neighbours used for localized methods" << std::endl;
			std::cout << "            The default is " << K << std::endl << std::endl;

			std::cout << " -r:        Constant used with Hardy method" << std::endl;
			std::cout << "            The default is " << R << std::endl << std::endl;
			
			exit(EXIT_SUCCESS);
		}

		if (cmdOptionExists(argv, argc + argv, "-input")){
			inputfilename = STRINGIFY(INPUT_DIR) +
				std::string(getCmdOption(argv, argv + argc, "-input"));
			std::cout << "	input= " << inputfilename << std::endl;
		}
		if (cmdOptionExists(argv, argc + argv, "-model")){
			model_name = std::string(getCmdOption(argv, argv + argc, "-model"));
			std::cout << "	model= " << model_name << std::endl;
		}		
		if (cmdOptionExists(argv, argc + argv, "-output")){
			output_image_name = std::string(getCmdOption(argv, argv + argc, "-output"));
		}
		if (cmdOptionExists(argv, argc + argv, "-samples")){
			samples_per_cell = atoi(getCmdOption(argv, argv + argc, "-samples"));
			std::cout << "	samples= " << samples_per_cell << std::endl;
		}

		if (cmdOptionExists(argv, argc + argv, "-data_fun")){
			data_fun = atoi(getCmdOption(argv, argv + argc, "-data_fun"));
		}

		if (cmdOptionExists(argv, argc + argv, "-slice")){
			slice_depth = atof(getCmdOption(argv, argv + argc, "-slice"));
		}

		if (cmdOptionExists(argv, argc + argv, "-xn")){
			n_grid[0] = atoi(getCmdOption(argv, argv + argc, "-xn"));
			std::cout << "	xn= " << n_grid[0] << std::endl;
		}
		if (cmdOptionExists(argv, argc + argv, "-yn")){
			n_grid[1] = atoi(getCmdOption(argv, argv + argc, "-yn"));
			std::cout << "	yn= " << n_grid[1] << std::endl;
		}
		if (cmdOptionExists(argv, argc + argv, "-zn")){
			n_grid[2] = atoi(getCmdOption(argv, argv + argc, "-zn"));
			std::cout << "	zn= " << n_grid[2] << std::endl;
		}

		if (cmdOptionExists(argv, argc + argv, "-x_lower")){
			grid_lower[0] = atoi(getCmdOption(argv, argv + argc, "-x_lower"));
		}

		if (cmdOptionExists(argv, argc + argv, "-y_lower")){
			grid_lower[1] = atoi(getCmdOption(argv, argv + argc, "-y_lower"));
		}

		if (cmdOptionExists(argv, argc + argv, "-z_lower")){
			grid_lower[2] = atoi(getCmdOption(argv, argv + argc, "-z_lower"));
			std::cout << "	lower corner= (" << grid_lower[0] << ", " << grid_lower[1] << ", "
				<< grid_lower[2] << ") " << std::endl;
		}

		if (cmdOptionExists(argv, argc + argv, "-x_upper")){
			grid_upper[0] = atoi(getCmdOption(argv, argv + argc, "-x_upper"));
		}

		if (cmdOptionExists(argv, argc + argv, "-y_upper")){
			grid_upper[1] = atoi(getCmdOption(argv, argv + argc, "-y_upper"));
		}

		if (cmdOptionExists(argv, argc + argv, "-z_upper")){
			grid_upper[2] = atoi(getCmdOption(argv, argv + argc, "-z_upper"));
			std::cout << "	upper corner= (" << grid_upper[0] << ", " << grid_upper[1] << ", "
				<< grid_upper[2] << ") " << std::endl;
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
			std::cout << "	background color= (" << bg_color[0] << ", " 
				<< bg_color[1] << ", "
				<< bg_color[2] << ", " << bg_color[3] << ") " << std::endl;
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
			std::cout << "	light color= (" << light_color[0] << ", "
				<< light_color[1] << ", "
				<< light_color[2] << ", " << light_color[3] << ") " << std::endl;
		}

		if (cmdOptionExists(argv, argc + argv, "-res_x")){
			image_res[0] = atoi(getCmdOption(argv, argv + argc, "-res_x"));
		}

		if (cmdOptionExists(argv, argc + argv, "-res_y")){
			image_res[1] = atoi(getCmdOption(argv, argv + argc, "-res_y"));
		}

		if (cmdOptionExists(argv, argc + argv, "-bits")){
			bits = atoi(getCmdOption(argv, argv + argc, "-bits"));
			std::cout << "	bits= " << bits << std::endl;
		}

		if (cmdOptionExists(argv, argc + argv, "-proj")){
			projection = atoi(getCmdOption(argv, argv + argc, "-proj"));
		}

		if (cmdOptionExists(argv, argc + argv, "-scat")){
			scat_data_interpol_method = 
				atoi(getCmdOption(argv, argv + argc, "-scat"));			
		}
		if (cmdOptionExists(argv, argc + argv, "-num_data")){
			num_data = atoi(getCmdOption(argv, argv + argc, "-num_data"));
		}
		if (cmdOptionExists(argv, argc + argv, "-k")){
			K = atoi(getCmdOption(argv, argv + argc, "-k"));
		}
		if (cmdOptionExists(argv, argc + argv, "-r")){
			R = atof(getCmdOption(argv, argv + argc, "-r"));
		}
	}

	std::string image_prefix;
	if (scat_data_interpol_method == 0){
		scat_data_interpol = INTERPOL_METHOD::S1_G;
		image_prefix = model_name + "S1_G_";
	}
	else if (scat_data_interpol_method == 1){
		scat_data_interpol = INTERPOL_METHOD::S1_L;
		image_prefix = model_name + "S1_L_";
	}
	else if (scat_data_interpol_method == 2){
		scat_data_interpol = INTERPOL_METHOD::S2_G;
		image_prefix = model_name + "S2_G_";
	}
	else if (scat_data_interpol_method == 3){
		scat_data_interpol = INTERPOL_METHOD::S2_L;
		image_prefix = model_name + "S2_L_";
	}
	else if (scat_data_interpol_method == 4){
		scat_data_interpol = INTERPOL_METHOD::H_G_MQ;
		image_prefix = model_name + "H_G_MQ_";
	}
	else if (scat_data_interpol_method == 5){
		scat_data_interpol = INTERPOL_METHOD::H_L_MQ;
		image_prefix = model_name + "H_L_MQ_";
	}
	else if (scat_data_interpol_method == 6){
		scat_data_interpol = INTERPOL_METHOD::H_G_RE;
		image_prefix = model_name + "H_G_RE";
	}
	else if (scat_data_interpol_method == 7){
		scat_data_interpol = INTERPOL_METHOD::H_L_RE;
		image_prefix = model_name + "H_L_RE";
	}
	else{
		PRINT_ERROR("Invalid scatter data interpolation method. It should be [0,7]");
	}

	
	if (model_name.size() == 0){		
		model_name = image_prefix;
	}
	else
	{
		model_name += image_prefix;
	}
		
	
	
	std::cout << "	output_image_name= " << output_image_name << std::endl;
	std::cout << "	projection= " << projection << std::endl;	
	std::cout << "	image_resolution= " << image_res[0] << " X " << image_res[1] << std::endl;
	std::cout << "	scat_data_interpol_method= " << scat_data_interpol_method << std::endl;
	std::cout << "	num_data= " << num_data << std::endl;
	std::cout << "	k= " << K << std::endl;
	std::cout << "	R= " << R << std::endl;
	std::cout << "	slice_depth= " << slice_depth << std::endl;
	std::cout << "	data_fun= " << data_fun << std::endl;
	

	//setup the image using the right projection
	data_t image_x0[3], image_xn[3], image_normal[3];		
	data_t frac = 0.1 *(grid_upper[0] - grid_lower[0]);
	
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

		image_prefix += ((projection < 0) ? "-z" : "+z");
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

		image_prefix += ((projection < 0) ? "-y" : "+y");
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

		image_prefix += ((projection < 0) ? "-x" : "+x");
	}
	else {
		PRINT_ERROR("Invalid projection direction");
	}

	bool flip_vertical(false), flip_horizontal(false);


	//Init data 
	data_t data_f_value_min(0), data_f_value_max(0);
	ScatData<index_t, data_t> *my_data= NULL;

	//==================================
	//lambda function for the analytical function 
	auto data_analytical_func = [data_fun](data_t x, data_t y, data_t z){
		data_t f;
		if (data_fun == 1){
			f = sqrt(x*x + y*y + z*z);
		}
		else if (data_fun == 2){
			//f = sqrt(sin(x) + cos(y));
			f = cos(y)*sin(x);
		}
		return f;
	};
	//==================================

	if (bits == 8){
		init_data<index_t, data_t, unsigned char>(
			(data_fun==0) ? inputfilename : "" , my_data, n_grid,
			grid_lower, grid_upper,
			data_f_value_min,
			data_f_value_max, num_data, data_analytical_func);
	}
	else if (bits == 16){
		init_data<index_t, data_t, unsigned short>(
			(data_fun == 0) ? inputfilename : "", my_data, n_grid,
			grid_lower, grid_upper,
			data_f_value_min,
			data_f_value_max, num_data, data_analytical_func);
	}
	else{
		PRINT_ERROR("Invalid bits size" + std::to_string(bits));
	}

	//my_data->export_data("data.csv");
	
	
	

	//Init grid 
	Grid<index_t, data_t, DIM> *my_grid = NULL;
	data_t f_value_min(0), f_value_max(0);	

	std::chrono::high_resolution_clock::time_point t0_init_grid =
		std::chrono::high_resolution_clock::now();

	my_data->precompute(scat_data_interpol, K, R);

	if (bits == 8){
		init_grid<index_t, data_t, unsigned char>("", my_grid,
			f_value_min, f_value_max, n_grid, grid_lower, grid_upper, bg_color,
			light_color, my_data, scat_data_interpol, K, R);
	}
	else if (bits == 16){
		init_grid<index_t, data_t, unsigned short>("", my_grid,
			f_value_min, f_value_max, n_grid, grid_lower, grid_upper, bg_color,
			light_color, my_data, scat_data_interpol, K, R);
	}
	else{
		PRINT_ERROR("Invalid bits size" + std::to_string(bits));
	}
	

	std::chrono::high_resolution_clock::time_point t1_init_grid =
		std::chrono::high_resolution_clock::now();

	double init_grid_time = 
		std::chrono::duration<double, std::milli>(t1_init_grid - t0_init_grid).count();

	//image for the trilinear case
	TGAImage tga_image_linear(image_res[0], image_res[1], TGAImage::RGBA);
	Image<index_t, int_t, data_t> my_image_linear(image_res, image_x0, image_xn, 
		image_normal, &tga_image_linear, 
		((output_image_name.size() >= 1) ? output_image_name : image_prefix) + ".tga",
		flip_vertical, flip_horizontal);

	//init renderer   
	Renderer<index_t, int_t, data_t> my_renderer(my_grid, samples_per_cell);
	

	//Color Transfer Function 
	data_t stent_min(50000), stent_max(-50000);
	auto color_transfer_func = [&stent_min, &stent_max, f_value_min, f_value_max,
		model_name](data_t f_value){

		stent_min = std::min(stent_min, f_value);
		stent_max = std::max(stent_max, f_value);

		color_t color;
		COLOR_MAP cm = COLOR_MAP::PLASMA;
		colormap(cm,
			(f_value - f_value_min) / (1980 - f_value_min) //f_value / f_value_max
			,color.r, color.g, color.b);
		color.a = 1.0;		
		return color;
	};

	
	//do slicing 
	my_renderer.slice(&my_image_linear, color_transfer_func,
		data_analytical_func, slice_depth, true);


	std::cout << "stent_min= " << stent_min << std::endl;
	std::cout << "stent_max= " << stent_max << std::endl;

	/*std::fstream error_file("error.txt", std::ios::app);
	error_file.precision(30);
	error_file << model_name << " " << K << "	" << my_image_linear.get_error()
		<<"	"<< init_grid_time << std::endl;
	error_file.close();*/

	system("pause");
	

	return 0;
}
