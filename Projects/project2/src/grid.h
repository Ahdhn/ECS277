#ifndef __GRID__
#define __GRID__

#include <stdint.h>
#include <cstdlib>
#include <stdio.h>
#include <algorithm>
#include <numeric>

#include "interpolator.h"


template <class T, class T_d, uint32_t DIM>
class Grid
{

public:	
	Grid(){};
	Grid(T n[3], T_d x0[3], T_d s[3], color_t background_color, 
		color_t light_color, T_d light_loc[3], T_d light_dir[3],T_d phong_power,
		T_d attenuation_const, T_d light_intensity, color_t ambient_light,
		T_d ambient_const, T_d diffuse_const, T_d specular_const);

	void fill_grid(T i, T j, T k, T_d data){		
		T my_flat_id = flat_id(i, j, k);
		m_data[my_flat_id] = data;
	};

	void fill_grid(T data_flat_id, T_d data){		
		if (data_flat_id >= m_grid_size){
			PRINT_ERROR("Grid::fill_grid() id is out of range!!!");
		}
		m_data[data_flat_id] = data;
	}
		

	T_d get_value(T i, T j, T k){
		T my_flat_id = flat_id(i,j,k);
		return m_data[my_flat_id];
	}

	bool get_ray_grid_intersect(T_d ray_org[3], T_d ray_dir[3], T_d seg_start[3],
		T_d seg_end[3]);

	bool is_inside_grid(T_d sample[3]){
		if (sample[0] + EPSILON > m_x1[0] ||		
			sample[1] + EPSILON > m_x1[1] ||
			sample[2] + EPSILON > m_x1[2] ){
			return false;
		}
		return true;
	}

	T_d get_f_value_at_sample(T_d sample[3], T_d grad[3], INTERPOL_TYPE 
		type = INTERPOL_TYPE::TRILINEAR);
	
	color_t get_background_color(){
		return m_background_color;
	}

	T get_resolution(T dir){
		return m_n[dir];
	}

	//shade a point in the volum based of Phone model 	
	color_t shade_point(const T_d sample[3], const T_d dist_sq,
		const color_t local_color, const T_d ray_dir[3], T_d N[3]);

	//export the grid as csv file 
	void export_grid(std::string filename);

	~Grid(){};

private:
	inline T flat_id(T i, T j, T k);
	bool ray_grid_plane_intersect(T plane_id, T_d ray_org[3], T_d ray_dir[3],
		T_d point[3]);
	inline void get_grid_plane(T plane_num, T plane[3]);
	inline void gird_corners_coordinates(T corner_num, T_d coord[3]);

	T_d*m_data = nullptr; //function data at each point 
	T m_n[DIM]; //number of points in each dim
	T m_grid_size;
	
	T_d m_s[DIM]; //spaceing

	T_d m_x0[DIM]; //lower most corner 
	T_d m_x1[DIM]; //upper most corner 
	T_d m_cell_f_value[8];//use this to pass cell function values to interpolator	
	T m_base_cell_index[3];
	T_d m_base_cell_corner[3];
	//approximate gradient from the function value via finite difference 
	//used with tricubic interpolation 
	T_d m_cell_approx_grad[8][3];


	color_t m_background_color;
	color_t m_light_color;
	color_t m_ambient_light;

	T_d m_light_loc[3];
	T_d m_light_dir[3];

	T_d m_phong_power;
	T_d m_attenuation_const;
	T_d m_light_intensity;
	T_d m_ambient_const;
	T_d m_diffuse_const;
	T_d m_specular_const;

	

	struct plane {
		T corners[4];
		T_d corners_coord[4][3];
		T_d normal[3];
		T_d lower[3];//min x,y,z
		T_d upper[3];//max x,y,z
	}grid_planes[6];
};

template<class T, class T_d,  uint32_t DIM>
Grid<T, T_d, DIM>::Grid(T n[3], T_d x0[3], T_d s[3], color_t background_color,
	color_t light_color, T_d light_loc[3], T_d light_dir[3], T_d phong_power,
	T_d attenuation_const, T_d light_intensity, color_t ambient_light,
	T_d ambient_const, T_d diffuse_const, T_d specular_const) :
	m_background_color(background_color), m_light_color(light_color),
	m_phong_power(phong_power), m_attenuation_const(attenuation_const),
	m_light_intensity(light_intensity), m_ambient_light(ambient_light),
	m_ambient_const(ambient_const),
	m_diffuse_const(diffuse_const), m_specular_const(specular_const){

	//TODO change this so we can take 2D and 3D grid 
	if (DIM != 3){
		printf("\n Grid::Grid DIM shoudl be 3!!!");
		exit(EXIT_FAILURE);
	}
	m_grid_size = 1;
	for (uint32_t i = 0; i < DIM; i++){
		m_n[i] = n[i];
		m_x0[i] = x0[i];
		m_s[i] = s[i];
		m_grid_size *= m_n[i];
		m_x1[i] = m_x0[i] + m_s[i] * (m_n[i] - 1);
	}
	
	T_d normal[3];
	for (uint32_t f = 0; f < 6; f++){
		get_grid_plane(f, grid_planes[f].corners);
		gird_corners_coordinates(grid_planes[f].corners[0], 
			grid_planes[f].corners_coord[0]);
		gird_corners_coordinates(grid_planes[f].corners[1], 
			grid_planes[f].corners_coord[1]);
		gird_corners_coordinates(grid_planes[f].corners[2], 
			grid_planes[f].corners_coord[2]);
		gird_corners_coordinates(grid_planes[f].corners[3],
			grid_planes[f].corners_coord[3]);

		
		for (uint32_t i = 0; i < 3; i++){
			grid_planes[f].lower[i] = std::numeric_limits<T_d>::max();
			grid_planes[f].upper[i] = -grid_planes[f].lower[i];
			for (uint32_t c = 0; c < 4; c++){
				grid_planes[f].lower[i] = std::min(grid_planes[f].lower[i], 
					grid_planes[f].corners_coord[c][i]);
				grid_planes[f].upper[i] = std::max(grid_planes[f].upper[i],
					grid_planes[f].corners_coord[c][i]);
			}
			
		}

		cross_product(grid_planes[f].corners_coord[0][0] - grid_planes[f].corners_coord[1][0],
			          grid_planes[f].corners_coord[0][1] - grid_planes[f].corners_coord[1][1],
					  grid_planes[f].corners_coord[0][2] - grid_planes[f].corners_coord[1][2],
					  grid_planes[f].corners_coord[2][0] - grid_planes[f].corners_coord[1][0],
					  grid_planes[f].corners_coord[2][1] - grid_planes[f].corners_coord[1][1],
					  grid_planes[f].corners_coord[2][2] - grid_planes[f].corners_coord[1][2],
					  grid_planes[f].normal);
		
		normalize_vector(grid_planes[f].normal[0], grid_planes[f].normal[1],
			grid_planes[f].normal[2]);		

	}


	m_light_loc[0] = light_loc[0];
	m_light_loc[1] = light_loc[1];
	m_light_loc[2] = light_loc[2];

	m_light_dir[0] = light_dir[0];
	m_light_dir[1] = light_dir[1];
	m_light_dir[2] = light_dir[2];

	

	m_data = (T_d*)malloc(m_grid_size*sizeof(T_d));

	m_background_color /= 255.0;
	m_light_color /= 255.0;
	m_ambient_light /= 255.0;

}

template<class T, class T_d, uint32_t DIM>
inline T Grid<T,T_d,DIM>::flat_id(T i, T j, T k){

	//T my_flat = i + m_n[1] * (j + m_n[2] * k);
	

	T my_flat = i + j*m_n[0] + k*m_n[0] * m_n[1]; //(i,j,k)
	//T my_flat = j + i*m_n[1] + k*m_n[0] * m_n[1]; //(j,i,k)
	//T my_flat = k + i*m_n[2] + j*m_n[0] * m_n[2]; //(k,i,j)
	//T my_flat = i + k*m_n[0] + j*m_n[0] * m_n[2]; //(i,k,j)
	//T my_flat = k + j*m_n[2] + i*m_n[1] * m_n[2]; //(k,j,i)
	//T my_flat = j + k*m_n[1] + i*m_n[1] * m_n[2]; //(j,k,i)
	
	
	if (my_flat >= m_grid_size){
		PRINT_ERROR("Grid::flat_id() invalid flat id with (" + std::to_string(i)
			+  ", " + std::to_string(j)+ ", " + std::to_string(k)+") gives " + 
			std::to_string(my_flat) + " where grid size is "+
			std::to_string(m_grid_size));
	}

	return my_flat;
}


template<class T, class T_d, uint32_t DIM>
inline void Grid<T, T_d, DIM>::get_grid_plane(T plane_num, T plane[4]){
	//plane are returned in ccw orientation when look at it from outside 
	if (plane_num == 0){
		plane[0] = 0;
		plane[1] = 1;
		plane[2] = 2;
		plane[3] = 3;
	}
	else if (plane_num == 1){
		plane[0] = 6;
		plane[1] = 5;
		plane[2] = 4;
		plane[3] = 7;
	}
	else if (plane_num == 2){
		plane[0] = 1;
		plane[1] = 5;
		plane[2] = 6;
		plane[3] = 2;
	}
	else if (plane_num == 3){
		plane[0] = 0;
		plane[1] = 3;
		plane[2] = 7;
		plane[3] = 4;
	}
	else if (plane_num == 4){
		plane[0] = 3;
		plane[1] = 2;
		plane[2] = 6;
		plane[3] = 7;
	}
	else if (plane_num == 5){
		plane[0] = 0;
		plane[1] = 4;
		plane[2] = 5;
		plane[3] = 1;
	}
	else{
		PRINT_ERROR("Grid::grid_plane plane_num should be [0,5]");
	}
}
template<class T, class T_d, uint32_t DIM>
inline void Grid<T, T_d, DIM>::gird_corners_coordinates(T corner_num,
	T_d coord[3]){
	if (corner_num == 0 || corner_num == 1 || corner_num == 2 || corner_num == 3){
		coord[2] = m_x0[2];
	}
	else{
		coord[2] = m_x1[2];
	}

	if (corner_num == 0 || corner_num == 1 || corner_num == 4 || corner_num == 5){
		coord[1] = m_x0[1];
	}
	else{
		coord[1] = m_x1[1];
	}

	if (corner_num == 0 || corner_num == 3 || corner_num == 4 || corner_num == 7){
		coord[0] = m_x0[0];
	}
	else{
		coord[0] = m_x1[0];
	}



}



template<class T, class T_d, uint32_t DIM>
bool Grid<T, T_d, DIM>::get_ray_grid_intersect(T_d ray_org[3], T_d ray_dir[3], 
	T_d seg_start[3], T_d seg_end[3]){
	//TODO intersect the ray with the grid six faces to get the start and 
	//end of the segment that penetrates the grid 
			
	uint32_t num_points = 0;
	T_d point[3];
	for (uint32_t f = 0; f < 6; f++){

		if (ray_grid_plane_intersect(f, ray_org, ray_dir, point)){
			if (num_points == 0){
				seg_start[0] = point[0];
				seg_start[1] = point[1];
				seg_start[2] = point[2];
				num_points++;
			}
			else if (num_points == 1){
				seg_end[0] = point[0];
				seg_end[1] = point[1];
				seg_end[2] = point[2];
				num_points++;
			}
			else {
				PRINT_ERROR("Grid::get_ray_gird_intersect() more than two intersection points with a plane!!!");
			}
		}
	}

	if (num_points != 0 && num_points != 2){
		PRINT_ERROR("Grid::get_ray_gird_intersect() ray-plane intersection went wrong!!!");
	}

	if (num_points == 0){
		return false;
	}

	//order the points 
	T_d d0 = Dist(seg_start[0], seg_start[1], seg_start[2], 
		ray_org[0], ray_org[1], ray_org[2]);
	T_d d1 = Dist(seg_end[0], seg_end[1], seg_end[2], 
		ray_org[0], ray_org[1], ray_org[2]);

	if (d1 < d0){
		std::swap(seg_start[0], seg_end[0]);
		std::swap(seg_start[1], seg_end[1]);
		std::swap(seg_start[2], seg_end[2]);
	}
	return true;
}


template<class T, class T_d, uint32_t DIM>
bool Grid<T, T_d, DIM>::ray_grid_plane_intersect(T plane_id, T_d ray_org[3],
	T_d ray_dir[3], T_d point[3]){
	
	if (plane_id > 5 || plane_id < 0){
		PRINT_ERROR("rid::ray_grid_plane_intersect plane_id should be in [0,5]");
	}
	
	if (!ray_plane_intersect<T_d>(
		grid_planes[plane_id].corners_coord[0][0],
		grid_planes[plane_id].corners_coord[0][1], 
		grid_planes[plane_id].corners_coord[0][2],
		grid_planes[plane_id].normal[0], grid_planes[plane_id].normal[1],
		grid_planes[plane_id].normal[2],			
		ray_dir[0], ray_dir[1], ray_dir[2], ray_org[0], ray_org[1], ray_org[2],
		point[0], point[1], point[2])){
		return false;
	}
	//make sure the point is actually inside the gird face 
	//TODO generalize this incase the grid is not axis-aligned 

	if (point[0] < grid_planes[plane_id].lower[0] ||
		point[0] > grid_planes[plane_id].upper[0] ||
		point[1] < grid_planes[plane_id].lower[1] ||
		point[1] > grid_planes[plane_id].upper[1] ||
		point[2] < grid_planes[plane_id].lower[2] ||
		point[2] > grid_planes[plane_id].upper[2]){
		return false;
	}

	return true;
	

}

template<class T, class T_d, uint32_t DIM>
T_d Grid<T, T_d, DIM>::get_f_value_at_sample(T_d sample[3], T_d grad[3],
	INTERPOL_TYPE type /*= INTERPOL_TYPE::TRILINEAR*/){

	//use interpolator to get this function

	//get the cell that contains the sample 
	
	for (T s = 0; s < 3; s++){
		m_base_cell_index[s] = T((sample[s] - m_x0[s]) / m_s[s]);		
		m_base_cell_corner[s] = m_x0[s] + m_base_cell_index[s] * m_s[s];
	}

	uint32_t i, j, k;
	for (uint32_t c = 0; c < 8; c++){
		if (c == 0){
			i = m_base_cell_index[0];
			j = m_base_cell_index[1];
			k = m_base_cell_index[2];			
		}
		else if (c == 1){
			i = m_base_cell_index[0] + 1;
			j = m_base_cell_index[1];
			k = m_base_cell_index[2];
		}
		else if (c == 2){
			i = m_base_cell_index[0] + 1;
			j = m_base_cell_index[1] + 1;
			k = m_base_cell_index[2];			
		}
		else if (c == 3){
			i = m_base_cell_index[0];
			j = m_base_cell_index[1] + 1;
			k = m_base_cell_index[2];			
		}
		else if (c == 4){
			i = m_base_cell_index[0];
			j = m_base_cell_index[1];
			k = m_base_cell_index[2] + 1;			
		}
		else if (c == 5){
			i = m_base_cell_index[0] + 1;
			j = m_base_cell_index[1];
			k = m_base_cell_index[2] + 1;			
		}
		else if (c == 6){
			i = m_base_cell_index[0] + 1;
			j = m_base_cell_index[1] + 1;
			k = m_base_cell_index[2] + 1;			
		}
		else if (c == 7){
			i = m_base_cell_index[0];
			j = m_base_cell_index[1] + 1;
			k = m_base_cell_index[2] + 1;			
		}
		T my_flat_id = flat_id(i, j, k);
		m_cell_f_value[c] = m_data[my_flat_id];

		if (type == INTERPOL_TYPE::TRICUBIC){
			uint32_t i_p, i_m, j_p, j_m, k_p, k_m;

			i_p = (i == m_n[0] - 1) ? i : i + 1;
			i_m = (i == 0) ? i : i - 1;

			j_p = (j == m_n[1] - 1) ? j : j + 1;
			j_m = (j == 0) ? j : j - 1;

			k_p = (k == m_n[2] - 1) ? k : k + 1;
			k_m = (k == 0) ? k : k - 1;
									
			m_cell_approx_grad[c][0] = 0.5*(get_value(i_p, j, k) -
				get_value(i_m, j, k));
			m_cell_approx_grad[c][1] = 0.5*(get_value(i, j_p, k) -
				get_value(i, j_m, k));
			m_cell_approx_grad[c][2] = 0.5*(get_value(i, j, k_p) -
				get_value(i, j, k_m));
		}
	}

	if (type == INTERPOL_TYPE::TRILINEAR){
		return interpolator_trilinear(sample, m_base_cell_corner, m_s,
			m_cell_f_value, grad);
	}
	else if (type == INTERPOL_TYPE::TRICUBIC){
		//compute the approximated graident 
		return interpolator_tricubic(sample, m_base_cell_corner, m_s,
			m_cell_f_value, m_cell_approx_grad, grad);

	}
	else{
		PRINT_ERROR("Grid::get_f_value_at_sample() unknown interpolation type");
	}

}

template <class T, class T_d, uint32_t DIM>
void Grid<T, T_d, DIM>::export_grid(std::string filename){
	std::fstream file(filename, std::ios::out);
	file.precision(12);
	file << "x coord,y coord,z coord,intensity" << std::endl;

	for (uint32_t i = 50; i < 100; i++){
		T_d xx = m_x0[0] + i*m_s[0];
		for (uint32_t j = 150; j < m_n[1]; j++){
			T_d yy = m_x0[1] + j*m_s[1];
			for (uint32_t k = 0; k < m_n[2]; k++){
				T_d zz = m_x0[2] + k*m_s[2];
				
				file << xx << ", " << yy << ", " << zz << ", "
					<< get_value(i, j, k) << std::endl;				
			}
		}
	}

	file.close();
}

template <class T, class T_d, uint32_t DIM>
color_t Grid<T, T_d, DIM>::shade_point(const T_d sample[3], const T_d dist_sq,
	const color_t local_color, const T_d ray_dir[3], T_d N[3]){
	//dist_to_ff is distance to the eye

	T_d L[3]; //vector from sample to the light position
	L[0] = m_light_loc[0] - sample[0];
	L[1] = m_light_loc[1] - sample[1];
	L[2] = m_light_loc[2] - sample[2];
	normalize_vector(L);

	T_d V[3]; //vector from sample to the eye
	V[0] = -ray_dir[0];
	V[1] = -ray_dir[1];
	V[2] = -ray_dir[2];
	normalize_vector(V);

	if (abs(N[0]) < EPSILON && abs(N[1]) < EPSILON && abs(N[2]) < EPSILON){
		//constant gradient at the sample will lead to inaccurate calc
		//thus we return the diffuse color 
		color_t col;
		T_d scale = 0.1;
		col = local_color;
		col *= scale;
		return col;
	}

	//vector normal to the surface that sample lies on 
	normalize_vector(N);	

	//https://math.stackexchange.com/a/13263
	T_d R[3];
	T_d LN_dot = Dot(-L[0], -L[1], -L[2], N[0], N[1], N[2]);	
	R[0] = -L[0] - 2 * LN_dot * N[0];
	R[1] = -L[1] - 2 * LN_dot * N[1];
	R[2] = -L[2] - 2 * LN_dot * N[2];
	normalize_vector(R);
	

	//attemp 1
	/*color_t my_color_ambient;	
	my_color_ambient = m_ambient_light;
	my_color_ambient *= m_ambient_const;

	color_t my_color_diffuse = local_color;
	my_color_diffuse *= m_diffuse_const*std::max(Dot(N, L), 0.0);

	color_t my_color_specular = m_light_color;	
	T_d specular = std::max(Dot(R, V), 0.0);
	my_color_specular *= m_specular_const * pow(specular, m_phong_power);
	
	color_t my_color;
	my_color = my_color_ambient;
	my_color += my_color_diffuse;
	my_color += my_color_specular;
	my_color.a = 1.0;*/

	
	//attemp 2
	color_t my_color;
	my_color = local_color;
	my_color *= m_ambient_light;
	my_color *= m_ambient_const;

	T_d diffuse = m_diffuse_const*std::max(Dot(N, L), 0.0);
	T_d specular;
	if (diffuse < 0.01){
		specular = 0.0f;
	}	
	else{
		specular = m_specular_const*std::max(Dot(R, V), 0.0);
		specular = pow(specular, m_phong_power);
	}

	T_d attenuation = 1.0 / (sqrt(dist_sq) + m_attenuation_const);

	color_t color_diff_spec = local_color;
	color_diff_spec *= m_light_color;
	color_diff_spec *= attenuation;
	color_diff_spec *= (diffuse + specular);
	my_color += color_diff_spec;

	
	T_d dif = std::max(Dot(N, L), 0.0);	
	color_t my_diff = local_color;
	my_diff *= dif;

	T_d spec = std::max(Dot(V, R), 0.0);
	spec = pow(spec, m_phong_power);
	color_t my_spec = m_light_color;
	my_spec *= spec;
		

	my_color += my_spec;
	my_color += my_diff;
	my_color *= m_light_intensity / (dist_sq*m_attenuation_const + 1);
	
	
	return my_color;
}


#endif /*__GRID__*/