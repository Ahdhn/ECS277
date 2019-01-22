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
	Grid(T n[3], T_d x0[3], T_d s[3], color_t background_color);
	void fill_grid(T i, T j, T k, T_d data){		
		T my_flat_id = flat_id(i, j, k);
		m_data[my_flat_id] = data;
	};
	
	T_d get_value(T i, T j, T k){
		T my_flat_id = flat_id(i,j,k);
		return m_data[my_flat_id];
	}

	bool get_ray_grid_intersect(T_d ray_org[3], T_d ray_dir[3], T_d seg_start[3],
		T_d seg_end[3]);

	T_d get_f_value_at_sample(T_d sample[3]);

	~Grid(){
		
	};

private:
	inline T flat_id(T i, T j, T k);
	bool ray_grid_plane_intersect(T plane_id, T_d ray_org[3], T_d ray_dir[3],
		T_d point[3]);
	inline void get_grid_plane(T plane_num, T plane[3]);
	inline void gird_corners_coordinates(T corner_num, T_d coord[3]);

	T_d*m_data = nullptr; //function data at each point 
	T m_n[DIM]; //number of points in each dim
	T grid_size;
	
	T_d m_s[DIM]; //spaceing

	T_d m_x0[DIM]; //lower most corner 
	T_d m_x1[DIM]; //upper most corner 
	T_d m_cell_f_value[6];//use this to pass cell function values to interpolator
	T m_base_cell_index[3];
	T_d m_base_cell_corner[3];

	color_t m_background_color;

	struct plane {
		T corners[4];
		T_d corners_coord[4][3];
		T_d normal[3];
		T_d lower[3];//min x,y,z
		T_d upper[3];//max x,y,z
	}grid_planes[6];
};

template<class T, class T_d,  uint32_t DIM>
Grid<T, T_d, DIM>::Grid(T n[3], T_d x0[3], T_d s[3], color_t background_color):
	m_background_color(background_color)
{
	//TODO change this so we can take 2D and 3D grid 
	if (DIM != 3){
		printf("\n Grid::Grid DIM shoudl be 3!!!");
		exit(EXIT_FAILURE);
	}
	grid_size = 1;
	for (uint32_t i = 0; i < DIM; i++){
		m_n[i] = n[i];
		m_x0[i] = x0[i];
		m_s[i] = s[i];
		grid_size *= m_n[i];
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

	m_data = (T_d*)malloc(grid_size*sizeof(T_d));
}

template<class T, class T_d, uint32_t DIM>
inline T Grid<T,T_d,DIM>::flat_id(T i, T j, T k){
	return i + m_n[1] * (j + m_n[2] * k);
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
T_d Grid<T, T_d, DIM>::get_f_value_at_sample(T_d sample[3]){
	//use interpolator to get this function

	//get the cell that contains the sample 
	

	for (T i = 0; i < 3; i++){
		m_base_cell_index[i] = T((sample[i] - m_x0[i]) / m_s[i]);		
		m_base_cell_corner[i] = m_x0[i] + m_base_cell_index[i] * m_s[i];
	}

	//c=0 --> m_base_cell_index[0], m_base_cell_index[1], m_base_cell_index[2]
	//c=1 --> m_base_cell_index[0] + 1, m_base_cell_index[1], m_base_cell_index[2]
	//c=2 --> m_base_cell_index[0], m_base_cell_index[1] + 1, m_base_cell_index[2]
	//c=3 --> m_base_cell_index[0] + 1, m_base_cell_index[1] + 1, m_base_cell_index[2]
	//c=4 --> m_base_cell_index[0], m_base_cell_index[1], m_base_cell_index[2] + 1
	//c=5 --> m_base_cell_index[0] + 1, m_base_cell_index[1], m_base_cell_index[2] + 1
	//c=6 --> m_base_cell_index[0], m_base_cell_index[1] + 1, m_base_cell_index[2] + 1
	//c=7 --> m_base_cell_index[0] + 1, m_base_cell_index[1] + 1, m_base_cell_index[2] + 1

	T i, j, k;
	for (uint32_t c = 0; c < 6; c++){
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
			i = m_base_cell_index[0];
			j = m_base_cell_index[1] + 1;
			k = m_base_cell_index[2];
		}
		else if (c == 3){
			i = m_base_cell_index[0] + 1;
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
			i = m_base_cell_index[0];
			j = m_base_cell_index[1] + 1;
			k = m_base_cell_index[2] + 1;
		}
		else if (c == 7){
			i = m_base_cell_index[0] + 1;
			j = m_base_cell_index[1] + 1;
			k = m_base_cell_index[2] + 1;
		}
		T my_flat_id = flat_id(i, j, k);
		m_cell_f_value[i] = m_data[my_flat_id];
	}

	return interpolator_trilinear(sample, m_base_cell_corner, m_s, m_cell_f_value);

}


#endif /*__GRID__*/