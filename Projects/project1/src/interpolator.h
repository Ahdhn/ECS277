#ifndef __INTERPOLATOR__
#define __INTERPOLATOR__


template<typename T_d>
T_d interpolator_trilinear(const T_d sample[3], const T_d base_cell_corner[3],
	const T_d spacing [3], const T_d cell_f_value[6]){
	//interpolate function at sample using function at six cell corners 
	
	T d_x, d_y, d_z;
	d_x = sample[0] - base_cell_corner[0];
	d_y = sample[1] - base_cell_corner[1];
	d_z = sample[2] - base_cell_corner[2];

	T_d f = (1 - d_x)*(1 - d_y)*(1 - d_z)cell_f_value[0] + 
		    (d_x)*(1 - d_y)*(1 - d_z)cell_f_value[1] +
			(1 - d_x)*(d_y)*(1 - d_z)cell_f_value[2] +
			(1 - d_x)*(1 - d_y)*(d_z)cell_f_value[3] +
			(d_x)*(d_y)*(1 - d_z)cell_f_value[3] +
}


#endif /*__INTERPOLATOR__*/