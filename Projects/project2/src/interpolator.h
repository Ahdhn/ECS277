#ifndef __INTERPOLATOR__
#define __INTERPOLATOR__


template<typename T_d>
T_d interpolator_trilinear(const T_d sample[3], const T_d base_cell_corner[3],
	const T_d spacing [3], const T_d cell_f_value[8], T_d grad[3],
	bool skip_gradient = false){

	//Trilinear interpolate function at sample using function at six cell corners 	

	T_d x, y, z;
	x = (sample[0] - base_cell_corner[0]) / spacing[0];
	y = (sample[1] - base_cell_corner[1]) / spacing[1];
	z = (sample[2] - base_cell_corner[2]) / spacing[2];
		
	
	T_d c000 = cell_f_value[0];
	T_d c100 = cell_f_value[1];
	T_d c110 = cell_f_value[2];
	T_d c010 = cell_f_value[3];

	T_d c001 = cell_f_value[4];
	T_d c101 = cell_f_value[5];
	T_d c111 = cell_f_value[6];
	T_d c011 = cell_f_value[7];
	

	T_d f = (1 - x)*(1 - y)*(1 - z)*c000 +
		    (x)*(1 - y)*(1 - z)*c100 +
		    (1 - x)*(y)*(1 - z)*c010 +
		    (x)*(y)*(1 - z)*c110 +
		    (1 - x)*(1 - y)*(z)*c001 +
		    (x)*(1 - y)*(z)*c101 +
		    (1 - x)*(y)*(z)*c011 +
		    (x)*(y)*(z)*c111;

	if (!skip_gradient){
		grad[0] = (1 - y)*(1 - z)*(c100 - c000) +
			      (y)*(1 - z)*(c110 - c010) +
			      (1 - y)*(z)*(c101 - c001) +
			      (y)*(z)*(c111 - c011);

		grad[1] = (1 - x)*(1 - z)*(c010 - c000) +
			      (x)*(1 - z)*(c110 - c100) +
			      (1 - x)*(z)*(c011 - c001) +
			      (x)*(z)*(c111 - c101);

		grad[2] = (1 - x)*(1 - y)*(c001 - c000) +
			      (x)*(1 - y)*(c101 - c100) +
			      (1 - x)*(y)*(c011 - c010) +
			      (x)*(y)*(c111 - c110);
	}
	return f;
}

template<typename T_d, typename T>
T_d BB3(const T_d u, const T i){
	//Bernstein basis polynomials for n = 3
	//(3) (1-u)^(3-i) * u^(i) 
	//(i)
	T_d b_u = pow(1.0 - u, 3.0 - i)* pow(u, i);

	if (i == 0 || i == 3){
		return b_u;
	}
	else if (i == 1 || i == 2){
		return b_u*3.0;
	}
	else {
		PRINT_ERROR("BB3() invalid i");
	}
}
template<typename T_d, typename T>
T_d BB2(const T_d u, const T i){
	//Bernstein basis polynomials for n = 2
	//(2) (1-u)^(2-i) * u^(i) 
	//(i)
	T_d b_u = pow(1.0 - u, 2.0 - i)* pow(u, i);

	if (i == 0 || i == 2){
		return b_u;
	}
	else if (i == 1){
		return b_u*2.0;
	}
	else {
		PRINT_ERROR("BB2() invalid i");
	}
}
template<typename T_d>
T_d interpolator_tricubic(const T_d sample[3], const T_d base_cell_corner[3],
	const T_d spacing[3], const T_d cell_f_value[8], 
	const T_d cell_approx_g[8][3], T_d grad[3],
	bool skip_gradient = false){

	//Tricubic interpolate function at sample using function at six cell corners 	

	T_d x, y, z;
	x = (sample[0] - base_cell_corner[0]) / spacing[0];
	y = (sample[1] - base_cell_corner[1]) / spacing[1];
	z = (sample[2] - base_cell_corner[2]) / spacing[2];


	T_d c000 = cell_f_value[0];
	T_d c100 = cell_f_value[1];
	T_d c110 = cell_f_value[2];
	T_d c010 = cell_f_value[3];

	T_d c001 = cell_f_value[4];
	T_d c101 = cell_f_value[5];
	T_d c111 = cell_f_value[6];
	T_d c011 = cell_f_value[7];

	T_d g000[3] = { cell_approx_g[0][0], cell_approx_g[0][1], cell_approx_g[0][2] };
	T_d g100[3] = { cell_approx_g[1][0], cell_approx_g[1][1], cell_approx_g[1][2] };
	T_d g110[3] = { cell_approx_g[2][0], cell_approx_g[2][1], cell_approx_g[2][2] };
	T_d g010[3] = { cell_approx_g[3][0], cell_approx_g[3][1], cell_approx_g[3][2] };

	T_d g001[3] = { cell_approx_g[4][0], cell_approx_g[4][1], cell_approx_g[4][2] };
	T_d g101[3] = { cell_approx_g[5][0], cell_approx_g[5][1], cell_approx_g[5][2] };
	T_d g111[3] = { cell_approx_g[6][0], cell_approx_g[6][1], cell_approx_g[6][2] };
	T_d g011[3] = { cell_approx_g[7][0], cell_approx_g[7][1], cell_approx_g[7][2] };


	T_d b[4][4][4];
	//Ever b_{i,j,k} can be written as 
	//F_{i,j,k} +/- 0.3*GX_{i,j,k} +/- 0.3*GY_{i,j,k} +/- 0.3*GZ_{i,j,k}

	//The three indices of F can be derived such that is dimesion D is 0 or 1, 
	//then the corresponding index is 0 otherwise, the corresponding index is 1
	//e.g. b000 -> F000, b301 -> F100, b222-> F111

	//if the dimension D is 0 or 3, then the corrsponding G will be eliminated 
	//e.g. b000 = F000 (no G's), b101 = F000 + 0.3*GX000 + 0.3GZ000, 
	//b222 = F111 - 0.3GX111 - 0.3GY111 - 0.3GZ111, b333 = F111

	//The indices of the G's follows the same rule as the indices of F. The sign 
	//of the 0.3 is derived such that if the index is 0 or 1, it is +ve, otherwise
	//it is -ve 
	double X = 1.0 / 3.0;

	//z = 0
	b[0][0][0] = c000;//
	b[1][0][0] = c000 + X*g000[0];//+x
	b[2][0][0] = c100 - X*g100[0];//-x
	b[3][0][0] = c100;//

	b[0][1][0] = c000 + X*g000[1];//+y
	b[1][1][0] = c000 + X*g000[0] + X*g000[1];//+x +y 
	b[2][1][0] = c100 - X*g100[0] + X*g100[1];//-x +y
	b[3][1][0] = c100 + X*g100[1];//+y

	b[0][2][0] = c010 - X*g010[1];//-y
	b[1][2][0] = c010 + X*g010[0] - X*g010[1];//+x -y
	b[2][2][0] = c110 - X*g110[0] - X*g110[1];//-x -y
	b[3][2][0] = c110 - X*g110[1];//-y

	b[0][3][0] = c010;//
	b[1][3][0] = c010 + X*g010[0];//+x 
	b[2][3][0] = c110 - X*g110[0];//-x
	b[3][3][0] = c110;//

	//z = 1
	b[0][0][1] = c000 + X*g000[2];//+z
	b[1][0][1] = c000 + X*g000[0] + X*g000[2];//+x +z
	b[2][0][1] = c100 - X*g100[0] + X*g100[2];//-x +z
	b[3][0][1] = c100 + X*g100[2];//+z

	b[0][1][1] = c000 + X*g000[1] + X*g000[2];//+y +z
	b[1][1][1] = c000 + X*g000[0] + X*g000[1] + X*g000[2];//+x +y +z
	b[2][1][1] = c100 - X*g100[0] - X*g100[1] + X*g100[2];//-x -y +z
	b[3][1][1] = c100 - X*g100[1] + X*g100[2];//-y +z

	b[0][2][1] = c010 - X*g010[1] + X*g010[2];//-y +z
	b[1][2][1] = c010 + X*g010[0] - X*g010[1] + X*g010[2];//+x -y +z
	b[2][2][1] = c110 - X*g110[0] - X*g110[1] + X*g110[2];//-x -y +z
	b[3][2][1] = c110 - X*g110[1] + X*g110[2];//-y +z

	b[0][3][1] = c010 + X*g010[2];//+z
	b[1][3][1] = c010 + X*g010[0] + X*g010[2];//+x +z
	b[2][3][1] = c110 - X*g110[0] + X*g110[2];//-x +z
	b[3][3][1] = c110 + X*g110[2];//+z

	//z = 2
	b[0][0][2] = c001 - X*g001[2];//-z
	b[1][0][2] = c001 + X*g001[0] - X*g001[2];//+x -z
	b[2][0][2] = c101 - X*g101[0] - X*g101[2];//-x -z
	b[3][0][2] = c101 - X*g101[2];//-z

	b[0][1][2] = c001 + X*g001[1] - X*g001[2];//+y -z
	b[1][1][2] = c001 + X*g001[0] + X*g001[1] - X*g001[2];//+x +y -z
	b[2][1][2] = c101 - X*g101[0] + X*g101[1] - X*g101[2];//-x +y -z
	b[3][1][2] = c101 + X*g101[1] - X*g101[2]; //+y -z

	b[0][2][2] = c011 - X*g011[1] - X*g011[2];//-y -z
	b[1][2][2] = c011 + X*g011[0] - X*g011[1] - X*g011[2];//+x -y -z
	b[2][2][2] = c111 - X*g111[0] - X*g111[1] - X*g111[2];//-x -y -z
	b[3][2][2] = c111 - X*g111[1] - X*g111[2];//-y -z

	b[0][3][2] = c011 - X*g011[2];//-z
	b[1][3][2] = c011 + X*g011[0] - X*g011[2];//+x -z
	b[2][3][2] = c111 - X*g111[0] - X*g111[2];//-x -z
	b[3][3][2] = c111 - X*g111[2];//-z


	//z = 3
	b[0][0][3] = c001;//
	b[1][0][3] = c001 + X*g001[0];//+x
	b[2][0][3] = c101 - X*g101[0];//-x 
	b[3][0][3] = c101;//

	b[0][1][3] = c001 + X*g001[1];//+y
	b[1][1][3] = c001 + X*g001[0] + X*g001[1];//+x +y
	b[2][1][3] = c101 - X*g101[0] + X*g101[1];//-x +y
	b[3][1][3] = c101 + X*g101[1];//+y

	b[0][2][3] = c011 - X*g011[1];//-y
	b[1][2][3] = c011 + X*g011[0] - X*g011[1];//+x -y
	b[2][2][3] = c111 - X*g111[0] - X*g111[1];//-x -y
	b[3][2][3] = c111 - X*g111[1];//-y

	b[0][3][3] = c011;//
	b[1][3][3] = c011 + X*g011[0];//+x
	b[2][3][3] = c111 - X*g111[0];//-x
	b[3][3][3] = c111;//
	


	//computation of the function 
	T_d f = 0;
	for (uint32_t i = 0; i < 4; i++){
		for (uint32_t j = 0; j < 4; j++){
			for (uint32_t k = 0; k < 4; k++){
				f += b[i][j][k] * BB3(x, i) *  BB3(y, j) * BB3(z, k);
			}
		}
	}

	if (!skip_gradient){
		//computation for the gradient 	
		//x
		grad[0] = 0;
		for (uint32_t i = 0; i < 3; i++){
			for (uint32_t j = 0; j < 4; j++){
				for (uint32_t k = 0; k < 4; k++){
					grad[0] += (b[i + 1][j][k] - b[i][j][k])*BB2(x, i)*BB3(y, j)*BB3(z, k);
				}
			}
		}
		grad[0] *= 3;

		//y
		grad[1] = 0;
		for (uint32_t i = 0; i < 4; i++){
			for (uint32_t j = 0; j < 3; j++){
				for (uint32_t k = 0; k < 4; k++){
					grad[1] += (b[i][j + 1][k] - b[i][j][k])*BB3(x, i)*BB2(y, j)*BB3(z, k);
				}
			}
		}
		grad[1] *= 3;

		//z
		grad[2] = 0;
		for (uint32_t i = 0; i < 4; i++){
			for (uint32_t j = 0; j < 4; j++){
				for (uint32_t k = 0; k < 3; k++){
					grad[2] += (b[i][j][k + 1] - b[i][j][k])*BB3(x, i)*BB3(y, j)*BB2(z, k);
				}
			}
		}
		grad[2] *= 3;
	}
	
	return f;

}
#endif /*__INTERPOLATOR__*/