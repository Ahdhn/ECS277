#ifndef __COMMON__
#define __COMMON__

#define EPSILON 10E-6

enum class INTERPOL_TYPE{
	TRILINEAR = 1,
	TRICUBIC = 2
};

template<typename T_d>
struct scat_data_t
{
	T_d x, y, z, f;	
};

enum class INTERPOL_METHOD {
	S1 = 0,
	S2 = 1,
	S3 = 2,
	H = 3,
};

struct kd_node_t{
	double*x;
	int myID;//the id as in the passed point cloud 
	struct kd_node_t *left, *right;
	double dist;//distance to quary point (used in FindNNearest to avoid redundant computation)
};

struct color_t {
	//the color range is from 0 to 255	
	double r, g, b, a;

	void clamp(){
		clamp_comp(r, 1.0);
		clamp_comp(g, 1.0);
		clamp_comp(b, 1.0);
		clamp_comp(a, 1.0);
	}

	color_t operator=(const color_t&copy){		
		this->r = copy.r;
		this->g = copy.g;
		this->b = copy.b;
		this->a = copy.a;
		return *this;
	}


	color_t operator+=(color_t&c1){
		this->r += c1.r;
		this->g += c1.g;
		this->b += c1.b;
		this->a += c1.a;
		return *this;
	}
	
	color_t operator*=(color_t&c1){
		this->r *= c1.r;
		this->g *= c1.g;
		this->b *= c1.b;
		this->a *= c1.a;
		return *this;
	}

	template<typename T>
	color_t operator*=(const T scale){
		this->r *= scale;
		this->g *= scale;
		this->b *= scale;
		this->a *= scale;
		return *this;
	}

	template<typename T>
	color_t operator/=(const T scale){
		this->r /= scale;
		this->g /= scale;
		this->b /= scale;
		//this->a /= scale;
		return *this;
	}



private:
	template<typename T>
	void clamp_comp(T&comp, const T val){
		comp = (comp > val) ? val : comp;
		comp = (comp < 0) ? 0 : comp;
	}
};


//********************** STRINGIFY
//http://www.decompile.com/cpp/faq/file_and_line_error_string.htm
#ifndef STRINGIFY
#define STRINGIFY(x) TOSTRING(x)
#define TOSTRING(x) #x
#endif
//******************************************************************************

//********************** Command Line Parser
inline char* getCmdOption(char ** begin, char ** end, const std::string &option)
{
	//https://stackoverflow.com/a/868894/1608232
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
	{
		return *itr;
	}
	return 0;
}
inline bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
	//https://stackoverflow.com/a/868894/1608232
	return std::find(begin, end, option) != end;
}
//******************************************************************************

//********************** PRINT_ERROR
#ifndef PRINT_ERROR
#include <stdio.h>
#include <string>
inline void Err(std::string err_line, const char *file, int line) {
	//Display the err_line 	
	printf("Error::%s \n Error generated in %s at line %d\n", err_line.c_str(), file, line);

#ifdef _WIN32
	system("pause");
#else
	exit(EXIT_FAILURE);
#endif

}
#define PRINT_ERROR( err_line ) (Err( err_line, __FILE__, __LINE__ ))
#endif 
//******************************************************************************

template<typename T>
inline T Dist(T x1, T y1, T z1, T x2, T y2, T z2)
{
	T dx, dy, dz;
	dx = x1 - x2;
	dy = y1 - y2;
	dz = z1 - z2;
	dx *= dx;
	dy *= dy;
	dz *= dz;

	return dx + dy + dz;

}

template<typename T>
inline void cross_product(T xv1, T yv1, T zv1, T xv2, T yv2, T zv2,
	T cp[3])
{
	cp[0] = yv1*zv2 - zv1*yv2;
	cp[1] = zv1*xv2 - xv1*zv2;
	cp[2] = xv1*yv2 - yv1*xv2;
}

template<typename T>
inline void normalize_vector(T&vector_x, T&vector_y, T&vector_z)
{
	T nn = sqrt(vector_x*vector_x + vector_y*vector_y + vector_z*vector_z);
	vector_x /= nn; vector_y /= nn; vector_z /= nn;
}

template<typename T>
inline void normalize_vector(T vect[3])
{
	T nn = sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]);
	vect[0] /= nn; vect[1] /= nn; vect[2] /= nn;
}

template<typename T>
inline T Dot(T xv1, T yv1, T zv1, T xv2, T yv2, T zv2)
{
	return xv1*xv2 + yv1*yv2 + zv1*zv2;
}
template<typename T>
inline T Dot(T v1[3], T v2[3])
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

template<typename T>
inline bool ray_plane_intersect(T pp_x, T pp_y, T pp_z, //point on plane 
	T pv_x, T pv_y, T pv_z, //normal to the plane 
	T ux, T uy, T uz, //ray direction
	T ip1_x, T ip1_y, T ip1_z, //ray starting point 
	T&point_x, T&point_y, T&point_z// returned point on ray (if exists)
	){

	//plane line intersection. plane define by normal vector (pv_x,pv_y,pv_z) 
	//and point on it(pp_x,pp_y,pp_z) and ray starting at (ip1_x,ip1_y,ip1_z) 
	//and in direction of vector (ux,uy,uz) and returns true and point on the ray
	//(point_x,point_y,point_z) if the ray intersects the plane. Otherwise, 
	//returns false 

	T dot = Dot(ux, uy, uz, pv_x, pv_y, pv_z);

	if (abs(dot) <= 0.0){
		return false;
	}

	
	T s = Dot(pv_x, pv_y, pv_z, pp_x - ip1_x, pp_y - ip1_y, pp_z - ip1_z) / dot;

	
	//if (s<-1.0*10E-12 || s>1.0 + 10E-12){
	//	return false;
	//}
	point_x = ip1_x + s*ux;
	point_y = ip1_y + s*uy;
	point_z = ip1_z + s*uz;

	return true;

}

template <typename T>
T vector_mag(T vec[3]){
	T mag = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	return mag;
}

template<typename T>
inline void compute_avg_stddev(const T*arr, uint32_t size, double &avg,
	double&stddev){
	if (size == 1){
		avg = arr[0];
		stddev = 0;
		return;
	}
	avg = 0;
	//compute avg
	for (uint32_t i = 0; i < size; i++){
		avg += double(arr[i]);
	}
	avg /= size;

	//compute stddev
	double sum = 0;
	for (uint32_t i = 0; i < size; i++){
		float diff = double(arr[i]) - avg;
		sum += diff*diff;
	}
	stddev = sqrt(sum / double(size - 1));
	return;

}

#endif /*__COMMON__*/