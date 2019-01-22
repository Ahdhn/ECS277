#ifndef __COMMON__
#define __COMMON__

struct color_t {
	float r, g, b, a;
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
inline T Dist(T x1, double y1, T z1, T x2, T y2, T z2)
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
#endif /*__COMMON__*/