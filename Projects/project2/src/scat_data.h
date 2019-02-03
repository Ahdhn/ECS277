#ifndef __SCAT_DATA__
#define __SCAT_DATA__

#include "common.h"

#include "KdTree.h"

template <class T, class T_d>
class ScatData
{
public:
	
	ScatData();

	ScatData(bool is_analytical, int num_data);

	void fill_data(T data_id, scat_data_t<T_d> data){
		
		if (data_id >= m_num_data){
			PRINT_ERROR("ScatData::fill_data() id is out of range!!!");
		}
		m_data[data_id].x = data.x;
		m_data[data_id].y = data.y;
		m_data[data_id].z = data.z;
		m_data[data_id].f = data.f;

		set_max_min(data_id);
	}

	void fill_data(T data_id, T_d x, T_d y, T_d z, T_d f){
		if (data_id >= m_num_data){
			PRINT_ERROR("ScatData::fill_data() id is out of range!!!");
		}
		m_data[data_id].x = x;
		m_data[data_id].y = y;
		m_data[data_id].z = z;
		m_data[data_id].f = f;

		set_max_min(data_id);
	}

	scat_data_t<T_d> get_data(const T data_id){		
		return m_data[data_id];
	}

	void precomute(INTERPOL_METHOD method, uint32_t K = 5){
		if (!m_kd_built){
			m_kd_tree->BuildTree(m_data);
			m_kd_built = true;
		}

		if (method == INTERPOL_METHOD::S1){
			//nothing to do here
			m_s1_precomputed = true;
			return;
		}
		else if (method == INTERPOL_METHOD::S2){
			//nothing to do here
			m_s2_precomputed = true;
			return;
			
		}
		else if (method == INTERPOL_METHOD::S3){
			//precompute the gradient at each point using KNN data point 

		}
		else if (method == INTERPOL_METHOD::H){
			//precompute the coefficient C in Hardy's 

		}
		else {
			PRINT_ERROR("ScatData::precompute() Unknown interpolation method");
		}
	}
		
	void get_data(const T data_id, T_d&x, T_d&y, T_d&z, T_d&f){
		if (data_id >= m_num_data){
			PRINT_ERROR("ScatData::get_data() data_id out of range");
		}
		x = m_data[data_id]->x;
		y = m_data[data_id]->y;
		z = m_data[data_id]->z;
		f = m_data[data_id]->f;
	}
		
	T get_num_data(){
		return m_num_data;
	}

	T_d interpolate(T_d x, T_d y, T_d z, INTERPOL_METHOD method);
		
	void export_data(std::string filename);
	~ScatData(){};

private:
	inline T_d shepard(const T_d x, const T_d y, const T_d z,
		const T shepard_num);

	inline T_d hardy(T_d&x, T_d&y, T_d&z);

	inline void set_max_min(const T id);
	inline bool check_point(const T_d x, const T_d y, const T_d z);

	bool m_is_analytical;
	T m_num_data;	
	scat_data_t<T_d> * m_data;
	
	bool m_s1_precomputed, m_s2_precomputed,m_s3_precomputed, m_h_precomputed;

	KdTree *m_kd_tree;
	double *m_point;
	bool m_kd_built;
	std::vector<int>m_knn_container;//used to contain the returned indeices of knn
	T_d m_x0[3]; //lower most corner 
	T_d m_x1[3]; //upper most corner 
};

template <class T, class T_d>
ScatData<T, T_d>::ScatData(bool is_analytical, int num_data) :
m_is_analytical(is_analytical), m_num_data(num_data){

	m_data = (scat_data_t<T_d>*)malloc(m_num_data*sizeof(scat_data_t<T_d>));

	m_s1_precomputed = m_s2_precomputed= m_s3_precomputed = 
		m_h_precomputed = false;

	static KdTree kd(3, m_num_data);
	m_kd_tree = &kd;	

	m_point = new double[3];

	m_kd_built = false;

	m_x0[0] = m_x0[1] = m_x0[1] = std::numeric_limits<T_d>::max();
	m_x1[0] = m_x1[1] = m_x1[1] = -std::numeric_limits<T_d>::max();
}

template <class T, class T_d>
inline void ScatData<T, T_d>::set_max_min(const T id){
	m_x0[0] = std::min(m_x0[0], m_data[id].x);
	m_x0[1] = std::min(m_x0[1], m_data[id].y);
	m_x0[2] = std::min(m_x0[2], m_data[id].z);

	m_x1[0] = std::max(m_x1[0], m_data[id].x);
	m_x1[1] = std::max(m_x1[1], m_data[id].y);
	m_x1[2] = std::max(m_x1[2], m_data[id].z);
}

template <class T, class T_d>
inline bool ScatData<T, T_d>::check_point(const T_d x, const T_d y, const T_d z){

	if (x - EPSILON > m_x1[0] || y - EPSILON > m_x1[1] || z - EPSILON> m_x1[2] ||
		x + EPSILON < m_x0[0] || y + EPSILON < m_x0[1] || z + EPSILON < m_x0[2]){
		return false;
	}
	return true;
}

template<class T, class T_d>
inline T_d ScatData<T, T_d>::interpolate(T_d x, T_d y, T_d z, 
	INTERPOL_METHOD method){

	//if (!check_point(x, y, z)){
	//	PRINT_ERROR("ScatData::interpolate() input point outside the bounding box ("+
	//		std::to_string(x)+ ", "+ std::to_string(y)+ ", "+ std::to_string(z) + ")");
	//}

	if (method == INTERPOL_METHOD::S1){
		return shepard(x, y, z, 1);
	}
	else if (method == INTERPOL_METHOD::S2){
		return shepard(x, y, z, 2);
	}
	else if (method == INTERPOL_METHOD::S3){
		return shepard(x, y, z, 3);
	}
	else if (method == INTERPOL_METHOD::H){

	}
	else {
		PRINT_ERROR("ScatData::interpolate() Unknown interpolation method");
	}

}

template<class T, class T_d>
inline T_d ScatData<T, T_d>::shepard(const T_d x, const T_d y, const T_d z, 
	const T shepard_num){
	if (!m_s1_precomputed){
		PRINT_ERROR("ScatData::shepard_one() precompute function should be called first");
	}

	m_point[0] = x;
	m_point[1] = y;
	m_point[2] = z;
	int nearest_id = m_kd_tree->FindNearest(m_point);
	if (Dist(x, y, z, m_data[nearest_id].x, m_data[nearest_id].y,
		m_data[nearest_id].z) < EPSILON){
		return m_data[nearest_id].f;
	}

	T_d ff_num = 0;
	T_d ff_den = 0;
	for (uint32_t i = 0; i < m_num_data; i++){
		
		T_d dist = Dist(x, y, z, m_data[i].x, m_data[i].y, m_data[i].z);

		if (shepard_num == 1){
			dist = sqrt(dist);
		}

		dist = 1.0 / dist;
		T_d fi;

		if (shepard_num == 1 || shepard_num == 2){
			fi = m_data[i].f;
		}		
		else if (shepard_num == 3){
			
		}
		ff_num += dist*fi;
		ff_den += dist;
	}

	return ff_num / ff_den;
}


template<class T, class T_d>
inline T_d ScatData<T, T_d>::hardy(T_d&x, T_d&y, T_d&z){
	if (!m_h_precomputed){
		PRINT_ERROR("ScatData::hardy() precompute function should be called first");
	}
}

template <class T, class T_d>
void ScatData<T, T_d>::export_data(std::string filename){
	//export data as points in csv format
	std::fstream file(filename, std::ios::out);
	file.precision(12);
	file << "x coord,y coord,z coord,intensity" << std::endl;

	for (uint32_t i = 0; i < m_num_data; i++){
		file << m_data[i].x << ", " << m_data[i].y << ", " << m_data[i].z
			<< ", " << m_data[i].f << std::endl;
	}

	file.close();
}
#endif/*__SCAT_DATA__*/