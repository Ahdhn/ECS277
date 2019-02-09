#ifndef __SCAT_DATA__
#define __SCAT_DATA__

#include "common.h"

#include "KdTree.h"

#define EIGEN_DONT_ALIGN
#include <Eigen/Dense>

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

	void precompute(INTERPOL_METHOD method, uint32_t K = 5, double R = 0.1){
		if (!m_kd_built){
			m_kd_tree->BuildTree(m_data);
			m_kd_built = true;
		}						
		if (method == INTERPOL_METHOD::H_G_MQ || 
			method == INTERPOL_METHOD::H_G_RE){
			m_h_precomputed = true;
			//precompute the coefficient C in Hardy's 

			for (T i = 0; i < m_num_data; i++){
				for (T j = 0; j < m_num_data; j++){
					
					m_M_hardy(i, j) = sqrt(R*R + Dist(m_data[i].x, m_data[i].y,
						m_data[i].z, m_data[j].x, m_data[j].y, m_data[j].z));
					if (method == INTERPOL_METHOD::H_G_RE){
						m_M_hardy(i, j) = 1.0 / m_M_hardy(i, j);
					}
				}
				m_F_hardy(i) = m_data[i].f;
			}
			//solve the system 
			m_C_hardy = m_M_hardy.colPivHouseholderQr().solve(m_F_hardy);
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

	T_d interpolate(T_d x, T_d y, T_d z, INTERPOL_METHOD method, const T K = 5,
		const T_d R = 0.1);
		
	void export_data(std::string filename);
	~ScatData(){};

private:
	inline T_d shepard_global(const T_d x, const T_d y, const T_d z,
		const T shepard_num);
	
	inline T_d shepard_local(const T_d x, const T_d y, const T_d z,
		const T shepard_num, const T K = 5);
	
	inline T_d hardy_global(T_d&x, T_d&y, T_d&z, bool is_reciprocal = false,
		const T_d R = 0.1);

	inline T_d hardy_local(T_d&x, T_d&y, T_d&z, bool is_reciprocal = false,
		const T K = 5, const T_d R = 0.1);


	inline void set_max_min(const T id);
	inline bool check_point(const T_d x, const T_d y, const T_d z);

	bool m_is_analytical;
	T m_num_data;	
	scat_data_t<T_d> * m_data;
	
	bool m_h_precomputed;

	KdTree *m_kd_tree;
	double *m_point;
	bool m_kd_built;
	std::vector<int>m_knn_container;//used to contain the returned indeices of knn
	T_d m_x0[3]; //lower most corner 
	T_d m_x1[3]; //upper most corner 

	std::vector<int> m_k_nearest;

	//Hardy matrices 
	Eigen::MatrixXd m_M_hardy;
	Eigen::VectorXd m_C_hardy;
	Eigen::VectorXd m_F_hardy;
};

template <class T, class T_d>
ScatData<T, T_d>::ScatData(bool is_analytical, int num_data) :
m_is_analytical(is_analytical), m_num_data(num_data){

	m_data = (scat_data_t<T_d>*)malloc(m_num_data*sizeof(scat_data_t<T_d>));

	m_h_precomputed = false;

	static KdTree kd(3, m_num_data);
	m_kd_tree = &kd;	

	m_point = new double[3];

	m_kd_built = false;

	m_x0[0] = m_x0[1] = m_x0[1] = std::numeric_limits<T_d>::max();
	m_x1[0] = m_x1[1] = m_x1[1] = -std::numeric_limits<T_d>::max();

	m_M_hardy.resize(m_num_data, m_num_data);
	m_C_hardy.resize(m_num_data);
	m_F_hardy.resize(m_num_data);

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
	INTERPOL_METHOD method, const T K /* = 5*/, const T_d R /*= 0.1*/ ){

	//if (!check_point(x, y, z)){
	//	PRINT_ERROR("ScatData::interpolate() input point outside the bounding box ("+
	//		std::to_string(x)+ ", "+ std::to_string(y)+ ", "+ std::to_string(z) + ")");
	//}

	if (method == INTERPOL_METHOD::S1_G){
		return shepard_global(x, y, z, 1);
	}else if (method == INTERPOL_METHOD::S1_L){
		return shepard_local(x, y, z, 1, K);
	}
	else if (method == INTERPOL_METHOD::S2_G){
		return shepard_global(x, y, z, 2);
	}
	else if (method == INTERPOL_METHOD::S2_L){
		return shepard_local(x, y, z, 2, K);
	}
	else if (method == INTERPOL_METHOD::H_G_MQ){
		return hardy_global(x, y, z, false);
	}
	else if (method == INTERPOL_METHOD::H_L_MQ){
		return hardy_local(x, y, z, false, K);
	}
	else if (method == INTERPOL_METHOD::H_G_RE){
		return hardy_global(x, y, z, true);
	}
	else if (method == INTERPOL_METHOD::H_L_RE){
		return hardy_local(x, y, z, true, K);
	}
	else {
		PRINT_ERROR("ScatData::interpolate() Unknown interpolation method");
	}

}

template<class T, class T_d>
inline T_d ScatData<T, T_d>::shepard_global(const T_d x, const T_d y, const T_d z,
	const T shepard_num){
	
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
inline T_d ScatData<T, T_d>::shepard_local(const T_d x, const T_d y, const T_d z,
	const T shepard_num, const T K /* = 5*/){

	m_point[0] = x;
	m_point[1] = y;
	m_point[2] = z;
	m_kd_tree->FindNNearest(m_point, K, m_k_nearest);
	if (m_k_nearest.size() != K){
		PRINT_ERROR("ScatData::shepard_local() can not find the K nearest neighbours!!");
	}
	if (Dist(x, y, z, m_data[m_k_nearest[0]].x, m_data[m_k_nearest[0]].y,
		m_data[m_k_nearest[0]].z) < EPSILON){
		return m_data[m_k_nearest[0]].f;
	}

	T_d ff_num = 0;
	T_d ff_den = 0;
	for (uint32_t p = 0; p < m_k_nearest.size(); p++){
		uint32_t i = m_k_nearest[p];

		T_d dist = Dist(x, y, z, m_data[i].x, m_data[i].y, m_data[i].z);

		if (shepard_num == 1){
			dist = sqrt(dist);
		}
		dist = 1.0 / dist;

		T_d fi = m_data[i].f;
				
		ff_num += dist*fi;
		ff_den += dist;
	}

	return ff_num / ff_den;

}


template<class T, class T_d>
inline T_d ScatData<T, T_d>::hardy_global(T_d&x, T_d&y, T_d&z, 
	bool is_reciprocal /*= false*/, const T_d R /*= 0.1*/){

	if (!m_h_precomputed){
		PRINT_ERROR("ScatData::hardy() precompute function should be called first");
	}

	T_d val = 0;
	for (T i = 0; i < m_num_data; i++){
		T_d tt = R*R + Dist(x, y, z, m_data[i].x, m_data[i].y, m_data[i].z);
		tt = sqrt(tt);
		if (is_reciprocal){
			tt = 1.0 / tt;
		}		
		val += m_C_hardy(i)*tt;
	}
	return val;
}

template<class T, class T_d>
inline T_d ScatData<T, T_d>::hardy_local(T_d&x, T_d&y, T_d&z,
	bool is_reciprocal /*= false*/, const T K /* = 5*/, const T_d R /*= 0.1*/){
	




	return 0;
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