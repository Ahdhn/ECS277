#ifndef __SCAT_DATA__
#define __SCAT_DATA__

#include <assert.h>

#include "common.h"

#include "KdTree.h"


template <class T, class T_d>
class ScatData2D
{
public:

	//Handle 2d scatter data points only 
	ScatData2D(T K) : m_kd_tree(KdTree()), m_K(K){
		m_point = new T_d[3];
		m_x0[0] = m_x0[1] = val_min = std::numeric_limits<T_d>::max();
		m_x1[0] = m_x1[1] = val_max = -val_min;	

	};

	void add_data(T_d x, T_d y, T_d f){
		scat_data2d_t<T_d> dd;
		dd.x = x;		
		dd.y = y;		
		dd.f = f;
		m_data.push_back(dd);
		set_max_min(m_data.size() - 1);
	}

	void fill_data(T data_id, scat_data2d_t<T_d> data){
		
		if (data_id >= m_data.size()){
			PRINT_ERROR("ScatData2D::fill_data() id is out of range!!!");
		}
		m_data[data_id].x = data.x;
		m_data[data_id].y = data.y;		
		m_data[data_id].f = data.f;
		set_max_min(data_id);
	}

	void fill_data(T data_id, T_d x, T_d y, T_d z, T_d f){
		if (data_id >= m_data.size()){
			PRINT_ERROR("ScatData2D::fill_data() id is out of range!!!");
		}
		m_data[data_id].x = x;
		m_data[data_id].y = y;		
		m_data[data_id].f = f;

		set_max_min(data_id);
	}

	void fill_data(T data_id, T_d f){
		if (data_id >= m_data.size()){
			PRINT_ERROR("ScatData2D::fill_data() id is out of range!!!");
		}		
		m_data[data_id].f = f;	
		set_max_min(data_id);
	}

	scat_data2d_t<T_d> get_data(const T data_id){		
		return m_data[data_id];
	}

	T_d get_data_coord(const T data_id, const T dim){
		if (dim == 0){
			return m_data[data_id].x;
		}
		else if (dim == 1){
			return m_data[data_id].y;
		}
		PRINT_ERROR("Error at ScatData2D::get_data_coord() invalid dim");		
	}
	
	T_d get_data_value(const T data_id){
		return m_data[data_id].f;
	}
		
	T get_num_data(){
		return m_data.size();
	}
				
	T_d get_max_val(){
		return val_max;
	}

	T_d get_min_val(){
		return val_min;
	}

	T_d get_min_x(){
		return m_x0[0];
	}

	T_d get_min_y(){
		return m_x0[1];
	}

	T_d get_max_x(){
		return m_x1[0];
	}

	T_d get_max_y(){
		return m_x1[1];
	}

	void export_data(std::string filename){
		//TODO
	};

	void build_kdtree(){
		m_kd_tree.BuildTree(m_data);
	}

	void get_knn(const T pnt_id, const T K, std::vector<int>& k_nearest){
		m_point[0] = m_data[pnt_id].x;
		m_point[1] = m_data[pnt_id].y;
		m_point[2] = 0;
		m_kd_tree.FindNNearest(m_point, K, k_nearest);
	}
	
	void normalize_values(){
		//normalize the input function to be (0,1)		
		T_d l = val_max - val_min;
		for (T i = 0; i < m_data.size(); i++){
			m_data[i].f -= val_min;			
			m_data[i].f /= l;
		}
		val_min = 0;
		val_max -= val_min;
		val_max /= l;
	}

	void run_multi_res();

	~ScatData2D(){};

private:
	
	inline void set_max_min(const T id);

	inline bool check_point(const T_d x, const T_d y);
		
	T_d val_min, val_max;

	std::vector<scat_data2d_t<T_d>>m_data;
	std::vector<T>m_data_current;
	std::vector<T>m_data_leftout;
	std::vector<T>m_interpol;
	std::vector<T>m_error;

	KdTree m_kd_tree;
	T m_K;
	
	T_d m_x0[2]; //lower most corner 
	T_d m_x1[2]; //upper most corner 

	std::vector<int> m_k_nearest;

	double *m_point;
	
};

template <class T, class T_d>
inline void ScatData2D<T, T_d>::set_max_min(const T id){

	m_x0[0] = std::min(m_x0[0], m_data[id].x);
	m_x0[1] = std::min(m_x0[1], m_data[id].y);	

	m_x1[0] = std::max(m_x1[0], m_data[id].x);
	m_x1[1] = std::max(m_x1[1], m_data[id].y);

	val_min = std::min(val_min, m_data[id].f);
	val_max = std::max(val_max, m_data[id].f);
	
}

template <class T, class T_d>
inline bool ScatData2D<T, T_d>::check_point(const T_d x, const T_d y){

	if (x - EPSILON > m_x1[0] || y - EPSILON > m_x1[1] || 
		x + EPSILON < m_x0[0] || y + EPSILON < m_x0[1] ){
		return false;
	}
	return true;
}

template <class T, class T_d>
void ScatData2D<T, T_d>::run_multi_res(){
	
	//prepare containers 
	m_data_current.clear();
	m_data_current.reserve(m_data.size());
	m_data_leftout.reserve(m_data.size());
	m_interpol.reserve(m_data.size());//stores the interpolation of the leftout
	m_error.reserve(m_data.size());//stores the error of the current data points
	
	//initial vorodel (to setup the boundary points)
	VoroDel<T, T_d> vorodel(m_K, m_data.size());
	vorodel.construct_del(this, m_data_current);
	
	//insert all boundary point as starting point 
	/*for (T i = 0; i < m_data.size(); ++i){
		if (vorodel.is_boundary(i)){
			m_data_current.push_back(i);
		}
		else{
			if (i != 4){
				m_data_leftout.push_back(i);
			}
		}
	}
	m_data_current.push_back(4);*/

	for (T i = 0; i < m_data.size(); ++i){
		if (i == 4){ continue; }
		m_data_current.push_back(i);		
		
	}
	m_data_leftout.push_back(4);

	assert(m_data_current.size() + m_data_leftout.size() == m_data.size());

	//Voronoi shading lambda function 
	auto voro_shading = [this](index_t p){
		data_t f = this->get_data_value(p);
		color_t color;
		COLOR_MAP cm = COLOR_MAP::RAINBOW;
		colormap(cm,
			(f - this->get_min_val()) /
			(this->get_max_val() - this->get_min_val())
			, color.r, color.g, color.b);
		color.a = 1.0;
		return color;

	};
	vorodel.plot("test.ps", this, false, true, voro_shading, true);

	for (;;){
		assert(m_data_current.size() + m_data_leftout.size() == m_data.size());
		//build new kd-tree on the current data
		m_kd_tree.~KdTree();
		m_kd_tree.BuildTree(m_data, m_data_current);

		//build new vorodel on the current data
		vorodel.update_k(std::min(T(m_data_current.size() - 2), m_K));
		vorodel.construct_del(this, m_data_current);
		vorodel.plot("test.ps", this, false, true, voro_shading, true);
		
		//compute interpolation of the left out 
		for (T i = 0; i < m_data_leftout.size(); ++i){
			vorodel.sibson_interpol(this, m_data_current.size(), m_data_leftout[i]);
		}


		//compute the error of each voronoi cell 
		

	}
}

#endif/*__SCAT_DATA__*/