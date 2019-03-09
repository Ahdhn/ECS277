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
		if (K == 1){
			int nearest = m_kd_tree.FindNearest(m_point);
			k_nearest.clear();
			k_nearest.push_back(nearest);
		}
		else{
			m_kd_tree.FindNNearest(m_point, K, k_nearest);
		}
		
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
	std::vector<T_d>m_interpol;
	std::vector<T_d>m_error;
	std::vector<T_d>m_error_max;
	std::vector<T>m_tile_count;

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
	m_interpol.resize(m_data.size());//stores the interpolation of the leftout
	m_error.resize(m_data.size());//stores the error of the current data points(RMS)
	m_error_max.resize(m_data.size());//stores the error of the current data points(max)
	m_tile_count.resize(m_data.size());//number of left out in this tile 

	//initial vorodel (to setup the boundary points)
	VoroDel<T, T_d> vorodel(m_K, m_data.size());
	vorodel.construct_del(this, m_data_current);
	
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

	
	vorodel.plot("voro.ps", this, false, true, voro_shading, false);

	//insert all boundary point as starting point 
	
	for (T i = 0; i < m_data.size(); ++i){
		if (vorodel.is_boundary(i)){
			m_data_current.push_back(i);
		}
		else{
			if (i != 8){
				m_data_leftout.push_back(i);
			}
		}
	}
	m_data_current.push_back(8);

	int step = 0;
	while (m_data_current.size() <= m_data.size()){		
		std::cout << " step = " << step << " current = " <<m_data_current.size()
			<< " leftout = " << m_data_leftout.size() << std::endl;
		assert(m_data_current.size() + m_data_leftout.size() == m_data.size());
		memset(m_interpol.data(), 0, m_interpol.size()*sizeof(T_d));
		memset(m_error.data(), 0, m_error.size()*sizeof(T_d));
		memset(m_error_max.data(), 0, m_error_max.size()*sizeof(T_d));
		memset(m_tile_count.data(), 0, m_tile_count.size()*sizeof(T));


		//build new kd-tree on the current data
		m_kd_tree.~KdTree();
		m_kd_tree.BuildTree(m_data, m_data_current);

		//build new vorodel on the current data
		vorodel.update_k(std::min(T(m_data_current.size() - 2), m_K));
		vorodel.construct_del(this, m_data_current);

		if (step % 100 == 0){
			std::string name_v = "voro" + std::to_string(step) + ".ps";
			//std::string name_v = "voro.ps";
			vorodel.plot(name_v, this, false, true, voro_shading, false);
		}
		
		if (m_data_leftout.size() == 0){
			std::string name_v = "voro" + std::to_string(step) + ".ps";			
			vorodel.plot(name_v, this, false, true, voro_shading, false);
			break;
		}

		//compute interpolation of the left out 
		for (T i = 0; i < m_data_leftout.size(); ++i){
			
			m_interpol[m_data_leftout[i]] = vorodel.sibson_interpol(this,
				m_data_current.size(), m_data_leftout[i]);			
		}


		//compute the error of each voronoi cell 
		//each left out point will contribute to the error of the cell it lies 
		//inside i.e., its closest current data so we can use the kdtree we have 
		//for this 		
 		for (T i = 0; i < m_data_leftout.size(); ++i){
			get_knn(m_data_leftout[i], 1, m_k_nearest);

			T tile = m_k_nearest[0];

			T_d f_tile = get_data_value(tile);
						
			T_d diff = f_tile - m_interpol[m_data_leftout[i]];

			m_error[tile] += diff*diff;

			//m_error_max[tile] = std::max(abs(diff), m_error_max[tile]);

			m_tile_count[tile]++;
		}


		T_d max_error(-100000), min_error(100000);
		T max_error_tile_id;
		for (T i = 0; i < m_data_current.size();i++){
			T tile = m_data_current[i];
			if (m_tile_count[tile] > 0){
				m_error[tile] /= m_tile_count[tile];

				if (m_error[tile] > max_error){
					max_error = m_error[tile];
					max_error_tile_id = tile;
				}

				if (m_error[tile] < min_error){
					min_error = m_error[tile];
				}
			}
		}		

		//shading lambda function 
		auto voro_shading_error = [this, &max_error, &min_error](index_t p){
			data_t e = this->m_error[p];			
			color_t color;
			COLOR_MAP cm = COLOR_MAP::MAGMA;
			colormap(cm,
				(e - min_error) /
				(max_error - min_error)
				, color.r, color.g, color.b);
			color.a = 1.0;
			return color;
		};

		if (step % 100 == 0 || m_data_leftout.size() == 1){
			std::string name_e = "voro_err" + std::to_string(step) + ".ps";
			//std::string name_e = "voro_err.ps";
			vorodel.plot(name_e, this, false, true, voro_shading_error, false);
		}

		//insert one of the leftout points that lie inside the tile_max_error 		
		for (T i = 0; i < m_data_leftout.size(); i++){
			get_knn(m_data_leftout[i], 1, m_k_nearest);
			if (m_k_nearest[0] == max_error_tile_id){
				//insert 
				std::cout << "Max Error = " << max_error << " at tile " <<
					max_error_tile_id << "at ("
					<< get_data_coord(max_error_tile_id, 0) 
					<< ", "
					<< get_data_coord(max_error_tile_id, 1) 
					<< ") --> ";
					std::cout << m_data_leftout[i] << " inserted at ("
					<< get_data_coord(m_data_leftout[i], 0) << ", "
					<< get_data_coord(m_data_leftout[i], 1) << ") " 
					<< std::endl << std::endl;

				m_data_current.push_back(m_data_leftout[i]);
				m_data_leftout.erase(m_data_leftout.begin() + i);				
				break;
			}
		}
		
		step++;
	}
}

#endif/*__SCAT_DATA__*/