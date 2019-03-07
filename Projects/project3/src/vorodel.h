#ifndef __VORO__
#define __VORO__

#include <assert.h>

#include "scat_data_2d.h"

template<class T, class T_d>
class VoroDel
{
public:

	VoroDel(const T K, const T num_points_all): m_K(K) {
		m_angles.resize(K + 10);
		m_pseudo_del.reserve(K + 10);

		for (T i = 0; i < num_points_all; i++){
			std::vector<T> dd;
			m_del.push_back(dd);
		}

		m_boundary.resize(num_points_all);
		for (T p = 0; p < num_points_all; p++){
			m_boundary[p] = false;
		}

	};
		
	void construct_del(ScatData2D<T, T_d>*scat_data, 
		const std::vector<T>&current_data);

	void reset(){
		m_del.clear();
	}
	
	template<typename VoroCol>
	void plot(std::string filename, ScatData2D<T, T_d>*scat_data,
		const bool delaunay, const bool voronoi, VoroCol voronoi_cell_shading,
		const bool skip_boundary);
		
	T_d sibson_interpol(ScatData2D<T, T_d>*scat_data, const T num_points_current,
		const T new_pnt);

	bool is_boundary(T id){
		return m_boundary[id];
	}

	void update_k(T new_k){
		m_K = new_k;
	}

	~VoroDel(){}

private:

	bool is_delaunay(const T p, const T q, ScatData2D<T, T_d>*scat_data, T& w);
	void sort_delaunay(const T p, std::vector<T>&p_del, 
		ScatData2D<T, T_d>*scat_data);
	inline void sort(T left, T right, std::vector<T>&p_del);
	inline void partition(T &i, T &j, std::vector<T>&p_del);	
	void pseudo_insert(ScatData2D<T, T_d>*scat_data, const T num_points_current,
		const T new_pnt);
	
		
	T find_common_vertex(ScatData2D<T, T_d>*scat_data, const T p,
		const T m, const T q, T_d&xx, T_d&yy);

	

	inline void plot_single_voronoi(ScatData2D<T, T_d>*scat_data, const T p,
		const bool with_del);

	std::vector<T> m_pseudo_del;

	std::vector<int>m_knn_container;
	std::vector<int>m_knn_search;
	std::vector<std::vector<T>> m_del;
	std::vector<bool> m_boundary;
	std::vector<T_d> m_angles;
	T m_K;
};


template<class T, class T_d>
void VoroDel<T, T_d>::construct_del(ScatData2D<T, T_d>*scat_data, 
	const std::vector<T>&current_data){
	//construct delauny from the data. current_data holds the ids of the 
	//points to build the delaunay on 
	//if current_data is empty vector, we will build delaunay for everybody 
	
	T num_points_current = current_data.size();
	T num_points_all = m_del.size();
	
	if (num_points_current == 0){
		num_points_current = num_points_all;
	}


	for (T p = 0; p < num_points_all; p++){
		m_boundary[p] = false;
	}
		
	for (T p = 0; p < num_points_all; p++){
		m_del[p].clear();
	}


	for (T p = 0; p < num_points_current; p++){
		
		//for each point get it K closest point 				
		T p_id = (current_data.size() == 0) ? p : current_data[p];
		
		scat_data->get_knn(p_id, m_K, m_knn_container);
		assert(m_knn_container.size() == m_K);
				
		scat_data->get_knn(p_id, 
			((2 * m_K) <  num_points_current -1 ? (2 * m_K) : num_points_current - 1),
			m_knn_search);
		
		
		if (current_data.size() != 0){
			for (T i = 0; i < m_knn_search.size(); i++){
				uint32_t id = get_index(int(m_knn_search[i]), current_data);
				if (id == std::numeric_limits<uint32_t>::max()){
					PRINT_ERROR("Can find the knn search in current data");
				}
			}
		}
		for (T i = 0; i < m_knn_container.size(); ++i){
			T q = m_knn_container[i];

			if (q == p_id){ continue; }

			//fist check if there are connected 
			
			uint32_t id = get_index(q, m_del[p_id]);
			if (id == std::numeric_limits<uint32_t>::max()){				

				T w(0);				
				if (is_delaunay(p_id, q, scat_data, w)){
					//if is boundary, it stays boundary 
															
					m_del[p_id].push_back(q);
					m_del[q].push_back(p_id);

					uint32_t iwp = get_index(w, m_del[p_id]);
					if (iwp == std::numeric_limits<uint32_t>::max()){
						m_del[w].push_back(p_id);
						m_del[p_id].push_back(w);
					}

					uint32_t iwq = get_index(w, m_del[q]);
					if (iwq == std::numeric_limits<uint32_t>::max()){
						m_del[q].push_back(w);
						m_del[w].push_back(q);
					}					
				}
			}			
		}
	}

	//sort delaunay 
	for (T p = 0; p < num_points_current; p++){
		//T p_id = current_data[p];
		T p_id = (current_data.size() == 0) ? p : current_data[p];
		if (m_del[p_id].size() == 0) { continue; }
		sort_delaunay(p_id, m_del[p_id], scat_data);
	}

	//set m_boundary

	T_d max_x = 1.1 * scat_data->get_max_x();
	T_d max_y = 1.1 * scat_data->get_max_y();

	T_d min_x = scat_data->get_min_x() - 0.1;
	T_d min_y = scat_data->get_min_y() - 0.1;	

	for (T p = 0; p < num_points_current; p++){
		//T p_id = current_data[p];
		T p_id = (current_data.size() == 0) ? p : current_data[p];
		if (m_del[p_id].size() == 0){ continue; }

		for (T i = 0; i < m_del[p_id].size(); i++){
			T q = m_del[p_id][i];
			T m = (i == m_del[p_id].size() - 1) ? m_del[p_id][0] :
				m_del[p_id][i + 1];
			if (p_id < q && p_id < m){
				T_d x_c(0), y_c(0);
				T_d r_2 = TriCircumcenter2d(scat_data->get_data_coord(p_id, 0),
					scat_data->get_data_coord(p_id, 1),
					scat_data->get_data_coord(q, 0),
					scat_data->get_data_coord(q, 1),
					scat_data->get_data_coord(m, 0),
					scat_data->get_data_coord(m, 1),
					x_c, y_c);

				if (x_c > max_x || x_c < min_x ||
					y_c > max_y || y_c < min_y){
					m_boundary[p_id] = true;
					m_boundary[q] = true;
					m_boundary[m] = true;
				}
			}
		}
	}


}

template<class T, class T_d>
bool VoroDel<T, T_d>::is_delaunay(const T p, const T q,
	ScatData2D<T, T_d>*scat_data, T& w){

	//check if the p-q form a delaunay edge by checking if p-q is a part of 
	//any delaunay triangle. If true, return w as the third point of this 
	//triangle

	//expect m_knn_constainer to be filled with the K closest neighbour to p 
		
	//construct the circumcircle of p,q,m 

	for (T i = 0; i < m_knn_container.size(); ++i){
		T m = m_knn_container[i];

		if (m == p || m == q){ continue; }

		T_d x_cir(0), y_cir(0);

		T_d r_2 = TriCircumcenter2d(
			scat_data->get_data_coord(p, 0), scat_data->get_data_coord(p, 1),
			scat_data->get_data_coord(q, 0), scat_data->get_data_coord(q, 1),
			scat_data->get_data_coord(m, 0), scat_data->get_data_coord(m, 1),
			x_cir, y_cir);

		bool is_del = true;

		for (T j = 0; j < m_knn_search.size(); ++j){
			T s = m_knn_search[j];

			if (s == p || s == q || s == m){ continue; }

			T_d dist = Dist<T_d>(scat_data->get_data_coord(s, 0),
				scat_data->get_data_coord(s, 1), T_d(0), x_cir, y_cir, T_d(0));

			if (dist < r_2){				
				is_del = false;
				break;				
			}
		}

		if (is_del){			
			w = m;
			return true;
		}
	}
	return false;
}

template<class T, class T_d>
void VoroDel<T, T_d>::sort_delaunay(const T p, std::vector<T>&p_del,
	ScatData2D<T, T_d>*scat_data){
	//for each point, sort its delaunay based of the angle 
		
	T num_del = p_del.size();
	T_d xp = scat_data->get_data_coord(p, 0);
	T_d yp = scat_data->get_data_coord(p, 1);
	for (T i = 0; i < num_del; i++){
		
		T q = p_del[i];

		T_d dx = scat_data->get_data_coord(q, 0) - xp;
		T_d dy = scat_data->get_data_coord(q, 1) - yp;

		T_d theta = atan2(dy, dx);
		if (theta < 0){ theta += TwoPI; }

		m_angles[i] = theta;
	}

	sort(0, num_del - 1, p_del);

}

template<class T, class T_d>
inline void VoroDel<T, T_d>::sort(T left, T right, std::vector<T>&p_del){
	T i, j;
	i = left;
	j = right;
	partition(i, j, p_del);
	if (left < j){
		sort(left, j, p_del);
	}
	if (i < right){
		sort(i, right, p_del);
	}
}
template<class T, class T_d>
inline void VoroDel<T, T_d>::partition(T &i, T &j, std::vector<T>&p_del){

	T pivot = (i + j) / 2;
	T tmp;
	T_d tmp1;
	T_d pivot_value = m_angles[pivot];

	while (i <= j) {
		while (pivot_value < m_angles[j]){ j--; }
		while (pivot_value > m_angles[i]){ i++; }
		if (i <= j) {
			std::swap(p_del[i], p_del[j]);
			std::swap(m_angles[i], m_angles[j]);
			j--;
			i++;
		}
	}
}


template<class T, class T_d>
template<typename VoroCol>
void VoroDel<T, T_d>::plot(std::string filename, ScatData2D<T, T_d>*scat_data,
	const bool delaunay, const bool voronoi, VoroCol voronoi_cell_shading,
	const bool skip_boundary){

	if (scat_data->get_num_data() != m_del.size()){
		PRINT_ERROR("VoroDel::plot_del() the number of scatter data does not match the number of delaunay points");
	}


	filename = STRINGIFY(OUTPUT_DIR) + filename;

	std::fstream file(filename, std::ios::out);

	file.precision(30);

	file << "%!PS-Adobe-3.0" << std::endl;
	file << "72 72 scale     % one unit = one inch" << std::endl;

	double scale_x, scale_y, scale, shift_x, shift_y;

	double len = 1.2;

	scale_x = 6.5 / len;
	scale_y = 9.0 / len;

	if (scale_x < scale_y){
		scale = scale_x;
		shift_x = 0.9;
		shift_y = (11.0 - (len*scale)) / 2.0;
	}
	else {
		scale = scale_y;
		shift_x = (8.3 - (len*scale)) / 2.0;
		shift_y = 1.35;
	}

	file << shift_x << " " << shift_y << " translate" << std::endl;

#pragma region shapes
	file << "/quad      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/redquad      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << " 0.01 setlinewidth" << std::endl;
	file << "1 0 0 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/quad_bold      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0 0 setrgbcolor" << std::endl;
	//	file << "fill" << std::endl;
	file << " 0.01 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/fquad_bold      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0 0 setrgbcolor" << std::endl;
	file << "fill" << std::endl;
	file << " 0.01 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/blue_tri      % stack: x0 y0 x1 y1 x2 y2 " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;	
	file << " closepath" << std::endl;
	file << " 0 0 1 setrgbcolor" << std::endl;
	//file << "fill" << std::endl;
	file << " 0.002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;


	file << "/point      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.01 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/red_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/black_point      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_point      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/yellow_point      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 1 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/disk      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 1 setrgbcolor" << std::endl;
	//	file << " fill" << std::endl;
	file << " 0.004 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/green_disk      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 1 0 setrgbcolor" << std::endl;
	file << " 0.003 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/red_line      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 1 0 0 setrgbcolor" << std::endl;
	file << " 0.004 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/blue_line      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0 0 1 setrgbcolor" << std::endl;
	file << " 0.004 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;


	file << "/orange_disk      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	file << " 0.008 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/white_disk      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 1 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.008 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/pink_disk      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << "  1 0.7529411764705882 0.796078431372549 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.008 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/black_disk      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	file << " 0.008 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/txt      % stack: x y " << std::endl;
	file << "{newpath" << std::endl;
	file << "/Times-Roman findfont" << std::endl;
	file << "0.09 scalefont" << std::endl;
	file << "0 0 0 setrgbcolor" << std::endl;
	file << "setfont" << std::endl;
	file << " moveto" << std::endl;
	file << "(" << std::endl;
	file << ") " << "show" << std::endl;
	file << "} def" << std::endl;
#pragma endregion

	if (delaunay){
		for (T p = 0; p < m_del.size(); p++){

			if (m_del[p].size() == 0){ continue; }

			for (T i = 0; i < m_del[p].size(); i++){
				T q = m_del[p][i];

				if (p < 4){
					//the first four points are the corner points and we draw 
					//them as lines (not triangles) to avoid looping around ourselves 
					//and end up with triangles that are not part of the mesh
					continue;
					file << scat_data->get_data_coord(p, 0)*scale << " "
						<< scat_data->get_data_coord(p, 1)*scale << " "
						<< scat_data->get_data_coord(q, 0)*scale << " "
						<< scat_data->get_data_coord(q, 1)*scale << " blue_line"
						<< std::endl;
					continue;
				}

				T m = (i == m_del[p].size() - 1) ? m_del[p][0] : m_del[p][i + 1];
				if (p < m && p < q){
					//this condition ensure drawing each triangle just once 
					file << scat_data->get_data_coord(p, 0)*scale << " "
						<< scat_data->get_data_coord(p, 1)*scale << " "
						<< scat_data->get_data_coord(q, 0)*scale << " "
						<< scat_data->get_data_coord(q, 1)*scale << " "
						<< scat_data->get_data_coord(m, 0)*scale << " "
						<< scat_data->get_data_coord(m, 1)*scale << " blue_tri"
						<< std::endl;
				}
			}
		}
	}

	if (voronoi){
	
		//bool skip_boundary = true;
		for (T p = 0; p < m_del.size(); p++){	

			if (skip_boundary && (p < 4 || m_boundary[p])){
				continue;
			}

			if (m_del[p].size() == 0){ continue; }

			color_t shade = voronoi_cell_shading(p);

			file << " newpath " << std::endl;
			file << " 0.01 setlinewidth" << std::endl;
			file << " " << shade.r << " " << shade.g << " " << shade.b 
				<< " setrgbcolor" << std::endl;
			//file << " 1.0 0.0 0.0 setrgbcolor" << std::endl;


			for (T i = 0; i < m_del[p].size(); i++){
				T q = m_del[p][i];
				T m = (i == m_del[p].size() - 1) ? m_del[p][0] : m_del[p][i + 1];
				
				T_d x_c(0), y_c(0);
				T_d r_2 = TriCircumcenter2d(scat_data->get_data_coord(p, 0),
					scat_data->get_data_coord(p, 1),
					scat_data->get_data_coord(q, 0),
					scat_data->get_data_coord(q, 1),
					scat_data->get_data_coord(m, 0),
					scat_data->get_data_coord(m, 1),
					x_c, y_c);

				file << x_c*scale << " " << y_c*scale;
				if (i == 0){
					file << " moveto " << std::endl;
				}
				else{
					file << " lineto " << std::endl;
				}
			}

			//file << "fill" << std::endl;
			file << " closepath" << std::endl;
			file << "stroke" << std::endl;
			//system("pause");
		}
		
		
	}

	for (T p = 0; p < m_del.size(); p++){
		T_d r_point = 0.001;
		
		if (m_boundary[p]){
			T_d r_disk = 0.003;
			file << scat_data->get_data_coord(p, 0)*scale << " "
				<< scat_data->get_data_coord(p, 1)*scale << " " << r_disk*scale
				<< " yellow_point" << std::endl;
		}
			

		file << scat_data->get_data_coord(p, 0)*scale << " "
			<< scat_data->get_data_coord(p, 1)*scale << " " << r_point*scale
			<< " point" << std::endl;
	}
}

template<class T, class T_d>
void VoroDel<T, T_d>::pseudo_insert(ScatData2D<T, T_d>*scat_data,
	const T num_points_current, const T new_pnt){

	//psuedo insert the new point i.e., build its delaunay 
	scat_data->get_knn(new_pnt,
		((2 * m_K) <  num_points_current - 1 ? (2 * m_K) : num_points_current - 1),
		m_knn_search);

	m_pseudo_del.clear();

	for (T i = 0; i < m_knn_container.size(); ++i){
		T q = m_knn_container[i];

		uint32_t id = get_index(q, m_pseudo_del);
		if (id == std::numeric_limits<uint32_t>::max()){

			T w(0);


			if (is_delaunay(new_pnt, q, scat_data, w)){

				uint32_t q_id = get_index(q, m_pseudo_del);
				if (q_id == std::numeric_limits<uint32_t>::max()){
					m_pseudo_del.push_back(q);
				}

				uint32_t w_id = get_index(w, m_pseudo_del);
				if (w_id == std::numeric_limits<uint32_t>::max()){
					m_pseudo_del.push_back(w);
				}
			}
		}
	}

	sort_delaunay(new_pnt, m_pseudo_del, scat_data);
}

template<class T, class T_d>
T_d VoroDel<T, T_d>::sibson_interpol(ScatData2D<T, T_d>*scat_data, 
	const T num_points_current, const T new_pnt){
	//compute sibson interpolation for a new point 
	//that should be called after constructing the delaunay 

	//find the closest point 
	scat_data->get_knn(new_pnt, m_K, m_knn_container);
	assert(m_knn_container.size() == m_K);

	T_d dist_to_nearest = Dist<T_d>(
		scat_data->get_data_coord(new_pnt, 0),
		scat_data->get_data_coord(new_pnt, 1), 0,
		scat_data->get_data_coord(m_knn_container[0], 0),
		scat_data->get_data_coord(m_knn_container[0], 1), 0);
	if (dist_to_nearest < EPSILON){
		return scat_data->get_data_value(m_knn_container[0]);
	}

	pseudo_insert(scat_data, num_points_current, new_pnt);

	plot_single_voronoi(scat_data, new_pnt, 1);

	//now we have a voronoi cell, we can disect it and do the interpolation 	
	for (T i = 0; i < m_pseudo_del.size(); i++){
		T q = m_pseudo_del[i];
		T m = (i == m_pseudo_del.size() - 1) ? m_pseudo_del[0] :
			m_pseudo_del[i + 1];

		T_d x_c(0), y_c(0);
		T_d r_2 = TriCircumcenter2d(scat_data->get_data_coord(new_pnt, 0),
			scat_data->get_data_coord(new_pnt, 1),
			scat_data->get_data_coord(q, 0),
			scat_data->get_data_coord(q, 1),
			scat_data->get_data_coord(m, 0),
			scat_data->get_data_coord(m, 1),
			x_c, y_c);
	}

}

template<class T, class T_d>
T VoroDel<T, T_d>::find_common_vertex(ScatData2D<T, T_d>*scat_data,
	const T p, const T m, const T q, T_d&xx, T_d&yy){
	//find the vertex v that is connected to m and q such that p lies inside 
	//the circumcircle of v-q-m and return the center of this circumcircle 
	//at xx,yy 
		
	for (T i = 0; i < m_del[q].size(); i++){
		T v = m_del[q][i];

		uint32_t id = get_index(v, m_del[m]);		

		if (id != std::numeric_limits<uint32_t>::max()){
			
			T_d r_2 = TriCircumcenter2d(scat_data->get_data_coord(v, 0),
				scat_data->get_data_coord(v, 1),
				scat_data->get_data_coord(q, 0),
				scat_data->get_data_coord(q, 1),
				scat_data->get_data_coord(m, 0),
				scat_data->get_data_coord(m, 1),
				xx, yy);

			T_d dist = Dist<T_d>(scat_data->get_data_coord(p, 0),
				scat_data->get_data_coord(p, 1), T_d(0), xx, yy, T_d(0));

			if (dist < r_2){				
				return v;
			}
		}
	}

	
	PRINT_ERROR("VoroDel::find_common_vertex() can not find the common vertex");
	

	
}
template<class T, class T_d>
inline void VoroDel<T, T_d>::plot_single_voronoi(ScatData2D<T, T_d>*scat_data,
	const T p, const bool with_del){

	std::string filename = STRINGIFY(OUTPUT_DIR)  "test.ps";

	std::fstream file(filename, std::ios::app);
	file.precision(30);

	double scale_x, scale_y, scale, shift_x, shift_y;

	double len = 1.2;

	scale_x = 6.5 / len;
	scale_y = 9.0 / len;

	if (scale_x < scale_y){
		scale = scale_x;
		shift_x = 0.9;
		shift_y = (11.0 - (len*scale)) / 2.0;
	}
	else {
		scale = scale_y;
		shift_x = (8.3 - (len*scale)) / 2.0;
		shift_y = 1.35;
	}


	T_d x_b(0), y_b(0), r_b(0);
	T_d xc_b(0), yc_b(0);
	
	T q0 = m_pseudo_del.back();
	T m0 = m_pseudo_del[0];

	r_b = TriCircumcenter2d(scat_data->get_data_coord(p, 0),
		scat_data->get_data_coord(p, 1),
		scat_data->get_data_coord(q0, 0),
		scat_data->get_data_coord(q0, 1),
		scat_data->get_data_coord(m0, 0),
		scat_data->get_data_coord(m0, 1),
		x_b, y_b);

	//find the common vertex between m0-q0
	T v0 = find_common_vertex(scat_data, p, q0, m0, xc_b, yc_b);

	for (T i = 0; i < m_pseudo_del.size(); i++){
		T q1 = m0;
		T m1 = (i == m_pseudo_del.size() - 1) ? m_pseudo_del[0] :
			m_pseudo_del[i + 1];

		T_d x_a(0), y_a(0), r_a(0);
		
		r_a = TriCircumcenter2d(scat_data->get_data_coord(p, 0),
			scat_data->get_data_coord(p, 1),
			scat_data->get_data_coord(q1, 0),
			scat_data->get_data_coord(q1, 1),
			scat_data->get_data_coord(m1, 0),
			scat_data->get_data_coord(m1, 1),
			x_a, y_a);

		T_d xc_a(0), yc_a(0);
		T v1 = find_common_vertex(scat_data, p, q1, m1, xc_a, yc_a);

		if (v1 == q0 || v0 == m1){
			//this section is triangle composed of (x_a,y_a), (x_b,y_b), (xc_a, yc_a)
			//since (xc_a, yc_a) := (xc_b, yc_b)
			if (!(v1 == q0 || v0 == m1)){
				PRINT_ERROR("VoroDel::plot_single_voronoi()????");
			}
			file << " newpath " << std::endl;
			file << " 0.02 setlinewidth" << std::endl;
			file << " " << double(rand()) / double(RAND_MAX)
				 << " " << double(rand()) / double(RAND_MAX) 
				 << " " << double(rand()) / double(RAND_MAX)  
				 << " setrgbcolor" << std::endl;
			file << x_a*scale << " " << y_a*scale;			
			file << " moveto " << std::endl;
			file << x_b*scale << " " << y_b*scale;
			file << " lineto " << std::endl;
			file << xc_a*scale << " " << yc_a*scale;
			file << " lineto " << std::endl;
			//file << "fill" << std::endl;
			file << " closepath" << std::endl;
			file << "stroke" << std::endl;

		}else{
			//this section is quad composed of (x_a,y_a), (x_b,y_b), (xc_b, yc_b),
			//(xc_a, yc_a)
			
			file << " newpath " << std::endl;
			file << " 0.02 setlinewidth" << std::endl;
			file << " " << double(rand()) / double(RAND_MAX)
				<< " " << double(rand()) / double(RAND_MAX)
				<< " " << double(rand()) / double(RAND_MAX)
				<< " setrgbcolor" << std::endl;
			file << x_a*scale << " " << y_a*scale;
			file << " moveto " << std::endl;
			file << x_b*scale << " " << y_b*scale;
			file << " lineto " << std::endl;
			file << xc_b*scale << " " << yc_b*scale;
			file << " lineto " << std::endl;
			file << xc_a*scale << " " << yc_a*scale;
			file << " lineto " << std::endl;
			file << "fill" << std::endl;
			file << " closepath" << std::endl;
			file << "stroke" << std::endl;
		}
		x_b = x_a;
		y_b = y_a;
		xc_b = xc_a;
		yc_b = yc_a;
		r_b = r_a;
		m0 = m1;
		q0 = q1;
		v0 = v1;
	}



	//Plot the boarder 
	file << " newpath " << std::endl;
	file << " 0.0001 setlinewidth" << std::endl;
	file << " 0 0 0  setrgbcolor" << std::endl;
	for (T i = 0; i < m_pseudo_del.size(); i++){
		T q = m_pseudo_del[i];
		T m = (i == m_pseudo_del.size() - 1) ? m_pseudo_del[0] :
			m_pseudo_del[i + 1];
		T_d x_c(0), y_c(0);
		T_d r_2 = TriCircumcenter2d(scat_data->get_data_coord(p, 0),
			scat_data->get_data_coord(p, 1),
			scat_data->get_data_coord(q, 0),
			scat_data->get_data_coord(q, 1),
			scat_data->get_data_coord(m, 0),
			scat_data->get_data_coord(m, 1),
			x_c, y_c);
		file << x_c*scale << " " << y_c*scale;
		if (i == 0){
			file << " moveto " << std::endl;
		}
		else{
			file << " lineto " << std::endl;
		}
	}

	//file << "fill" << std::endl;
	file << " closepath" << std::endl;
	file << "stroke" << std::endl;
	
	if (with_del){
		for (T i = 0; i < m_pseudo_del.size(); i++){
			T q = m_pseudo_del[i];

			file << " newpath " << std::endl;
			file << " 0.01 setlinewidth" << std::endl;
			file << " 0 0 0  setrgbcolor" << std::endl;

			file << scat_data->get_data_coord(p, 0)*scale << " " 
				 << scat_data->get_data_coord(p, 1)*scale;
			file << " moveto " << std::endl;

			file << scat_data->get_data_coord(q, 0)*scale << " " 
				 << scat_data->get_data_coord(q, 1)*scale;
			file << " lineto " << std::endl;

			file << " closepath" << std::endl;
			file << "stroke" << std::endl;

		}
	}

	file << scat_data->get_data_coord(p, 0)*scale << " "
		<< scat_data->get_data_coord(p, 1)*scale << " " << 0.006*scale
		<< " yellow_point" << std::endl;

}

#endif /*__VORO__*/
