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
	
	T find_next(const std::vector<T>& list, const T q0, const T v0);
		
	T find_common_vertex(ScatData2D<T, T_d>*scat_data, const T p,
		const T m, const T q, T_d&xx, T_d&yy, const T guide);
	
	T_d get_poly_area(const std::vector<T_d>poly);
	
	T_d triangle_area(T_d x1, T_d y1, T_d z1, T_d x2, T_d y2, T_d z2,
		T_d x3, T_d y3, T_d z3);

	inline void plot_pseudo(ScatData2D<T, T_d>*scat_data, const T p,
		const bool with_del);

	std::vector<T> m_pseudo_del;

	std::vector<int>m_knn_container;
	std::vector<int>m_knn_search;
	std::vector<std::vector<T>> m_del;
	std::vector<bool> m_boundary;
	std::vector<T_d> m_angles;
	std::vector<std::vector<T_d>> m_pseudo_voro;
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
		T p_id = (current_data.size() == 0) ? p : current_data[p];

		//if (p_id < 4 ){
		if (p_id < 8){
			//if it is connected to one of the four points of the bounding box
			m_boundary[p_id] = true;
		}

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

		//check co linearity 
		T_d arr = triangle_area(
			scat_data->get_data_coord(p, 0), scat_data->get_data_coord(p, 1), 0,
			scat_data->get_data_coord(q, 0), scat_data->get_data_coord(q, 1), 0,
			scat_data->get_data_coord(m, 0), scat_data->get_data_coord(m, 1), 0);

		if (abs(arr) < EPSILON){
			continue;
		}


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

	file << "/green_point      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 1 0 setrgbcolor" << std::endl;	
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
	//file << " fill" << std::endl;
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

	
	
	if (voronoi){

		std::vector<T_d> voro_vertices;	

		//bool skip_boundary = true;
		for (T p = 0; p < m_del.size(); p++){	

			//if (p < 4){
			if (p < 8){
				continue;
			}
			if (skip_boundary && m_boundary[p]){
				continue;
			}

			if (m_del[p].size() == 0){ continue; }
						
			voro_vertices.clear();

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
				voro_vertices.push_back(x_c);
				voro_vertices.push_back(y_c);

				file << x_c*scale << " " << y_c*scale;
				if (i == 0){
					file << " moveto " << std::endl;
				}
				else{
					file << " lineto " << std::endl;
				}
			}

			file << "fill" << std::endl;
			file << " closepath" << std::endl;
			file << "stroke" << std::endl;		

			//draw the cell boundaries 
			file << " newpath " << std::endl;
			file << " 0.01 setlinewidth" << std::endl;
			file << " 0 0 0 setrgbcolor" << std::endl;
			for (T i = 0; i < m_del[p].size(); i++){

				file << voro_vertices[i * 2] * scale << " "
					<< voro_vertices[i * 2 + 1] * scale;

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

			//draw dot 
			file << scat_data->get_data_coord(p, 0)*scale << " "
				<< scat_data->get_data_coord(p, 1)*scale << " " << 0.004*scale
				<< " green_point" << std::endl;
		}

	

	}


	if (delaunay){
		for (T p = 0; p < m_del.size(); p++){

			if (m_del[p].size() == 0){ continue; }

			for (T i = 0; i < m_del[p].size(); i++){
				T q = m_del[p][i];

				//if (p < 4){
				if (p < 8){
					//the first four points are the corner points and we draw 
					//them as lines (not triangles) to avoid looping around ourselves 
					//and end up with triangles that are not part of the mesh
					//continue;
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


	for (T p = 0; p < m_del.size(); p++){
		T_d r_point = 0.001;
		
		if (m_boundary[p]){		
			T_d r_disk = 0.004;
			file << scat_data->get_data_coord(p, 0)*scale << " "
				<< scat_data->get_data_coord(p, 1)*scale << " " << r_disk*scale
				<< " red_dot" << std::endl;
		}
			

		file << scat_data->get_data_coord(p, 0)*scale << " "
			<< scat_data->get_data_coord(p, 1)*scale << " " << r_point*scale
			<< " point" << std::endl;

		/*file << "/Times-Roman findfont" << std::endl;
		file << "0.04 scalefont" << std::endl;
		file << "0 0 0 setrgbcolor" << std::endl;
		file << "setfont" << std::endl;
		file << scat_data->get_data_coord(p, 0)*scale << " " 
			<< scat_data->get_data_coord(p, 1)*scale
			<< " moveto" << std::endl;
		file << "(" << p;
		file << ") " << "show" << std::endl;*/

	}
}

template<class T, class T_d>
void VoroDel<T, T_d>::pseudo_insert(ScatData2D<T, T_d>*scat_data,
	const T num_points_current, const T new_pnt){

	//psuedo insert the new point i.e., build its delaunay and store it in
	//m_pseudo_del

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

	m_pseudo_voro.clear();

	for (T i = 0; i < m_pseudo_del.size();i++){
		std::vector<T_d> vv;
		m_pseudo_voro.push_back(vv);
	}
	

	//now we have a voronoi cell, we can disect it and do the interpolation 	
	//m_pseudo_voro holds a disected voronoi cell of the pseudo-insert site	
	T q0 = m_pseudo_del.back();
	for (T i = 0; i < m_pseudo_del.size(); i++){
		T q = m_pseudo_del[i];
		T next_q = (i == m_pseudo_del.size() - 1) ? m_pseudo_del[0] :
			m_pseudo_del[i + 1];

		T next_q0 = q;

		//initial voronoi vertex
		T_d x_c0(0), y_c0(0);
		T_d r_c0 = TriCircumcenter2d(scat_data->get_data_coord(new_pnt, 0),
			scat_data->get_data_coord(new_pnt, 1),
			scat_data->get_data_coord(q0, 0),
			scat_data->get_data_coord(q0, 1),
			scat_data->get_data_coord(q, 0),
			scat_data->get_data_coord(q, 1),
			x_c0, y_c0);

		m_pseudo_voro[i].push_back(x_c0);
		m_pseudo_voro[i].push_back(y_c0);

		//intermediate voronoi vertices
		T_d x_c1(0), y_c1(0), r_c1;
		T guide = (m_pseudo_del.size() == 3) ? next_q : 
			scat_data->get_num_data();
		T v0 = find_common_vertex(scat_data, new_pnt, q0, q, x_c1, y_c1, guide);
		m_pseudo_voro[i].push_back(x_c1);
		m_pseudo_voro[i].push_back(y_c1);

		while (v0 != next_q){
			T n0 = find_next(m_del[q], v0, q0);
			if (n0 == std::numeric_limits<T>::max()){
				break;
			}
			q0 = v0;
			v0 = n0;

			r_c1 = TriCircumcenter2d(scat_data->get_data_coord(q, 0),
				scat_data->get_data_coord(q, 1),
				scat_data->get_data_coord(q0, 0),
				scat_data->get_data_coord(q0, 1),
				scat_data->get_data_coord(v0, 0),
				scat_data->get_data_coord(v0, 1),
				x_c1, y_c1);

			m_pseudo_voro[i].push_back(x_c1);
			m_pseudo_voro[i].push_back(y_c1);
		}


		//terminal voronoi vertex
		T_d x_c2(0), y_c2(0);
		T_d r_c2 = TriCircumcenter2d(scat_data->get_data_coord(new_pnt, 0),
			scat_data->get_data_coord(new_pnt, 1),
			scat_data->get_data_coord(q, 0),
			scat_data->get_data_coord(q, 1),
			scat_data->get_data_coord(next_q, 0),
			scat_data->get_data_coord(next_q, 1),
			x_c2, y_c2);

		m_pseudo_voro[i].push_back(x_c2);
		m_pseudo_voro[i].push_back(y_c2);

		q0 = next_q0;

	}
	
	if (false){
		auto voro_shade = [](T p){ color_t color; 
		color.r = color.g = color.b = 0;  return color;  };


		plot_pseudo(scat_data, new_pnt, false);
	}



	//actually compute the sibson interpolation 
	T_d sib = 0;
	T_d total_area = 0;
	for (T i = 0; i < m_pseudo_del.size(); i++){
		T q = m_pseudo_del[i];
		T_d fi = scat_data->get_data_value(q);

		T_d ai = get_poly_area(m_pseudo_voro[i]);
		total_area += ai;

		sib += ai*fi;

	}

	sib /= total_area;

	return sib;

}

template<class T, class T_d>
T VoroDel<T, T_d>::find_common_vertex(ScatData2D<T, T_d>*scat_data,
	const T p, const T m, const T q, T_d&xx, T_d&yy, const T guide){
	//find the vertex v that is connected to m and q such that p lies inside 
	//the circumcircle of v-q-m and return the center of this circumcircle 
	//at xx,yy 
	
	std::vector<T> candidate;
	std::vector<T_d> candidate_circ;
		
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
				candidate.push_back(v);
				candidate_circ.push_back(xx);
				candidate_circ.push_back(yy);
				//return v;
			}
		}
	}

	//if more than one circumcirlce that p lies inside, we return the one 
	//where the center is closer to p (heuristic but okay)
	if (candidate.size() == 1){
		return candidate[0];
	}
	else if (candidate.size() == 2){

		if (candidate[0] == guide){ return candidate[0]; }
		if (candidate[1] == guide){ return candidate[1]; }

		T_d d0 = Dist<T_d>(scat_data->get_data_coord(p, 0),
			scat_data->get_data_coord(p, 1), T_d(0),
			candidate_circ[0], candidate_circ[1], T_d(0));

		T_d d1 = Dist<T_d>(
			scat_data->get_data_coord(p, 0),
			scat_data->get_data_coord(p, 1), T_d(0),
			candidate_circ[2], candidate_circ[3], T_d(0));
		
		if (d0 > d1){
			return candidate[0];
		}
		else{
			return candidate[1];
		}		
	}
	
	//This werid configuration, delaunay says that and edge is shared 
	//two triangles only!!!
	//PRINT_ERROR("VoroDel::find_common_vertex() can not find the common vertex");
		
}
template<class T, class T_d>
T VoroDel<T, T_d>::find_next(const std::vector<T>& list, const T q0, 
	const T v0){
	//in the sort list list, find the next element, n0, after v0 such that 
	//q0 is before v0 and n0 is after n0
	//That means if list contains n0->v0->q0
	//it should returns n0 
	//Also if the list contains q0->v0->n0 
	//it should returns n0 
	//because in both cases v0 is in the middle and q0 is before v0 and n0 is
	//after v0 i.e., we define an order by n0, v0, q0 that is independent on 
	//the how list is sorted. The only case this funtion fail is when there is 
	//no link between q0 and v0 in list

	for (T i = 0; i < list.size(); i++){
		if (list[i] == v0){
			//there is only two cases:
			//I) q0 is behind v0
			//II) q0 is after v0
			T next_i = (i == list.size() - 1) ? 0 : i + 1;
			T before_i = (i == 0) ? list.size() - 1 : i - 1;

			if (q0 == list[next_i]){
				//I)
				T next_next_i = (next_i == list.size() - 1) ? 0 : next_i + 1;
				return list[next_next_i];
			
			}
			else if (q0 == list[before_i]){
				//II)
				T before_before_i = (before_i == 0) ? list.size() - 1 : before_i - 1;
				return list[before_before_i];
			}
			else{
				//ok i fail
				return std::numeric_limits<T>::max();
				PRINT_ERROR("VoroDel::find_next() q0 is not before nor after v0!!");
			}
			
		}
		
	}
	
	PRINT_ERROR("VoroDel::find_next() can not find v0!!");
}
template<class T, class T_d>
inline void VoroDel<T, T_d>::plot_pseudo(ScatData2D<T, T_d>*scat_data,
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

	
	for (T i = 0; i < m_pseudo_voro.size(); i++){
		file << " newpath " << std::endl;
		file << " 0.02 setlinewidth" << std::endl;
		file << " " << double(rand()) / double(RAND_MAX)
			<< " " << double(rand()) / double(RAND_MAX)
			<< " " << double(rand()) / double(RAND_MAX)
			<< " setrgbcolor" << std::endl;
		for (T j = 0; j < m_pseudo_voro[i].size() >> 1; ++j){
			file << m_pseudo_voro[i][j * 2] * scale << " "
				<< m_pseudo_voro[i][j * 2 + 1] * scale;
			if (j == 0){
				file << " moveto " << std::endl;
			}
			else{
				file << " lineto " << std::endl;;
			}
		}
		file << "fill" << std::endl;
		file << " closepath" << std::endl;
		file << "stroke" << std::endl;
	}
	

	//draw the voronoi site location
	file << scat_data->get_data_coord(p, 0)*scale << " "
		<< scat_data->get_data_coord(p, 1)*scale << " " << 0.004*scale
		<< " red_dot" << std::endl;


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

	

}

template<class T , class T_d>
T_d VoroDel<T, T_d>::get_poly_area(const std::vector<T_d>poly){
	//poly contain coordinates (sorted pairs of x,y) of a polygon 
	//we compute the area of the polygon as sum of triangles area assuming 
	//it is convex polygon 
	size_t num_points = poly.size() >> 1;
	if (num_points < 3){
		PRINT_ERROR("VoroDel::get_poly_area() poly has less three vertices!!");
	}	

	T_d area = 0;
	T p0 = 0;
	T p1 = 1;	
	for (T i = 2; i < poly.size() >> 1; i++){
		T p2 = i;
		area += triangle_area(poly[p0 * 2], poly[p0 * 2 + 1], 0,
			poly[p1 * 2], poly[p1 * 2 + 1], 0,
			poly[p2 * 2], poly[p2 * 2 + 1], 0);
		p1 = p2;
	}


	return area;

}

template<class T, class T_d>
T_d VoroDel<T, T_d>::triangle_area(T_d x1, T_d y1, T_d z1,
	                               T_d x2, T_d y2, T_d z2, 
					               T_d x3, T_d y3, T_d z3){
	T_d l1 = sqrt(Dist(x1, y1, T_d(0), x2, y2, T_d(0)));
	T_d l2 = sqrt(Dist(x3, y3, T_d(0), x2, y2, T_d(0)));
	T_d l3 = sqrt(Dist(x3, y3, T_d(0), x1, y1, T_d(0)));

	T_d p = (l1 + l2 + l3) / 2.0;

	if ((p - l1) <= 0 || (p - l2) <= 0 || (p - l3) <= 0){
		//std::cout << " Warning at TriArea = 0" << std::endl;
		return 0;
		//system("pause");
	}

	return sqrt(p*(p - l1)*(p - l2)*(p - l3));

}

#endif /*__VORO__*/
