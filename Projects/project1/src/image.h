#ifndef __IMAGE__
#define __IMAGE__

//T unsigned int 
//T_i signed int 
//T_d float/double 
template<class T, class T_i, class T_d>
class Image
{
public:
	Image();
	Image(T image_res[2], T_i image_x0[3], T_i image_xn[3]);
	
	void get_pixel(T i, T j, T_d pixel[3]){
		//TODO generalize this 

		pixel[0] = i*m_s + m_image_x0[0];
		pixel[1] = j*m_s + m_image_x0[1];
		pixel[2] = m_image_x0[2];
	}

	void get_image_res(T&res_x, T&res_y ){
		res_x = m_image_res[0];
		res_y = m_image_res[1];
	}

	void set_pixel_color(T i, T j, color_t color){
		T my_flat_id = flat_id(i, j);
		m_pixels_color[my_flat_id].r = color.r;
		m_pixels_color[my_flat_id].g = color.g;
		m_pixels_color[my_flat_id].b = color.b;
		m_pixels_color[my_flat_id].a = color.a;
	}

	void get_image_normal(T_d&p0, T_d&p1, T_d&p2){
		//TODO generlize this 
		p0 = 0;
		p1 = 0;
		p2 = 1;

		normalize_vector(p0, p1, p2);
	}

	~Image(){
		
	};

private:
	inline T flat_id(T i, T j);

	T m_image_res[2];
	T_i m_image_x0[3], m_image_xn[3];

	color_t* m_pixels_color;

	T_d m_s;//spacing
};

template<class T, class T_i, class T_d>
inline T Image<T, T_i, T_d>::flat_id(T i, T j){
	return i + m_image_res[1] * j;
}

template<class T, class T_i, class T_d>
Image<T, T_i, T_d>::Image(T image_res[2], T_i image_x0[3], T_i image_xn[3])
{
	m_image_res[0] = image_res[0];
	m_image_res[1] = image_res[1];

	m_image_x0[0] = image_x0[0];
	m_image_x0[1] = image_x0[1];
	m_image_x0[2] = image_x0[2];

	m_image_xn[0] = image_xn[0];
	m_image_xn[1] = image_xn[1];
	m_image_xn[2] = image_xn[2];

	m_s = T_d(m_image_xn[0] - m_image_x0[0]) / T_d(m_image_res[0]);

	m_pixels_color = 
		(color_t*)malloc(m_image_res[0] * m_image_res[1] * sizeof(color_t));


}


#endif /*__IMAGE__*/
