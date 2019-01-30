#ifndef __IMAGE__
#define __IMAGE__

#include "tgaimage.h"

//T unsigned int 
//T_i signed int 
//T_d float/double 
template<class T, class T_i, class T_d>
class Image
{
public:
	Image();
	Image(T image_res[2], T_d image_x0[3], T_d image_xn[3], T_d image_normal[3],
		TGAImage *tga_image, std::string filename, bool flip_vertical = false, 
		bool flip_horizontal = false);
	
	void get_pixel_location(T i, T j, T_d pixel[3]){
		//TODO generalize this 
		//we assume for that the image plane is axis aligned 
		if (abs(abs(m_image_normal[2]) - 1) < EPSILON &&
			abs(m_image_normal[1]) <EPSILON &&
			abs(m_image_normal[0]) <EPSILON){
			pixel[0] = T_d(i)*m_s + m_image_x0[0];
			pixel[1] = T_d(j)*m_s + m_image_x0[1];
			pixel[2] = m_image_x0[2];			
		}
		else if (abs(m_image_normal[2]) < EPSILON &&
			     abs(abs(m_image_normal[1]) - 1) <EPSILON &&
			     abs(m_image_normal[0]) <EPSILON){
			pixel[0] = T_d(i)*m_s + m_image_x0[0];
			pixel[1] = m_image_x0[1];
			pixel[2] = T_d(j)*m_s + m_image_x0[2];
		}
		else if (abs(m_image_normal[2]) < EPSILON &&
			     abs(m_image_normal[1]) <EPSILON &&
				 abs(abs(m_image_normal[0]) - 1) <EPSILON){
			pixel[0] = m_image_x0[0];
			pixel[1] = T_d(i)*m_s + m_image_x0[1];			
			pixel[2] = T_d(j)*m_s + m_image_x0[2];
		}
		else{
			PRINT_ERROR("Image::get_pixel() can not figure the pixel location!!");
		}		
	}

	color_t get_pixel_color(T i, T j){
		if (i >= m_image_res[0] || j >= m_image_res[1]){
			PRINT_ERROR("Image::get_pixel_color() asking for pixel out of image range!!");
		}
		TGAColor tga_color;
		tga_color = m_image->get(i, j);

		color_t color;
		color.r = tga_color.r;
		color.g = tga_color.g;
		color.b = tga_color.b;
		color.a = tga_color.a;

		color.r /= 255.0;
		color.g /= 255.0;
		color.b /= 255.0;
		color.a /= 255.0;

		return color;
	}

	void get_image_res(T&res_x, T&res_y){
		res_x = m_image_res[0];
		res_y = m_image_res[1];
	}

	void set_pixel_color(T i, T j, color_t color){
		T my_flat_id = flat_id(i, j);
		TGAColor tga_color;

		if (color.r > 1.0 || color.g > 1.0 || color.b > 1.0 ||
			color.a > 1.0 || color.r < 0.0 || color.g < 0.0 || color.b < 0.0 ||
			color.a < 0.0){
			PRINT_ERROR("Image::set_pixel_color() setting un-clamp color at ("
				+ std::to_string(i) + ", " + std::to_string(j) + ")");
		}

		tga_color.r = 255.0*color.r;
		tga_color.g = 255.0*color.g;
		tga_color.b = 255.0*color.b;
		tga_color.a = 255.0*color.a;

		if (!m_image->set(i, j, tga_color)){
			PRINT_ERROR("Image::set_pixel_color() can not set pixel color at ("
				+ std::to_string(i)+ ", " + std::to_string(j) + ")");
		}

		//m_pixels_color[my_flat_id].r = color.r;
		//m_pixels_color[my_flat_id].g = color.g;
		//m_pixels_color[my_flat_id].b = color.b;
		//m_pixels_color[my_flat_id].a = color.a;
	}

	void get_image_normal(T_d&p0, T_d&p1, T_d&p2){
		//TODO generlize this 
		p0 = m_image_normal[0];
		p1 = m_image_normal[1];
		p2 = m_image_normal[2];		
	}
		
	void export_image(){
		if (m_flip_vertical){
			m_image->flip_vertically();
		}
		if (m_flip_horizontal){
			m_image->flip_horizontally();
		}
		std::cout << " Exporting image to " << m_filename << std::endl;
		m_image->write_tga_file(m_filename.c_str());
	}

	~Image(){
		
	};

private:
	inline T flat_id(T i, T j);
	TGAImage *m_image = NULL;

	T m_image_res[2];
	T_d m_image_x0[3], m_image_xn[3];
	T_d m_image_normal[3];

	color_t* m_pixels_color;

	bool m_flip_vertical, m_flip_horizontal;

	T_d m_s;//spacing

	std::string m_filename;
};

template<class T, class T_i, class T_d>
inline T Image<T, T_i, T_d>::flat_id(T i, T j){
	if (i >= m_image_res[0] || j >= m_image_res[1]){
		PRINT_ERROR("Image::flat_id() invalid flat id with (" + std::to_string(i)
			+ ", " + std::to_string(j)+") where image resolution is ("+
			std::to_string(m_image_res[0]) + ", " +
			std::to_string(m_image_res[1]) + ")");
	}

	return i + m_image_res[1] * j;
}

template<class T, class T_i, class T_d>
Image<T, T_i, T_d>::Image(T image_res[2], T_d image_x0[3], T_d image_xn[3],
	T_d image_normal[3], TGAImage *tga_image, std::string filename,
	bool flip_vertical /*= false*/,	bool flip_horizontal /*= false*/) :
	m_image(tga_image), m_filename(filename), m_flip_vertical(flip_vertical),
	m_flip_horizontal(flip_horizontal){
	m_image_res[0] = image_res[0];
	m_image_res[1] = image_res[1];

	if (m_image_res[0] != m_image_res[1]){
		PRINT_ERROR("Image::Image() image resolution should be the same in both direction!!");
	}

	m_image_x0[0] = image_x0[0];
	m_image_x0[1] = image_x0[1];
	m_image_x0[2] = image_x0[2];

	m_image_xn[0] = image_xn[0];
	m_image_xn[1] = image_xn[1];
	m_image_xn[2] = image_xn[2];

	m_image_normal[0] = image_normal[0];
	m_image_normal[1] = image_normal[1];
	m_image_normal[2] = image_normal[2];

	normalize_vector(m_image_normal[0], m_image_normal[1], m_image_normal[2]);

	//TODO fix this (this assumes that the image is on the x-y plane and facing 
	//the grid)
	if (abs(abs(m_image_normal[2]) - 1) < EPSILON &&
		abs(m_image_normal[1]) <EPSILON &&
		abs(m_image_normal[0]) <EPSILON){
		m_s = T_d(m_image_xn[0] - m_image_x0[0]) / T_d(m_image_res[0] - 1);
	
	}
	else if (abs(m_image_normal[2]) < EPSILON &&
		abs(abs(m_image_normal[1]) - 1) <EPSILON &&
		abs(m_image_normal[0]) <EPSILON){
		m_s = T_d(m_image_xn[0] - m_image_x0[0]) / T_d(m_image_res[0] - 1);
	
	}
	else if (abs(m_image_normal[2]) < EPSILON &&
		abs(m_image_normal[1]) <EPSILON &&
		abs(abs(m_image_normal[0]) - 1) <EPSILON){
		m_s = T_d(m_image_xn[2] - m_image_x0[2]) / T_d(m_image_res[0] - 1);	
	}
	else{
		PRINT_ERROR("Image::Image() image plane is not axis aligned!!");
	}

	

	m_pixels_color = 
		(color_t*)malloc(m_image_res[0] * m_image_res[1] * sizeof(color_t));
	

	
	/*if (m_image == NULL){
		static TGAImage image(m_image_res[0], m_image_res[1], TGAImage::RGBA);
		m_image = &image;
	}*/


}


#endif /*__IMAGE__*/
