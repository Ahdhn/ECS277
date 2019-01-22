#ifndef __RC__
#define __RC__

#include "grid.h"
#include "image.h"
#include "common.h"

//T unsigned int 
//T_i signed int 
//T_d float/double 
template <class T, class T_i, class T_d>
class RC
{
public:
	RC(Grid<T, T_d, 3>* grid, Image<T, T_i, T_d>*image, T num_samples);
	void run_raycasting();
	~RC(){};

private:
	Grid<T, T_d, 3>* m_grid;
	Image<T, T_i, T_d>*m_image;
	T m_num_samples; //num_samples per ray (uniformly sampled)
};

template <class T, class T_i, class T_d>
RC<T, T_i, T_d>::RC(Grid<T, T_d, 3>* grid, Image<T, T_i, T_d>*image, 
	T num_samples) :
	m_grid(grid), m_image(image), m_num_samples(num_samples)
{	
	
}

template <class T, class T_i, class T_d>
void RC<T, T_i, T_d>::run_raycasting(){
		
	T image_res[2];
	m_image->get_image_res(image_res[0], image_res[1]);

	//ray direction is the image plane normal
	T_d ray_dir[3];
	m_image->get_image_normal(ray_dir[0], ray_dir[1], ray_dir[2]);

	T_d sample[3];

	for (T i = 0; i < image_res[0]; i++){
		for (T j = 0; j < image_res[1]; j++){
			T_d ray_org[3];
			m_image->get_pixel(i, j, ray_org);

			
			T_d seg_start[3], seg_end[3];
			if (!m_grid->get_ray_grid_intersect(ray_org, ray_dir, seg_start, 
				seg_end)){
				continue;
			}

			T_d seg_len = Dist(seg_start[0], seg_start[1], seg_start[2],
				seg_end[0], seg_end[1], seg_end[2]);

			T_d spacing = seg_len / m_num_samples;

			color_t pixel_color;
			pixel_color.r = 0;
			pixel_color.g = 0;
			pixel_color.b = 0;
			
			for (T s = 0; s < m_num_samples; s++){
				
				for (T p = 0; p < 3; p++){
					sample[p] = seg_start[p] + spacing*s*ray_dir[p];
				}
				
				
				T_d f_value = m_grid->get_f_value_at_sample(sample);


				//use transfer function to get color and opacity for 

			}


		}
	}
}


#endif /*__RC__*/