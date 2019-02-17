#ifndef __RENDERER__
#define __RENDERER__

#include "grid.h"
#include "image.h"
#include "common.h"
#include <vector>
//T unsigned int 
//T_i signed int 
//T_d float/double 
template <class T, class T_i, class T_d>
class Renderer
{
public:
	Renderer(Grid<T, T_d, 3>* grid, T samples_per_cell);

	template<typename ColorTF, typename AlphaTF, typename SkipF>
	void run_raycasting(Image<T, T_i, T_d>*image, 
		ColorTF color_trans_func, AlphaTF alpha_trans_func,
		SkipF threshold_skip_f, INTERPOL_TYPE type = INTERPOL_TYPE::TRILINEAR);


	template<typename ColorTF, typename AnalyticalFunc>

	void slice(Image<T, T_i, T_d>*image,
		ColorTF color_trans_func, 
		AnalyticalFunc analytical_func,
		T_d dist = 0.5,
		bool compute_error = true,
		INTERPOL_TYPE type = INTERPOL_TYPE::TRILINEAR);

	~Renderer(){};

private:
	Grid<T, T_d, 3>* m_grid;	
	T m_samples_per_cell; //num_samples per ray per cell (uniformly sampled)
};

template <class T, class T_i, class T_d>
Renderer<T, T_i, T_d>::Renderer(Grid<T, T_d, 3>* grid, T samples_per_cell) :
	m_grid(grid), m_samples_per_cell(samples_per_cell){	
	
}

template<class T, class T_i, class T_d>
template<typename ColorTF, typename AnalyticalFunc>
void Renderer<T, T_i, T_d>::slice(Image<T, T_i, T_d>*image, 
	ColorTF color_trans_func,
	AnalyticalFunc analytical_func,
	T_d dist/* = 0.5*/,
	bool compute_error/* = true*/, 
	INTERPOL_TYPE type /*= INTERPOL_TYPE::TRILINEAR*/){

	//dist is a normalized parameter that decide the 'depth' of the cutting 
	//dist should be betwee 0 and 1. For every pixel, we shoot a ray to check 
	//if the pixel can see the grid. We return the line segment inside the grid 
	//and use dist to be the (normalized) distance we walk to get the depth 
	//of the slicing plane

	color_t background_color = m_grid->get_background_color();
	background_color.a = 0.8;

	T_d grad[3];

	T image_res[2];
	image->get_image_res(image_res[0], image_res[1]);
		

	//view direction
	T_d view_dir[3];
	image->get_image_normal(view_dir[0], view_dir[1], view_dir[2]);

	T_d sample[3];

	for (T i = 0; i < image_res[0]; i++){
		for (T j = 0; j < image_res[1]; j++){

			T_d pixel[3];
			image->get_pixel_location(i, j, pixel);
			
			//we check if the pixel can see the grid by shooting a ray and check
			//if it is intersect the grid 
			T_d seg_start[3], seg_end[3];
			if (!m_grid->get_ray_grid_intersect(pixel, view_dir, seg_start,
				seg_end)){				
				image->set_pixel_color(i, j, background_color);

				if (compute_error){					
					image->accumelate_error(i, j, background_color);
				}
				continue;
			}


			T_d seg_len = Dist(seg_start[0], seg_start[1], seg_start[2],
				seg_end[0], seg_end[1], seg_end[2]);

			color_t pixel_color;
			pixel_color.r = 0;
			pixel_color.g = 0;
			pixel_color.b = 0;
			pixel_color.a = 0;


			for (T p = 0; p < 3; p++){
				//sample[p] = seg_start[p] + spacing*s*ray_dir[p];
				sample[p] = seg_end[p] - dist*view_dir[p];
			}


			if (!m_grid->is_inside_grid(sample)){
				//if the sample touches the edges of the grid 			
				image->set_pixel_color(i, j, background_color);
				continue;
			}

			T_d f_value = m_grid->get_f_value_at_sample(sample, grad, type,
				true);

			color_t local_color = color_trans_func(f_value);

			local_color.clamp();			

			image->set_pixel_color(i, j, local_color);
			

			if (compute_error){
				//ground truth 
				T_d f_value_gt = analytical_func(sample[0], sample[1], sample[2]);
				color_t gt_color = color_trans_func(f_value_gt);
				gt_color.clamp();
				image->accumelate_error(i, j, gt_color);
			}



		}
	}

	image->export_image();
}

template <class T, class T_i, class T_d>
template<typename ColorTF, typename AlphaTF, typename SkipF>
void Renderer<T, T_i, T_d>::run_raycasting(Image<T, T_i, T_d>*image,
	ColorTF color_trans_func, 
	AlphaTF alpha_trans_func, SkipF threshold_skip_f, 
	INTERPOL_TYPE type /*= INTERPOL_TYPE::TRILINEAR*/){
		
	std::cout << "\n Raycaster started with for ";
	if (type == INTERPOL_TYPE::TRILINEAR){
		std::cout << "Trilinear case" << std::endl;
	}
	else{
		std::cout << "Trilcubic case" << std::endl;
	}
		
	T image_res[2];
	image->get_image_res(image_res[0], image_res[1]);

	T image_size = image_res[0] * image_res[1];
	//ray direction is the image plane normal
	T_d ray_dir[3];
	image->get_image_normal(ray_dir[0], ray_dir[1], ray_dir[2]);

	T_d sample[3];
	const color_t background_color = m_grid->get_background_color();
	uint32_t num_samples_per_ray = 
		m_samples_per_cell* (m_grid->get_resolution(2) - 1);


	uint32_t num_processed = 0;
	std::cout << " Renderer::run_raycasting 0% ";
	std::vector<int> percentage(1, 0);


	for (T i = 0; i < image_res[0]; i++){
		for (T j = 0; j < image_res[1]; j++){
			
			num_processed++;
			int percetnage_processed = int(100 * double(num_processed) / 
				double(image_size));
			if (percetnage_processed % 10 == 0
				&& percetnage_processed != percentage.back()){
				std::cout << percetnage_processed << "%  ";
				percentage.push_back(percetnage_processed);
			}
			
			T_d ray_org[3];
			image->get_pixel_location(i, j, ray_org);

			//std::cout << "i = "<<i << " j= " << j << std::endl;
			T_d seg_start[3], seg_end[3];
			if (!m_grid->get_ray_grid_intersect(ray_org, ray_dir, seg_start, 
				seg_end)){
				color_t color = background_color;
				color.a = 0.8;
				image->set_pixel_color(i, j, color);
				continue;
			}
		
			T_d seg_len = Dist(seg_start[0], seg_start[1], seg_start[2],
				seg_end[0], seg_end[1], seg_end[2]);

			T_d spacing = seg_len / num_samples_per_ray;

			color_t pixel_color;
			pixel_color.r = 0;
			pixel_color.g = 0; 
			pixel_color.b = 0;
			pixel_color.a = 0;
												
			for (T s = 0; s < num_samples_per_ray; s++){
			//for (T s = 0; s < num_samples_per_ray; s++){
								
				for (T p = 0; p < 3; p++){
					//sample[p] = seg_start[p] + spacing*s*ray_dir[p];
					sample[p] = seg_end[p] - spacing*s*ray_dir[p];
				}

				if (!m_grid->is_inside_grid(sample)){
					//if the sample touches the edges of the grid 
					continue;
				}
				//interpolate to get the function value at this sample
				T_d grad[3];				
				T_d f_value = m_grid->get_f_value_at_sample(sample, grad, type);								
				
				if (threshold_skip_f(f_value)){
					continue;
				}

				

				//iso/imaginary surface diffuse color 
				color_t local_color = color_trans_func(f_value);

				//use transfer function to get color and opacity for 				
				local_color.a = alpha_trans_func(f_value);
								
				//shade the samples 
				T_d dist_to_image_sq = Dist(sample[0], sample[1], sample[2],
					                     ray_org[0], ray_org[1], ray_org[2]);
			
				color_t sample_color = m_grid->shade_point(sample, 
					dist_to_image_sq, local_color, ray_dir, grad);

				sample_color.a = local_color.a;

				//accumelate from beck to front								
				/*sample_color *= sample_color.a;
				pixel_color *= (1.0 - sample_color.a);
				pixel_color += sample_color;*/
								
							
				{
					//correcting opacity 
					sample_color.a = 1.0 - pow(1.0 - sample_color.a, 
						20.0/num_samples_per_ray);

					T_d acc_alpha = sample_color.a +
						(1.0 - sample_color.a)*pixel_color.a;
					pixel_color *= (1.0 - sample_color.a);
					pixel_color += sample_color;
					pixel_color.a = acc_alpha;

				}
				
				
				/*{
					if (s == 25){
						pixel_color = local_color;
						pixel_color.a = 1.0;

						break;
					}
				}*/

				//clamp 
				pixel_color.clamp();
												
				if (pixel_color.r > 0.99 &&
					pixel_color.g > 0.99 &&
					pixel_color.b > 0.99){
					//if the pixel is white enough, quit 
					//std::cout << ".";
					break;
				}
			}			
									
			image->set_pixel_color(i, j, pixel_color);
		}
	}

	std::cout << std::endl;
	image->export_image();
}


#endif /*__RENDERER__*/