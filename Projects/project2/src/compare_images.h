#ifndef __COMPARE_IMAGES__
#define __COMPARE_IMAGES__

template<class T, class T_i, class T_d>
T_d compare_images(Image<T, T_i, T_d>*image1, Image<T, T_i, T_d>*image2, 
	std::string filename){

	//the max error betrween two image is 1 i.e., a completely white image 

	T image_res1[2];
	image1->get_image_res(image_res1[0], image_res1[1]);

	T image_res2[2];
	image2->get_image_res(image_res2[0], image_res2[1]);

	if (image_res1[0] != image_res2[0] || image_res1[1] != image_res2[1]){
		PRINT_ERROR("compare_images() input images have different resolution!!");
	}

	TGAImage error_image(image_res1[0], image_res1[1], TGAImage::GRAYSCALE);

	T_d rms_error = 0;

	//T_d total_diff[image_res1[0]][image_res1[1]];
	std::vector<std::vector<T_d>>total_diff;

	T_d max_diff;
	//compute difference 
	for (T i = 0; i < image_res1[0]; i++){
		std::vector<T_d> diff_vec(image_res1[1]);
		for (T j = 0; j < image_res1[1]; j++){

			color_t c1 = image1->get_pixel_color(i, j);
			color_t c2 = image2->get_pixel_color(i, j);

			T_d dif_r = c1.r - c2.r;
			T_d dif_g = c1.g - c2.g;
			T_d dif_b = c1.b - c2.b;

			T_d diff = dif_r*dif_r + dif_g*dif_g + dif_b*dif_b;

			rms_error += diff;

			diff = sqrt(diff);

			//total_diff[i][j] = diff;
			diff_vec[j] = diff;
			
			max_diff = std::max(max_diff, diff);

		}
		total_diff.push_back(diff_vec);
	}

	//draw difference 
	for (T i = 0; i < image_res1[0]; i++){
		for (T j = 0; j < image_res1[1]; j++){
			TGAColor c;
			T_d diff = total_diff[i][j];
			c.r = 255.0 * diff / max_diff;
			c.g = 255.0 * diff / max_diff;
			c.b = 255.0 * diff / max_diff;
			c.a = 1.0;

			error_image.set(i, j, c);

		}
	}

	error_image.write_tga_file(filename.c_str());

	rms_error *= 1.0 / T_d(image_res1[0] * image_res1[1]);

	rms_error = sqrt(rms_error);

	return rms_error;
}

#endif /*__COMPARE_IMAGES__*/