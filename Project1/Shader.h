#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include<iostream>
#include"line.h"

#define my_pi 3.14159265358

using namespace std;

class tov
{
public:
	
	Vec3f pos;
	Vec3f normal;
	tov(Vec3f position,Vec3f n);
	tov();
	~tov();
private:
	
};

tov::tov(Vec3f position,Vec3f n)
{
	pos = position;
	normal = n;

}

tov::tov() {

}

tov::~tov()
{
}

class v2f
{
public:
	

	Vec3f normal;
	Vec3f color;
	Vec3f viewpos;
	Vec2f tex_coords;
	Vec3f pos;
	Vec3f trinormal[3];
	Vec2f triuv[3];
	Vec3f tripos[3];


	v2f(Vec3f n,Vec3f c,Vec3f view,Vec2f tex_coord,Vec3f position);
	v2f();
	~v2f();
private:
	
	
};

v2f::v2f(Vec3f n, Vec3f c, Vec3f view, Vec2f tex_coord, Vec3f position)
{
	normal = n;
	color = c;
	viewpos = view;
	tex_coords = tex_coord;
	pos = position;

}
v2f::v2f() {

}

v2f::~v2f()
{
}

class Shader
{
public:
	Shader();
	~Shader();
	virtual v2f vertex_shader(tov o,TGAImage& image);
	virtual float fragment_shader(v2f v);
	TGAColor normal_fragment_shader(v2f v,int x,int y,Model model,Vec3f inter);
	void shadow_fragment_shader(Vec3f* pts, v2f v, float* zbuffer, Model model, TGAImage& image);
	Matrix viewport;
	Matrix projection;
	Matrix view;
private:

};

Shader::Shader()
{
}
Shader::~Shader()
{
}

v2f Shader::vertex_shader(tov o,TGAImage& image) {
	v2f v;
	Vec3f eye(1, 1, 3);
	Vec3f center(0, 0, 0);
	view = view_matrix(eye, center, Vec3f(0, 1, 0));
	
	projection = Matrix::identity(4);
	projection[3][2] = -1 / (eye - center).norm();

	float width = image.get_width();
	float height = image.get_height();

	viewport = viewport_matrix(width / 8, height / 8, 3 * width / 4, 3 * height / 4, 255);

	v.pos = homodiv(viewport * projection * view * vec3tovec4(o.pos));
	v.normal = o.normal;
	return v;
}

float Shader::fragment_shader(v2f v) {
	Vec3f light_dir = Vec3f(1, -1, 1).normalize();
	float intensity = max(float(0), light_dir * v.normal);
	return intensity;
}


inline TGAColor Shader::normal_fragment_shader(v2f v,int x,int y,Model model,Vec3f inter) {
	///cout << "t0" << endl;
	Vec3f light_dir = Vec3f(1, -1, 1).normalize();
	//cout << "t1" << endl;

	
	//uv插值
	Vec2f uv;
	uv.x = v.triuv[0].x * inter.x + v.triuv[1].x * inter.y + v.triuv[2].x * inter.z;
	uv.y = v.triuv[0].y * inter.x + v.triuv[1].y * inter.y + v.triuv[2].y * inter.z;
	//cout << "t2" << endl;

	//法线插值
	Vec3f nnormal;
	nnormal.x = inter.x * v.trinormal[0].x + inter.y * v.trinormal[1].x + inter.z * v.trinormal[2].x;
	nnormal.y = inter.x * v.trinormal[0].y + inter.y * v.trinormal[1].y + inter.z * v.trinormal[2].y;
	nnormal.z = inter.x * v.trinormal[0].z + inter.y * v.trinormal[1].z + inter.z * v.trinormal[2].z;
	//cout << "t3" << endl;


	//空间转换
	Vec3f l1 = v.tripos[1] - v.tripos[0];
	Vec3f l2 = v.tripos[2] - v.tripos[0];
	float delta_u1 = v.triuv[1].x - v.triuv[0].x;
	float delta_u2 = v.triuv[2].x - v.triuv[0].x;
	float delta_v1 = v.triuv[1].y - v.triuv[0].y;
	float delta_v2 = v.triuv[2].y - v.triuv[0].y;
	//cout << "t4" << endl;


	Vec3f t = (l2 * delta_v1 - l1 * delta_v2) * (1 / (delta_v1 * delta_u2 - delta_v2 * delta_u1));
	Vec3f b = (l1 * delta_u2 - l2 * delta_u1) * (1 / (delta_v1 * delta_u2 - delta_v2 * delta_u1));
	//cout << "t5" << endl;


	//正交化
	Vec3f T = t - nnormal * (t * nnormal);
	Vec3f B = b - nnormal * (b * nnormal) - T * (b * T);
	//cout << "t6" << endl;


	nnormal = nnormal.normalize();
	T = T.normalize();
	B = B.normalize();
	//cout << "t7" << endl;

	//构建TBN矩阵
	Matrix TBN = Matrix::identity(3);
	for (int idx = 0; idx < 3; idx++) {
		TBN[idx][0] = T[idx];
		TBN[idx][1] = B[idx];
		TBN[idx][2] = nnormal[idx];
	}

	//法线查询
	Vec3f final_normal = model.normal(uv);
	final_normal = (TBN * final_normal).normalize();

	float diffuse_color = max(float(0), final_normal * light_dir);

	TGAColor result_color = TGAColor(diffuse_color * 255, diffuse_color * 255, diffuse_color * 255, 255);

	return result_color;

}

void triangle(Vec3f* pts, v2f v, float* zbuffer, Model model, TGAImage& image) {
	//cout << "ok" << endl;
	//make bounding box
	float  min_x = float(image.get_width());
	float  min_y = float(image.get_height());

	float  max_x = 0;
	float  max_y = 0;


	for (int i = 0; i < 3; i++) {
		min_x = min(min_x, pts[i].x);
		max_x = max(max_x, pts[i].x);
		min_y = min(min_y, pts[i].y);
		max_y = max(max_y, pts[i].y);
	}
	int final_min_x = floor(min_x);
	int final_max_x = floor(max_x) + 1;
	int final_min_y = floor(min_y);
	int final_max_y = floor(max_y) + 1;


	//itera through all the pixels
	for (int x = final_min_x; x <= final_max_x; x++) {
		for (int y = final_min_y; y <= final_max_y; y++) {
			if (isinside(x + 0.5, y + 0.5, pts)) {
				Vec3f inter = barycentric(pts, x, y);

				//透视校正插值
				float z_interpolated = 1 / (inter.x / pts[0].z + inter.y / pts[1].z + inter.z / pts[2].z);

				Shader my_shader;
			//	Vec3f light_dir = Vec3f(1, -1, 1).normalize();
			//
			//	//uv插值
			//	Vec2f uv;
			//	uv.x = v.triuv[0].x * inter.x + v.triuv[1].x * inter.y + v.triuv[2].x * inter.z;
			//	uv.y = v.triuv[0].y * inter.x + v.triuv[1].y * inter.y + v.triuv[2].y * inter.z;
			//
			//
			//	//法线插值
			//	Vec3f nnormal;
			//	nnormal.x = inter.x * v.trinormal[0].x + inter.y * v.trinormal[1].x + inter.z * v.trinormal[2].x;
			//	nnormal.y = inter.x * v.trinormal[0].y + inter.y * v.trinormal[1].y + inter.z * v.trinormal[2].y;
			//	nnormal.z = inter.x * v.trinormal[0].z + inter.y * v.trinormal[1].z + inter.z * v.trinormal[2].z;
			//
			//	//空间转换
			//	Vec3f l1 = v.tripos[1] - v.tripos[0];
			//	Vec3f l2 = v.tripos[2] - v.tripos[0];
			//	float delta_u1 = v.triuv[1].x - v.triuv[0].x;
			//	float delta_u2 = v.triuv[2].x - v.triuv[0].x;
			//	float delta_v1 = v.triuv[1].y - v.triuv[0].y;
			//	float delta_v2 = v.triuv[2].y - v.triuv[0].y;
			//
			//	Vec3f t = (l2 * delta_v1 - l1 * delta_v2) * (1 / (delta_v1 * delta_u2 - delta_v2 * delta_u1));
			//	Vec3f b = (l1 * delta_u2 - l2 * delta_u1) * (1 / (delta_v1 * delta_u2 - delta_v2 * delta_u1));
			//
			//	//正交化
			//	Vec3f T = t - nnormal * (t * nnormal);
			//	Vec3f B = b - nnormal * (b * nnormal) - T * (b * T);
			//
			//	nnormal = nnormal.normalize();
			//	T = T.normalize();
			//	B = B.normalize();
			//
			//	//构建TBN矩阵
			//	Matrix TBN = Matrix::identity(3);
			//	for (int idx = 0; idx < 3; idx++) {
			//		TBN[idx][0] = T[idx];
			//		TBN[idx][1] = B[idx];
			//		TBN[idx][2] = nnormal[idx];
			//	}
			//
			//	//法线查询
			//	Vec3f final_normal = model.normal(uv);
			//	TGAColor outlook = model.diffuse(uv);
			//
			//	final_normal = (TBN * final_normal).normalize();
			//	float diffuse_color = max(float(0), final_normal * light_dir);
			//	TGAColor result_color = TGAColor(diffuse_color * 255, diffuse_color * 255, diffuse_color * 255, 255);

				TGAColor result_color = my_shader.normal_fragment_shader(v, x, y, model, inter);

				int a = get_index(x, y, image.get_width());
				if (z_interpolated < zbuffer[a]) {

					image.set(x,y, result_color);
					zbuffer[get_index(x, y, image.get_width())] = z_interpolated;

				}
			}
		}
	}
}

class depth_shader:Shader
{
public:
	depth_shader();
	depth_shader(Vec3f l, Vec3f c, Vec3f u);
	~depth_shader();
	virtual v2f vertex_shader(tov o, TGAImage& image);
	virtual float fragment_shader(v2f v, Vec3f inter, Vec3f* pts);
	void shadow_fragment_shader(Vec3f* pts, v2f v, float* zbuffer, float* depth_buffer, Shader my_shader, Model model, TGAImage& image);
private:
	Vec3f light_dir;
	Vec3f center;
	Vec3f up;
};

depth_shader::depth_shader()
{
}
depth_shader::depth_shader(Vec3f l,Vec3f c,Vec3f u)
{
	light_dir = l;
	center = c;
	up = u;
}

depth_shader::~depth_shader()
{
}

v2f depth_shader::vertex_shader(tov o, TGAImage& image) {
	v2f v;
	view = view_matrix(light_dir, center, up);
	projection = Matrix::identity(4);
	projection[3][2] = 0;

	float width = image.get_width();
	float height = image.get_height();

	viewport = viewport_matrix(width / 8, height / 8, 3 * width / 4, 3 * height / 4, 255);

	v.pos = homodiv(viewport * projection * view * vec3tovec4(o.pos));
	v.normal = o.normal;
	return v;
}
float depth_shader::fragment_shader(v2f v,Vec3f inter,Vec3f* pts) {

	float z_interpolated = 1 / (inter.x / pts[0].z + inter.y / pts[1].z + inter.z / pts[2].z);
	//float z_interpolated = inter.x * pts[0].z + inter.y * pts[1].z + inter.z * pts[0].z;

	return z_interpolated;
}



void triangle(Vec3f* pts, v2f v, float* zbuffer,depth_shader my_depth_shader, TGAImage& image) {
	//cout << "ok" << endl;
	//make bounding box
	float  min_x = float(image.get_width());
	float  min_y = float(image.get_height());

	float  max_x = 0;
	float  max_y = 0;


	for (int i = 0; i < 3; i++) {
		min_x = min(min_x, pts[i].x);
		max_x = max(max_x, pts[i].x);
		min_y = min(min_y, pts[i].y);
		max_y = max(max_y, pts[i].y);
	}
	int final_min_x = floor(min_x);
	int final_max_x = floor(max_x) + 1;
	int final_min_y = floor(min_y);
	int final_max_y = floor(max_y) + 1;


	//itera through all the pixels
	for (int x = final_min_x; x <= final_max_x; x++) {
		for (int y = final_min_y; y <= final_max_y; y++) {
			if (isinside(x + 0.5, y + 0.5, pts)) {
				Vec3f inter = barycentric(pts, x, y);
				
				float depth = my_depth_shader.fragment_shader(v, inter, pts);
				//std::cout << depth << endl;
				TGAColor result_color = TGAColor(depth, depth, depth, 255);

				int a = get_index(x, y, image.get_width());
				if (depth < zbuffer[a]) {

					image.set(x, y, result_color);
					zbuffer[get_index(x, y, image.get_width())] = depth;

				}
			}
		}
	}
}

void depth_shader::shadow_fragment_shader(Vec3f* pts, v2f v, float* zbuffer, float* depth_buffer, Shader my_shader, Model model, TGAImage& image) {
	//cout << "ok" << endl;
	//make bounding box
	float  min_x = float(image.get_width());
	float  min_y = float(image.get_height());

	float  max_x = 0;
	float  max_y = 0;


	for (int i = 0; i < 3; i++) {
		min_x = min(min_x, pts[i].x);
		max_x = max(max_x, pts[i].x);
		min_y = min(min_y, pts[i].y);
		max_y = max(max_y, pts[i].y);
	}
	int final_min_x = floor(min_x);
	int final_max_x = floor(max_x) + 1;
	int final_min_y = floor(min_y);
	int final_max_y = floor(max_y) + 1;


	//itera through all the pixels
	for (int x = final_min_x; x <= final_max_x; x++) {
		for (int y = final_min_y; y <= final_max_y; y++) {
			if (isinside(x + 0.5, y + 0.5, pts)) {
				Vec3f inter = barycentric(pts, x, y);
				
				//透视校正插值
				float z_interpolated = 1 / (inter.x / pts[0].z + inter.y / pts[1].z + inter.z / pts[2].z);
				//z_interpolated = 1;
				Vec3f light_dir = Vec3f(1, -1, 1).normalize();

				//uv插值
				Vec2f uv;
				uv.x = v.triuv[0].x * inter.x + v.triuv[1].x * inter.y + v.triuv[2].x * inter.z;
				uv.y = v.triuv[0].y * inter.x + v.triuv[1].y * inter.y + v.triuv[2].y * inter.z;


				//法线插值
				Vec3f nnormal;
				nnormal.x = inter.x * v.trinormal[0].x + inter.y * v.trinormal[1].x + inter.z * v.trinormal[2].x;
				nnormal.y = inter.x * v.trinormal[0].y + inter.y * v.trinormal[1].y + inter.z * v.trinormal[2].y;
				nnormal.z = inter.x * v.trinormal[0].z + inter.y * v.trinormal[1].z + inter.z * v.trinormal[2].z;

				//空间转换
				Vec3f l1 = v.tripos[1] - v.tripos[0];
				Vec3f l2 = v.tripos[2] - v.tripos[0];
				float delta_u1 = v.triuv[1].x - v.triuv[0].x;
				float delta_u2 = v.triuv[2].x - v.triuv[0].x;
				float delta_v1 = v.triuv[1].y - v.triuv[0].y;
				float delta_v2 = v.triuv[2].y - v.triuv[0].y;

				Vec3f t = (l2 * delta_v1 - l1 * delta_v2) * (1 / (delta_v1 * delta_u2 - delta_v2 * delta_u1));
				Vec3f b = (l1 * delta_u2 - l2 * delta_u1) * (1 / (delta_v1 * delta_u2 - delta_v2 * delta_u1));

				//正交化
				Vec3f T = t - nnormal * (t * nnormal);
				Vec3f B = b - nnormal * (b * nnormal) - T * (b * T);

				nnormal = nnormal.normalize();
				T = T.normalize();
				B = B.normalize();

				//构建TBN矩阵
				Matrix TBN = Matrix::identity(3);
				for (int idx = 0; idx < 3; idx++) {
					TBN[idx][0] = T[idx];
					TBN[idx][1] = B[idx];
					TBN[idx][2] = nnormal[idx];
				}
				//法线查询
				Vec3f final_normal = model.normal(uv);
				TGAColor outlook = model.diffuse(uv);

				Vec3f lig_dir(1, 1, 0);

				final_normal = (TBN * final_normal).normalize();
				float diffuse_color = max(float(0), final_normal * light_dir);

				//阴影判定
				tov my_temp_tov;

				Matrix temp_vec4(4, 1);
				temp_vec4[0][0] = x;
				temp_vec4[1][0] = y;
				temp_vec4[2][0] = z_interpolated;
				temp_vec4[3][0] = 1;

				temp_vec4 = (my_shader.viewport * my_shader.projection * my_shader.view).inverse()*temp_vec4;
				my_temp_tov.pos = Vec3f(temp_vec4[0][0], temp_vec4[1][0], temp_vec4[2][0]);

				Vec3f depth_pos = vertex_shader(my_temp_tov, image).pos;

				int depth_idx = get_index(depth_pos.x,depth_pos.y ,image.get_width());
				int count = 0;
				if (depth_idx < 0) {
					count += 1;
					cout << count << endl;
					continue;
				}
					
				float shadow = 0.3 + 0.7 * (depth_buffer[depth_idx] > depth_pos.z - 10);
				diffuse_color *= shadow;

				TGAColor result_color = TGAColor(diffuse_color * 255, diffuse_color * 255, diffuse_color * 255, 255);

				int a = get_index(x, y, image.get_width());
				if (z_interpolated < zbuffer[a]) {

					image.set(x, y, result_color);
					zbuffer[get_index(x, y, image.get_width())] = z_interpolated;

				}
			}
		}
	}
}

class ao_shader:Shader
{
public:
	ao_shader();
	~ao_shader();
	void ao_fragment_shader(Vec3f* pts, float* depth_buffer,TGAImage& image);
private:

};

ao_shader::ao_shader()
{
}

ao_shader::~ao_shader()
{
}

void ao_shader::ao_fragment_shader(Vec3f* pts, float* depth_buffer, TGAImage& image) {
	//cout << "ok" << endl;
//make bounding box
	float  min_x = float(image.get_width());
	float  min_y = float(image.get_height());

	float  max_x = 0;
	float  max_y = 0;


	for (int i = 0; i < 3; i++) {
		min_x = min(min_x, pts[i].x);
		max_x = max(max_x, pts[i].x);
		min_y = min(min_y, pts[i].y);
		max_y = max(max_y, pts[i].y);
	}
	int final_min_x = floor(min_x);
	int final_max_x = floor(max_x) + 1;
	int final_min_y = floor(min_y);
	int final_max_y = floor(max_y) + 1;


	//itera through all the pixels
	for (int x = final_min_x; x <= final_max_x; x++) {
		for (int y = final_min_y; y <= final_max_y; y++) {
			if (isinside(x + 0.5, y + 0.5, pts)) {
				Vec3f inter = barycentric(pts, x, y);

				//透视校正插值
				float z_interpolated = 1 / (inter.x / pts[0].z + inter.y / pts[1].z + inter.z / pts[2].z);
				//z_interpolated = 1;

				//ao计算
				float total = 0;

				for (float angel = 0; angel < 2 * my_pi - 1e-4; angel += my_pi / 8) {
					Vec2f dir = Vec2f(cos(angel), sin(angel));
					Vec2f uv = Vec2f(x, y);

					//total += get_max_angel(dir, uv, depth_buffer, image);

					Vec2f current_pos;
					float distance;
					int idx;
					float vertical_dis;
					float max_angel = 0;
					for (int i = 1; i < 10; i++) {
						current_pos = uv + dir * i;
						if (current_pos.x > image.get_width() || current_pos.y > image.get_height() || current_pos.x < 0 || current_pos.y < 0) break;
						distance = Vec3f((current_pos - uv).x, (current_pos - uv).y, 0).norm();
						idx = get_index(int(current_pos.x), int(current_pos.y), image.get_width());

						vertical_dis = depth_buffer[idx] - depth_buffer[get_index(int(uv.x), int(uv.y), image.get_width())];

						max_angel = max(max_angel, atanf(-vertical_dis / distance));
					}

					total += max_angel;
				}
					float final_intensity = (8 * my_pi - total) / (8 * my_pi);
					final_intensity = pow(final_intensity, 1.1);

					TGAColor result_color = TGAColor(final_intensity * 255, final_intensity * 255, final_intensity * 255, 255);

					int a = get_index(x, y, image.get_width());

					if (z_interpolated <= depth_buffer[a]) {
						image.set(x, y, result_color);
					}
				}
			}
		}
	}
