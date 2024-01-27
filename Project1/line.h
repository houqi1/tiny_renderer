#include"tgaimage.h"
#include"geometry.h"
#include"model.h"

void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color);

void triangle(Vec2i* pts, TGAImage& image, TGAColor color);
void triangle(Vec3f* pts,float* zbuffer, TGAImage& image, TGAColor color);
void triangle(Vec3f* pts, float* zbuffer, Vec3f intensity, TGAImage& image, TGAColor color);

int get_index(int x, int y, int width);

bool isinside(float x, float y, Vec3f* pts);
bool isinside(float x, float y, Vec2i* pts);
Vec3f barycentric(Vec3f* pts, float x, float y);

Vec3f world2screen(Vec3f world, TGAImage& image);

Matrix vec3tovec4(Vec3f v);
Matrix viewport_matrix(float x,float y,float width, float height,float depth);

Matrix get_projection(float z_near, float zfar);
Vec3f homodiv(Matrix m);

Matrix view_matrix(Vec3f eye, Vec3f center, Vec3f up);

Matrix model_matrix();

float get_max_angel(Vec2f dir, Vec2f uv, float* depth_buffer,TGAImage image);
