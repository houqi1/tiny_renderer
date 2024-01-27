#include"line.h"
#include<stdio.h>
#include <vector>
#include <cmath>
#include"geometry.h"
#include<iostream>
//#include"Shader.h"
using namespace std;


void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
    bool steep = false;
    if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    for (int x = x0; x <= x1; x++) {
        float t = (x - x0) / (float)(x1 - x0);
        int y = y0 * (1. - t) + y1 * t;
        if (steep) {
            image.set(y, x, color);
        }
        else {
            image.set(x, y, color);
        }
    }
}
bool isinside(float x, float y, Vec2i* pts) {
    Vec3i cuurent_point = { int(x),int(y),1 };

    int sum_score = 0;
    for (int i = 0; i < 3; i++) {
        Vec3i temp_point1 = { pts[i].x,pts[i].y,1};
        Vec3i temp_point2 = { pts[(i + 1) % 3].x,pts[(i + 1) % 3].y,1 };

        Vec3i v1 = cuurent_point - temp_point1;
        Vec3i v2 = temp_point2 - temp_point1;

        int cp = (v2 ^ v1).z;

        if (cp > 0) sum_score += 1;
    }
    if (sum_score == 3) return true;
    return false;
}

bool isinside(float x, float y, Vec3f* pts) {
    Vec3f cuurent_point = { x,y,1 };

    int sum_score = 0;
    for (int i = 0; i < 3; i++) {

        Vec3f v1 = cuurent_point - pts[i];
        Vec3f v2 = pts[(i+1)%3] - pts[i];

        int cp = (v2 ^ v1).z;

        if (cp >= 0) sum_score += 1;
    }
    if (sum_score == 3) return true;
    return false;
}

int get_index(int x, int y, int width) {
    return x + y * width;
}

//重心插值
Vec3f barycentric(Vec3f* pts, float x,float y) {

    
    Vec3f barycentric_value;
    float alpha = (-(x - pts[1].x) * (pts[2].y - pts[1].y) + (y - pts[1].y) * (pts[2].x - pts[1].x)) / (-(pts[0].x - pts[1].x) * (pts[2].y - pts[1].y) + (pts[0].y - pts[1].y) * (pts[2].x - pts[1].x));
    float beta = (-(x - pts[2].x) * (pts[0].y - pts[2].y) + (y - pts[2].y) * (pts[0].x - pts[2].x)) / (-(pts[1].x - pts[2].x) * (pts[0].y - pts[2].y) + (pts[1].y - pts[2].y) * (pts[0].x - pts[2].x));
    float gamma = 1 - alpha - beta;
    
    barycentric_value = { alpha,beta,gamma };
    
    return barycentric_value;

}

void triangle(Vec2i* pts, TGAImage& image, TGAColor color) {
//make bounding box
    int  min_x = image.get_width();
    int  min_y = image.get_height();

    int  max_x = 0;
    int  max_y = 0;

    
    for (int i = 0; i < 3; i++) {
        min_x = min(min_x, pts[i].x);
        max_x = max(max_x, pts[i].x);
        min_y = min(min_y, pts[i].y);
        max_y = max(max_y, pts[i].y);
    }
    int* zbuffer =new int[image.get_width() * image.get_height()];

    for (int i = 0; i < image.get_width() * image.get_height(); ++i) {
        zbuffer[i] = std::numeric_limits<int>::infinity();
    }


   //itera through all the pixels
    for (int x = min_x; x <= max_x; x++) {
        for (int y = min_y; y <= max_y; y++) {
            if (isinside(x + 0.5, y + 0.5, pts)) {
                image.set(x, y, color);
            }
        }
    }
}
void triangle(Vec3f* pts, float* zbuffer,TGAImage& image, TGAColor color) {
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
    for (int x = min_x; x <= max_x; x++) {
        for (int y = min_y; y <= max_y; y++) {
            if (isinside(x + 0.5, y + 0.5, pts)) {
                Vec3f inter = barycentric(pts, x, y);
                //透视校正插值

                float z_interpolated = 1/(inter.x /pts[0].z + inter.y /pts[1].z + inter.z /pts[2].z);
                z_interpolated *= (inter.x * pts[0].z + inter.y * pts[1].z + inter.z * pts[2].z);
                //z_interpolated = 1;
                
               //cout << z_interpolated << endl;
                int a = get_index(x, y, image.get_width());
                if (z_interpolated < zbuffer[a]) {
                    image.set(x, y, color);
                    zbuffer[get_index(x, y, image.get_width())] = z_interpolated;
                }
            }
        }
    }
}

void triangle(Vec3f* pts, float* zbuffer, Vec3f intensity,TGAImage& image, TGAColor color) {
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
                z_interpolated *= (inter.x * pts[0].z + inter.y * pts[1].z + inter.z * pts[2].z);
                //z_interpolated = 1;
                

                //光照强度插值
                float real_intensity = inter.x * intensity.x + inter.y * intensity.y + inter.z * intensity.z;

               //cout << z_interpolated << endl;
                int a = get_index(x, y, image.get_width());
                if (z_interpolated < zbuffer[a]) {
                    image.set(x, y, TGAColor(real_intensity*255,real_intensity*255,real_intensity*255,255));
                    zbuffer[get_index(x, y, image.get_width())] = z_interpolated;

                }
            }
        }
    }
}



Vec3f world2screen(Vec3f world,TGAImage& image) {
    Vec3f screen_coord;
    float width = image.get_width();
    float height = image.get_height();
    screen_coord = { (world.x + 1) * width / 2,(world.y + 1) * height / 2,world.z };
    return screen_coord;
}

Matrix vec3tovec4(Vec3f v) {
    Matrix m(4, 1);

    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1;

    return m;

}

Matrix viewport_matrix(float x,float y,float width, float height,float depth) {
    Matrix m = Matrix::identity(4);
    m[0][0] = width/2;
    m[1][1] = height/2;
    m[2][2] = depth/ 2;

    m[0][3] = x+width / 2;
    m[1][3] = y+height /2;
    m[2][3] = depth / 2;
    

    return m;
}

//不含正交投影的透视矩阵
Matrix get_projection(float z_near, float zfar) {
    Matrix projection = Matrix::identity(4);

    projection[0][0] = z_near;
    projection[1][1] = z_near;
    projection[2][2] = z_near + zfar;
    projection[2][3] = -zfar * z_near;
    projection[3][2] = 1;
    projection[3][3] = 0;

    return projection;
}

Vec3f homodiv(Matrix m) {
    Vec3f homo;
    homo = { m[0][0] / m[3][0],m[1][0] / m[3][0],m[2][0] / m[3][0] };
    return homo;
}

Matrix view_matrix(Vec3f eye, Vec3f center, Vec3f up) {

    Vec3f temp_look = (center - eye).normalize();//z
    Vec3f temp_up = (temp_look^up).normalize();//x
    Vec3f right = (temp_up ^ temp_look).normalize();//y

    Matrix view_transform = Matrix::identity(4);
    for (int i = 0; i < 3; i++) {
        view_transform[0][i] = temp_up[i];
        view_transform[1][i] = right[i];
        view_transform[2][i] = temp_look[i];
        view_transform[i][3] = -center[i];
    }
    return view_transform;
}

//模型矩阵（暂时）
Matrix model_matrix() {
    return Matrix::identity(4);
}

float get_max_angel(Vec2f dir, Vec2f uv, float* depth_buffer,TGAImage image) {
    Vec2f current_pos;
    float distance;
    int idx;
    float vertical_dis;
    float max_angel = 0;
    for (int i = 1; i < 1000; i++) {
        current_pos = uv + dir * i;
        if (current_pos.x > image.get_width() || current_pos.y > image.get_height() || current_pos.x < 0 || current_pos.y < 0) return max_angel;
        distance = Vec3f((current_pos - uv).x, (current_pos - uv).y, 0).norm();
        idx = get_index(int(current_pos.x), int(current_pos.y), image.get_width());
    
        vertical_dis = depth_buffer[idx] - depth_buffer[get_index(int(uv.x), int(uv.y), image.get_width())];
    
        max_angel = max(max_angel, -atan(vertical_dis / distance));
    }
    return max_angel;

}
