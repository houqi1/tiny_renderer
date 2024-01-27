#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include<iostream>
#include"line.h"
#include"Shader.h"

using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
Model* model = NULL;
Vec3f eye(1, 1, 3);
Vec3f center(0, 0, 0);
Vec3f light_dir(1, 1, 0);



int main(int argc, char** argv) {
    if (2 == argc) {
        model = new Model(argv[1]);
    }
    else {
        model = new Model("obj/african_head/african_head.obj");
    }
    float width = 800;
    float height = 800;

    Shader my_shader;
    depth_shader my_depth_shader(eye, center, Vec3f(0, 1, 0));

    TGAImage image(width, height, TGAImage::RGB);
    TGAImage depth_image(width, height, TGAImage::RGB);

    float* zbuffer = new float[width * height];
   // for (int i = 0; i < image.get_width() * image.get_height(); ++i) {
   //     zbuffer[i] = std::numeric_limits<float>::infinity();
   // }
   //     for (int j = 0; j < model->nfaces(); j++) {
   //         std::vector<int> face = model->face(j);
   //         Vec3f pts1[3];
   //
   //         // Vec3f world_coord[3];
   //         Vec3f intensity;
   //         tov my_tov;
   //         v2f my_v2f;
   //         for (int i = 0; i < 3; i++) {
   //
   //             my_tov.pos = model->vert(face[i]);
   //
   //             pts1[i] = my_depth_shader.vertex_shader(my_tov, image).pos;
   //
   //
   //         }
   //         //cout << pts1[1] << endl;
   //         triangle(pts1, my_v2f, zbuffer, my_depth_shader, depth_image);
   //     }
   //     depth_image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
   //     depth_image.write_tga_file("depth_buffer.tga");
   //
   //     ao_shader my_ao_shader;
   //
   //     float* zbuffer2 = new float[width * height];
   //     for (int i = 0; i < image.get_width() * image.get_height(); ++i) {
   //         zbuffer2[i] = std::numeric_limits<float>::infinity();
   //     }
   //
   //     for (int j = 0; j < model->nfaces(); j++) {
   //         std::vector<int> face = model->face(j);
   //         Vec3f pts1[3];
   //
   //         // Vec3f world_coord[3];
   //         Vec3f intensity;
   //         tov my_tov;
   //         for (int i = 0; i < 3; i++) {
   //
   //             my_tov.pos = model->vert(face[i]);
   //
   //             pts1[i] = my_depth_shader.vertex_shader(my_tov, image).pos;
   //
   //
   //         }
   //         //cout << pts1[1] << endl;
   //         cout << j << endl;
   //         my_ao_shader.ao_fragment_shader(pts1, zbuffer,image);
   //     }
   //
   //     image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
   //     image.write_tga_file("ao_output.tga");
   //
   //     return 0;

        float* zbuffer1 = new float[width * height];
        for (int i = 0; i < image.get_width() * image.get_height(); ++i) {
            zbuffer1[i] = std::numeric_limits<float>::infinity();
        }

        for (int j = 0; j < model->nfaces(); j++) {
            std::vector<int> face = model->face(j);
            Vec3f pts1[3];

            tov my_tov;
            v2f my_v2f;
            for (int i = 0; i < 3; i++) {

                my_tov.pos = model->vert(face[i]);
                pts1[i] = my_shader.vertex_shader(my_tov, image).pos;

                my_v2f.tripos[i] = model->vert(face[i]);
                my_v2f.triuv[i] = model->uv(j, i);

                my_v2f.trinormal[i] = model->normal(j, i);
            }
            cout << j << endl;
            triangle(pts1, my_v2f, zbuffer1, *model, image);
        }

        image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
        image.write_tga_file("output.tga");
        delete model;

        return 0;
    }
