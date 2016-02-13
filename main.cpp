#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <limits>
#include <iostream>
#include "tgaimage.h"
#include "geometry.h"

const int width  = 800;
const int height = 800;
const TGAColor white(255, 255, 255, 255);
Vec3f light_dir(0,0,-1);

std::vector<Vec3f> verts;
std::vector<std::vector<int> > faces;

void load_obj(const char *filename) {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v[i];
            verts.push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            std::vector<int> f;
            int idx;
            iss >> trash;
            while (iss >> idx) {
                idx--; // in wavefront obj all indices start at 1, not zero
                f.push_back(idx);
            }
            faces.push_back(f);
        }
    }
    std::cerr << "# v# " << verts.size() << " f# "  << faces.size() << std::endl;
}

Vec3f barycentric(Vec2i *pts, Vec2i P) {
    Vec3f u = cross(Vec3f(pts[2][0]-pts[0][0], pts[1][0]-pts[0][0], pts[0][0]-P[0]), Vec3f(pts[2][1]-pts[0][1], pts[1][1]-pts[0][1], pts[0][1]-P[1]));
    if (std::abs(u[2])<1) return Vec3f(-1,1,1); // triangle is degenerate, in this case return smth with negative coordinates
    return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
}

void triangle(Vec2i *pts, TGAImage &image, TGAColor color) {
    Vec2i bboxmin(image.get_width()-1,  image.get_height()-1);
    Vec2i bboxmax(0, 0);
    Vec2i clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0,        std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    Vec2i P;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f bc_screen  = barycentric(pts, P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            image.set(P.x, P.y, color);
        }
    }
}


int main() {
    load_obj("face.obj");
    TGAImage frame(width, height, TGAImage::RGB);

    for (unsigned int i=0; i<faces.size(); i++) {
        std::vector<int> &face = faces[i];
        Vec2i screen_coords[3];
        Vec3f world_coords[3];

        for (int j=0; j<3; j++) {
            world_coords[j]  = verts[face[j]];
            screen_coords[j] = Vec2i((world_coords[j].x+1.)*width/2., (world_coords[j].y+1.)*height/2.);
        }

        Vec3f n = cross(world_coords[2]-world_coords[0], world_coords[1]-world_coords[0]);
        n.normalize();
        float intensity = n*light_dir;
        if (intensity>0) {
            triangle(screen_coords, frame, white*intensity);
        }
    }

    frame.flip_vertically();
    frame.write_tga_file("framebuffer.tga");
    return 0;
}

