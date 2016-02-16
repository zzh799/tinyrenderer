#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <limits>
#include <iostream>

#include "OpenNL_psm.h"
#include "tgaimage.h"
#include "geometry.h"

const int width  = 800;
const int height = 800;
const TGAColor white(255, 255, 255, 255);
const TGAColor red  (255,   0,   0, 255);
Vec3f light_dir(0,0,-1);

std::vector<Vec3f> verts;
std::vector<std::vector<int> > faces;
std::vector<bool> border;

// fills the border array, verts and faces are supposed to be already loaded (see load_obj())
void find_border_verts() {
    border = std::vector<bool>(verts.size(), false);
    std::vector<std::vector<int> > adj(verts.size());
    for (int i=0; i<(int)faces.size(); i++) {
        for (int k=0; k<3; k++) {
            adj[faces[i][k]].push_back(i*3+k);
        }
    }
    for (int i=0; i<(int)verts.size(); i++) {
        for (int j=0; j<(int)adj[i].size(); j++) { // loop through all faces incident to the vertex i and search for a halfedge without an opposite
            bool flag = false;
            int v1 = faces[adj[i][j]/3][(adj[i][j]%3 + 1)%3];
            for (int k=0; k<(int)adj[i].size(); k++) {
                int v2 = faces[adj[i][k]/3][(adj[i][k]%3 + 2)%3];
                if (v1==v2) {
                    flag = true;
                    break;
                }
            }
            if (!flag) {
                border[i] = true;
                break;
            }
        }
    }
}

// fills verts and faces arrays, supposes .obj file to have "f " entries without slashes
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

// well, a line is a line, no tricks here
void line(Vec2i a, Vec2i b, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(a.x-b.x)<std::abs(a.y-b.y)) {
        std::swap(a.x, a.y);
        std::swap(b.x, b.y);
        steep = true;
    }
    if (a.x>b.x) {
        std::swap(a, b);
    }
    for (int x=a.x; x<=b.x; x++) {
        float t = (x-a.x)/(float)(b.x-a.x);
        int y = a.y*(1.-t) + b.y*t;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

// given a pts triangle and a point P, find barycentric coordinates of P wrt to pts
Vec3f barycentric(Vec2i *pts, Vec2i P) {
    Vec3f u = cross(Vec3f(pts[2][0]-pts[0][0], pts[1][0]-pts[0][0], pts[0][0]-P[0]), Vec3f(pts[2][1]-pts[0][1], pts[1][1]-pts[0][1], pts[0][1]-P[1]));
    if (std::abs(u[2])<1) return Vec3f(-1,1,1); // triangle is degenerate, in this case return smth with negative coordinates
    return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
}

// triangle rasterizer with clamping to the screen size
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
    find_border_verts();
    TGAImage frame(width, height, TGAImage::RGB);

    std::vector<Vec3f> orig_verts = verts;

    // this loop smoothes the mesh separately for x, y and z dimensions
    for (int d=0; d<3; d++) {
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, verts.size());
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlBegin(NL_SYSTEM);
        nlBegin(NL_MATRIX);

        for (int i=0; i<(int)verts.size(); i++) {
            float scale = border[i] ? 100 : 1;
            nlBegin(NL_ROW);
            nlCoefficient(i, 1);
            nlRightHandSide(verts[i][d]);
            nlScaleRow(scale);
            nlEnd(NL_ROW);
        }

        for (unsigned int i=0; i<faces.size(); i++) {
            std::vector<int> &face = faces[i];
            for (int j=0; j<3; j++) {
                nlBegin(NL_ROW);
                nlCoefficient(face[ j     ],  1);
                nlCoefficient(face[(j+1)%3], -1);
                nlEnd(NL_ROW);
            }
        }

        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        nlSolve();

        for (int i=0; i<(int)verts.size(); i++) {
            verts[i][d] = nlGetVariable(i);
        }
    }

    // draw boundary of the original mesh
    for (unsigned int i=0; i<faces.size(); i++) {
        std::vector<int> &face = faces[i];
        Vec2i screen_coords[3];
        Vec3f world_coords[3];

        for (int j=0; j<3; j++) {
            world_coords[j]  = orig_verts[face[j]];
            screen_coords[j] = Vec2i((world_coords[j].x+1.)*width/2., (world_coords[j].y+1.)*height/2.);
        }

        for (int j=0; j<3; j++) {
            if (!border[face[j]] || !border[face[(j+1)%3]]) continue;
            line(screen_coords[j], screen_coords[(j+1)%3], frame, red);
        }
    }

    // draw smoothed mesh
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

    frame.write_tga_file("framebuffer.tga", false);
    return 0;
}

