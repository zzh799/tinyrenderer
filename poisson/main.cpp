#include <cassert>
#include "OpenNL_psm.h"
#include "tgaimage.h"

void nlScaleRow(NLdouble s);
int main() {
    TGAImage a,b,m,result;
    a.read_tga_file("a.tga");
    b.read_tga_file("b.tga");
    m.read_tga_file("m.tga");

    assert(a.get_width()  == b.get_width()  && m.get_width()  == a.get_width());
    assert(a.get_height() == b.get_height() && m.get_height() == a.get_height());

    int w = a.get_width();
    int h = a.get_height();
    result = TGAImage(w, h, TGAImage::GRAYSCALE);

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, w*h);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);

    for (int i=0; i<w*h; i++) {
        if (m.get(i%w, i/w)[0]<128) continue;
        nlBegin(NL_ROW);
        nlCoefficient(i, 1);
        nlRightHandSide(a.get(i%w, i/w)[0]);
        nlScaleRow(100);
        nlEnd(NL_ROW);
    }

    for (int i=0; i<w-1; i++) {
        for (int j=0; j<h-1; j++) {
            for (int d=0; d<2; d++) {
                nlBegin(NL_ROW);
                int v1 = b.get(i,   j    )[0];
                int v2 = b.get(i+d, j+1-d)[0];
                nlCoefficient( i   + j     *w,  1);
                nlCoefficient((i+d)+(j+1-d)*w, -1);
                nlRightHandSide(v1-v2);
                nlScaleRow(10);
                nlEnd(NL_ROW);
            }
        }
    }

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    nlSolve();

    for (int i=0; i<w*h; i++) {
        float v = std::max(0., std::min(255., nlGetVariable(i)));
        result.set(i%w, i/w, TGAColor((unsigned char)v));
    }

    result.write_tga_file("result.tga", true);
    return 0;
}

