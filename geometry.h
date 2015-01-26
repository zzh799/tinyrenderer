#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <iomanip>

template<size_t DimRows,size_t DimCols,typename Number> class mat;

template <size_t Dim,typename Number> struct vec {
    Number items[Dim];
    static const size_t DimN=Dim;
    vec<Dim,Number> normalize() const { return (*this)/norm(); }
    Number norm() const { return std::sqrt((*this)*(*this)); }

    size_t maxElementPos() const { //найти номер максимального элемента
        size_t ret=0;
        for (size_t i=Dim;--i;) {
            if (items[i]>items[ret]) {
                ret=i;
            }
        }
        return ret;
    }

    static vec<Dim,Number > fill(const Number& val=0) { // построить вектор, заполненный константой TODO переделать в конструктор
        vec<Dim, Number> ret;
        for (size_t i=Dim; i--; ret[i]=val);
        return ret;
    }

    bool operator!=(const vec<Dim,Number>&v) { // TODO переделать в >=, пригодится для растеризатора
        for (size_t i=Dim; i--; ) {
            if(v[i]!=items[i])
                return true;
        }
        return false;
    }

    Number& operator [](size_t index) {
        assert(index<Dim);
        return items[index];
    }

    const Number& operator [](size_t index) const {
        assert(index<Dim);
        return items[index];
    }
};

template<size_t Dim,typename Number>vec<Dim,Number > operator+(vec<Dim,Number > lhs, const vec<Dim,Number >& rhs) {
    for (size_t i=Dim; i--; lhs[i]+=rhs[i]);
    return lhs;
}

template<size_t Dim,typename Number> vec<Dim,Number> operator^(const vec<Dim,Number >&lhs, const vec<Dim,Number >& rhs) { // векторное умножение TODO переделать в cross()
    vec<Dim,Number> ret;
    for (size_t i=Dim; i--; ) {
        mat<Dim,3,Number> temp;
        temp[0]=vec<Dim,Number>::fill(0);
        temp[0][i]=1;
        temp[1]=lhs;
        temp[2]=rhs;
        ret[i]=temp.det();
    }
    return ret;
}

template<size_t Dim,typename Number> Number operator*(const vec<Dim,Number >&lhs, const vec<Dim,Number >& rhs) {
    Number ret=0;
    for (size_t i=Dim; i--; ret+=lhs[i]*rhs[i]);
    return ret;
}

template<size_t Dim,typename Number>vec<Dim,Number > operator-(vec<Dim,Number > lhs, const vec<Dim,Number >& rhs) {
    for (size_t i=Dim; i--; lhs[i]-=rhs[i]);
    return lhs;
}

template<size_t Dim,typename Number>vec<Dim,Number > operator*(vec<Dim,Number > lhs, const Number& rhs) {
    for (size_t i=Dim; i--; lhs[i]*=rhs);
    return lhs;
}

template<size_t Dim,typename Number>vec<Dim,Number > operator/(vec<Dim,Number > lhs, const Number& rhs) {
    for (size_t i=Dim; i--; lhs[i]/=rhs);
    return lhs;
}

template<size_t len,size_t Dim, typename Number> vec<len,Number > embed(const vec<Dim,Number> &v,const Number& fill=1) { // погружение вектора
    vec<len,Number> ret = vec<len,Number>::fill(fill);
    for (size_t i=Dim; i--; ret[i]=v[i]);
    return ret;
}

template<size_t len,size_t Dim, typename Number> vec<len,Number > proj(const vec<Dim,Number> &v) { //проекция вектора
    vec<len,Number> ret;
    for (size_t i=len; i--; ret[i]=v[i]);
    return ret;
}

template<size_t Dim,typename Number> std::ostream& operator<<(std::ostream& out,const vec<Dim,Number >& v) {
    out<<"{ ";
    for (size_t i=0; i<Dim; i++) {
        out<<std::setw(6)<<v[i]<<" ";
    }
    out<<"} ";
    return out;
}

#include "vec-simple.h"

/////////////////////////////////////////////////////////////////////////////////

typedef vec<2,float> Vec2f;
typedef vec<2,int>   Vec2i;
typedef vec<3,float> Vec3f;
typedef vec<3,int>   Vec3i;
typedef mat<4,4,float> Matrix;


#endif //__GEOMETRY_H__
