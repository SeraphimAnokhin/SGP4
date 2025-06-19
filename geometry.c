#pragma once

typedef struct Matrix {
    double a11, a12, a13, a21, a22, a23, a31, a32, a33;
} Matrix;


typedef struct Vector {
    double x, y, z;
} Vector;


Vector mul_mat_vec(Matrix M, Vector v) {
    Vector ans = {M.a11 * v.x + M.a12 * v.y + M.a13 * v.z,
                  M.a21 * v.x + M.a22 * v.y + M.a23 * v.z,
                  M.a31 * v.x + M.a32 * v.y + M.a33 * v.z};
    return ans;
}


Matrix mul_mat_mat(Matrix A, Matrix B) {
    Matrix ans = {A.a11 * B.a11 + A.a12 * B.a21 + A.a13 * B.a31,    A.a11 * B.a12 + A.a12 * B.a22 + A.a13 * B.a32,    A.a11 * B.a13 + A.a12 * B.a23 + A.a13 * B.a33,
                  A.a21 * B.a11 + A.a22 * B.a21 + A.a23 * B.a31,    A.a21 * B.a12 + A.a22 * B.a22 + A.a23 * B.a32,    A.a21 * B.a13 + A.a22 * B.a23 + A.a23 * B.a33,
                  A.a31 * B.a11 + A.a32 * B.a21 + A.a33 * B.a31,    A.a31 * B.a12 + A.a32 * B.a22 + A.a33 * B.a32,    A.a31 * B.a13 + A.a32 * B.a23 + A.a33 * B.a33};
    return ans;
}


Vector mul_vec_num(Vector v, double n) {
    Vector ans = {v.x * n, v.y * n, v.z * n};
    return ans;
}


Vector add_vec(Vector a, Vector b) {
    Vector ans = {a.x + b.x, a.y + b.y, a.z + b.z};
    return ans;
}


Vector sub_vec(Vector a, Vector b) {
    Vector ans = {a.x - b.x, a.y - b.y, a.z - b.z};
    return ans;
}


Vector cross_product(Vector a, Vector b) {
    Vector ans = {a.y * b.z - a.z * b.y,
                  a.z * b.x - a.x * b.z,
                  a.x * b.y - a.y * b.x};
    return ans;
}


double dot_product(Vector a, Vector b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}


double vec_len(Vector v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}