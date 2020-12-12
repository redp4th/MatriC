#include "matric.h"
#include <cmath>
#include <random>

double dot(const Vector &a, const Vector &b) {
    if (a.length() != b.length())
        throw std::length_error("incompatible dimension");
    double result = 0;
    for (int i = 0; i < a.length(); ++i)
        result += a[i] * b[i];
    return result;
}

double euclidean_norm(const Vector &vec) {
    return std::sqrt(dot(vec, vec));
}

double max_norm(const Vector &vec) {
    double max = 0.0;
    for (int i = 0; i < vec.length(); ++i)
        max = std::abs(vec[i]) > max ? std::abs(vec[i]) : max;
    return max;
}

Vector one(int n) {
    return Vector(n, 1.0);
}

Vector elementary(int n, int j) {
    Vector v = Vector(n);
    v[j] = 1;
    return v;
}

Vector random_vector(int n) {
    Vector v(n);
    std::knuth_b engine;
    std::uniform_real_distribution<double> dist(-5, 5);
    for (int i = 0; i < n; ++i)
        v[i] = dist(engine);
    return v;
}

Vector operator*(double d, const Vector &vec) {
    return vec * d;
}

std::ostream &operator<<(std::ostream &out, const Vector &vec) {
    out << "[";
    for (int i = 0; i < vec.length(); ++i)
        out << vec[i] << (i == vec.length()-1 ? "" : ", ");
    out << "]";
    return out;
}

std::ostream &operator<<(std::ostream &out, Vector &&vec) {
    out << "[";
    for (int i = 0; i < vec.length(); ++i)
        out << vec[i] << (i == vec.length()-1 ? "" : ", ");
    out << "]";
    return out;
}

Matrix operator*(double d, const Matrix &m) {
    Matrix r(m);
    for (int i = 0; i < m.row; ++i)
        r[i] *= d;
    return r;
}

Matrix eye(int n) {
    Matrix m(n);
    for (int i = 0; i < n; ++i)
        m[i] = elementary(n, i);
    return m;
}

std::ostream &operator<<(std::ostream &out, const Matrix &m) {
    for (int i = 0; i < m.row; ++i)
        out << (i == 0 ? "[" : " ") << m[i] << (i == m.row - 1 ? "]" : ",\n");
    return out;
}

std::ostream &operator<<(std::ostream &out, Matrix &&m) {
    for (int i = 0; i < m.row; ++i)
        out << (i == 0 ? "[" : " ") << m[i] << (i == m.row - 1 ? "]" : ",\n");
    return out;
}

double rayleigh_quotient(const Matrix &matrix, const Vector &vec) {
    return dot(vec, matrix * vec) / dot(vec, vec);
}

Matrix givens(int n, int i, int j, double theta) {
    Matrix g = eye(n);
    g[i][i] = g[j][j] = std::cos(theta);
    g[i][j] = std::sin(theta);
    g[j][i] = -g[i][j];
    return g;
}

bool is_zero(double n, double tol) {
    return absolute(n) < tol;
}