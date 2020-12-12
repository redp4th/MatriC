/*
 *  +------------------------------------------+
 *  |     __  ___        __         _  ______  |
 *  |    /  |/  /____ _ / /_ _____ (_)/ ____/  |
 *  |   / /|_/ // __ `// __// ___// // /       |
 *  |  / /  / // /_/ // /_ / /   / // /___     |
 *  | /_/  /_/ \__,_/ \__//_/   /_/ \____/     |
 *  |                                          |
 *  +------------------------------------------+
 *  A simple matrix library by Jiangnan Tang @ bupt.
 *
 *  Created at: 2020/12/10
 *
 *  Written for use of Numerical Analysis experiment.
 */

#ifndef MATRIC_H
#define MATRIC_H

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <utility>
#include <iostream>
#include <exception>
#include <initializer_list>

/* 
 * Forward declaration
 */
struct Vector;
struct Matrix;

double dot(const Vector &, const Vector &);
double euclidean_norm(const Vector &);
double max_norm(const Vector &);

inline double identity(double n) { return n; }
inline double absolute(double n) { return n > 0 ? n : -n; }

struct Vector {
    /*
        For now, it's just a real vector.
        Support for entries been complex is on its way.
        (but not until the author have finished advanced courses
        such as complex analysis.)
    */
    double *data;
    int len;

    Vector(): data(nullptr), len(0) {}

    Vector(int n): data(new double[n]), len(n) {
        for (int i = 0; i < n; ++i)
            this->data[i] = 0;
    }

    Vector(int n, double d): data(new double[n]), len(n) {
        for (int i = 0; i < n; ++i)
            this->data[i] = d;
    }

    Vector(std::initializer_list<double> l): data(new double[l.size()]), len(l.size()) {
        int j = 0;
        for (auto d: l)
            this->data[j++] = d;
    }

    Vector(const Vector &vec): len(vec.length()), data(new double[vec.length()]) {
        std::memcpy(this->data, vec.data, sizeof(double) * this->len);
    }

    ~Vector() {
        delete[] data;
    }

    Vector &operator=(const Vector &other) {
        if (len != other.length()) {
            delete[] this->data;
            this->data = new double[other.length()];
            this->len = other.length();
        }
        std::memcpy(this->data, other.data, sizeof(double) * this->len);
        return *this;
    }

    Vector &operator=(std::initializer_list<double> l) {
        if (len != l.size()) {
            delete[] this->data;
            this->data = new double[l.size()];
            this->len = l.size();
        }
        int j = 0;
        for (auto d: l)
            this->data[j++] = d;
        return *this;
    }

    double &operator[](int i) {
        if (!(0 <= i < this->len))
            throw std::out_of_range("index out of bound");
        return this->data[i];
    }

    const double &operator[](int i) const {
        if (!(0 <= i < this->len))
            throw std::out_of_range("index out of bound");
        return this->data[i];
    }
    
    Vector operator+(const Vector &other) const {
        if (this->len != other.length())
            throw std::length_error("incompatible dimension");
        Vector result(len);
        for (int i = 0; i < len; ++i)
            result[i] = this->data[i] + other[i];
        return result;
    }

    Vector operator-() const {
        Vector v(this->len);
        for (int i = 0; i < this->len; ++i)
            v[i] = -this->data[i];
        return v;
    }

    Vector operator-(const Vector &other) const {
        if (this->len != other.length())
            throw std::length_error("incompatible dimension");
        Vector result(this->len);
        for (int i = 0; i < this->len; ++i)
            result[i] = this->data[i] - other[i];
        return result;
    }

    double operator*(const Vector &other) const {
        return dot(*this, other);
    }

    Vector operator*(double n) {
        Vector v(this->length());
        for (int i = 0; i < len; ++i)
            v[i] = this->data[i] * n;
        return v;
    }

    Vector operator/(double n) {
        if (n == 0)
            throw std::invalid_argument("division by 0");
        Vector v(this->length());
        for (int i = 0; i < len; ++i)
            v[i] = this->data[i] / n;
        return v;
    }

    Vector operator^(double n) {
        Vector v(this->length());
        for (int i = 0; i < this->len; ++i)
            v[i] = std::pow(this->data[i], n);
        return v;
    }

    Vector &operator+=(const Vector &other) {
        if (this-> len != other.length())
            throw std::length_error("incompatible dimension");
        for (int i = 0; i < len; ++i)
            this->data[i] += other[i];
        return *this;
    }


    Vector &operator-=(const Vector &other) {
        if (this-> len != other.length())
            throw std::length_error("incompatible dimension");
        for (int i = 0; i < len; ++i)
            this->data[i] -= other[i];
        return *this;
    }

    Vector &operator*=(double d) {
        for (int i = 0; i < len; ++i)
            this->data[i] *= d;
        return *this;
    }

    Vector &operator/=(double d) {
        if (d == 0)
            throw std::invalid_argument("division by 0");
        for (int i = 0; i < len; ++i)
            this->data[i] /= d;
        return *this;
    }

    Vector &operator^=(double n) {
        for (int i = 0; i < this->len; ++i)
            this->data[i] = std::pow(this->data[i], n);
        return *this;
    }

    int length() const {
        return len;
    }

    double norm(double (* f) (const Vector &) = &euclidean_norm) const {
        return f(*this);
    }

    /*
        Normalize vector according to some norm `f`.
        The magnitude of this vector is thence 1,
        with respect to metric defined by `f`.
    */
    Vector &normalize(double (* f) (const Vector &) = &euclidean_norm) {
        auto magnitude = this->norm();
        for (int i = 0; i < this->length(); ++i)
            this->data[i] /= magnitude;
        return *this;
    }

    double mean() const {
        if (!this->len)
            return 0.0;
        double avg = 0.0;
        for (int i = 0; i < this->len; ++i)
            avg += this->data[i];
        return avg / this->len;
    }

    double var() const {
        if (!(this->len >> 1))
            return 0.0;
        double v = 0.0;
        double sum = 0.0;
        for (int i = 0; i < this->len; ++i) {
            sum += this->data[i];
            v += this->data[i] * this->data[i];
        }
        return (v - sum * sum / this->len) / (this->len - 1);
    }

    int max_index(double (* metric)(double) = &identity) const {
        int index = 0;
        for (int i = 1; i < this->len; ++i)
            if (metric(this->data[i]) > metric(this->data[index]))
                index = i;
        return index;
    }

    double max_element(double (* metric)(double)) const {
        if (this->len == 0)
            return 0.0;
        return this->data[this->max_index(metric)];
    }

    friend Vector operator*(double, const Vector &);
    friend std::ostream &operator<<(std::ostream &, const Vector &);
    friend std::ostream &operator<<(std::ostream &, Vector &&);
};

Vector one(int);

Vector elementary(int, int);

Vector random_vector(int);

std::ostream &operator<<(std::ostream &, const Vector &);
std::ostream &operator<<(std::ostream &, Vector &&);
Vector operator*(double, const Vector &);

Matrix eye(int);

struct Matrix {
    Vector *rows;
    Vector *cols;
    /*
        A matrix have row vector of dimension `col`,
        for total `row` of them.
    */
    int row, col;

    Matrix(int m, int n): rows(new Vector[m]), row(m), col(n) {
        for (int i = 0; i < row; ++i)
            rows[i] = Vector(col);
    }

    Matrix(int n): rows(new Vector[n]), row(n), col(n) {
        for (int i = 0; i < n; ++i)
            rows[i] = Vector(n);
    }

    Matrix(std::initializer_list<Vector> l): rows(new Vector[l.size()]), row(l.size()) {
        int j = 0;
        for (auto v: l)
            this->rows[j++] = v;
        this->col = rows[0].length();
    }

    Matrix(const Matrix &other): rows(new Vector[other.row]), row(other.row), col(other.col) {
        for (int i = 0;i < this->row; ++i)
            this->rows[i] = other[i];
    }

    ~Matrix() {
        delete[] rows;
    }

    Vector &operator[](int i) {
        return this->rows[i];
    }

    const Vector &operator[](int i) const {
        return this->rows[i];
    }

    Matrix &operator=(const Matrix &other) {
        if (this->row != other.row) {
            delete[] this->rows;
            this->row = other.row;
            this->rows = new Vector(this->row);
        }
        this->col = other.col;
        for (int i = 0; i < row; ++i)
            this->rows[i] = other[i];
        return *this;
    }

    Matrix operator+(const Matrix &other) const {
        if (this->row != other.row || this->col != other.col)
            throw std::length_error("incompatible dimension");
        Matrix r(*this);
        for (int i = 0; i < this->row; ++i)
            r.rows[i] += other.rows[i];
        return r;
    }

    Matrix operator-(const Matrix &other) const {
        if (this->row != other.row || this->col != other.col)
            throw std::length_error("incompatible dimension");
        Matrix r(*this);
        for (int i = 0; i < this->row; ++i)
            r.rows[i] -= other.rows[i];
        return r;
    }

    Matrix operator-() const {
        Matrix r(*this);
        for (int i = 0; i < this->row; ++i)
            r.rows[i] = -r.rows[i];
        return r;
    }

    /*
        Overloads 3 kinds of multiplication operations:
        1) multiply by scalar, which is a ring action
        2) multiply by vector, which is a linear space homomorphism
        3) multiply by matrix, which is a composition of linear space homomorphisms
    */

    Matrix operator*(double n) const {
        Matrix m(*this);
        for (int i = 0; i < row; ++i)
            m.rows[i] *= n;
        return m;
    }

    Vector operator*(const Vector &vec) const {
        /*
            Ordinary approach
        */
        if (this->col != vec.length())
            throw std::length_error("incompatible dimension");
        Vector v(this->row);
        for (int i = 0; i < col; ++i)
            v[i] = dot(this->rows[i], vec);
        return v;
    }

    Matrix operator*(const Matrix &other) const {
        /*
            Naive method, time complexity is O(n^3).
            View matrix-matrix multiplication as a change of basis,
            with tricks by treating multiplication of matrix as
            combinations of rows.
        */
        if (this->col != other.row)
            throw std::length_error("incompatible dimension");
        Matrix m(this->row, other.col);
        for (int i = 0; i < this->row; ++i)
            for (int j = 0; j < other.row; ++j)
                m[i] += other.rows[j] * this->rows[i][j];
        return m;
    }

    std::pair<Matrix, Matrix> gaussian() {
        /**
         * Performs gaussian elimination and returns callee's LU factorization.
         */
        if (this->col != this->row)
            throw std::invalid_argument("cannot perform LU factorize");
        Matrix l = eye(this->row);
        Matrix u = Matrix(*this);
        for (int i = 0; i < u.row - 1; ++i) {
            for (int j = i + 1; j < u.row; ++j) {
                if (u.rows[j][i] == 0)
                    continue;
                double c = u.rows[j][i] / u.rows[i][i];
                u.rows[j] -= (u.rows[i] * c);
                l[j][i] = c;
            }
        }
        return std::make_pair(l, u);
    }

    Vector solve(const Vector &b) {
        /**
         * Caution: This method is quite naive,
         * it does not check whether a matrix is invertible.
         */
        if (this->col != b.length())
            throw std::length_error("incompatible dimension");
        auto A = this->gaussian();
        Matrix L = A.first, U = A.second;
        /**
         * Solving for Lc = b.
         */
        Vector c(b);
        for (int i = 0; i < b.length(); ++i)
            for (int j = 0; j < i; ++j)
                c[i] -= L[i][j] * c[j];
        
        /**
         * Solving for Ux = c.
         */
        Vector x(c);
        for (int i = c.length() - 1; i >= 0; --i) {
            for (int j = c.length() - 1; j > i; --j)
                x[i] -= U[i][j] * x[j];
            x[i] /= U[i][i];
        }
        return x;
    }

    Vector diag() const {
        if (this->col != this->row)
            throw std::length_error("non-square matrix diagonal is ill-defined");
        Vector v(this->row);
        for (int i = 0; i < this->row; ++i)
            v[i] = this->rows[i][i];
        return v;
    }

    Matrix transpose() const {
        Matrix t(this->col, this->row);
        for (int i = 0; i < this->row; ++i)
            for (int j = 0; j < this->col; ++j)
                t[j][i] = this->rows[i][j];
        return t;
    }

    std::pair<int, int> max_index(double (* metric)(double) = &identity) const {
        int m = 0, n = 0;
        for (int i = 0; i < this->row; ++i) {
            int tmp = this->rows[i].max_index(metric);
            if (metric(this->rows[m][n]) < metric(this->rows[i][tmp])) {
                m = i;
                n = tmp;
            }
        }
        return std::make_pair(m, n);
    }

    /** 
     * Finds maximum off-diagonal element by `metric`
     * and return its index.
     */
    std::pair<int, int> pivot_index(double (* metric)(double) = &identity) const {
        if (this->row < 2 )
            return std::make_pair(-1, -1);
        int m = 0, n = 1;
        for (int i = 0; i < this->row; ++i) {
            int tmp = i == 0 ? 1 : 0;
            for (int j = 0; j < this->col; ++j) 
                if (i ^ j && metric(this->rows[i][j]) > metric(this->rows[i][tmp]))
                    tmp = j;
            if (metric(this->rows[m][n]) < metric(this->rows[i][tmp])) {
                m = i;
                n = tmp;
            }
        }
        return std::make_pair(m, n);
    }

    friend std::ostream &operator<<(std::ostream &, const Matrix &);
    friend std::ostream &operator<<(std::ostream &, Matrix &&);
    friend Matrix operator*(double, const Matrix &);
};

double rayleigh_quotient(const Matrix &, const Vector &);

Matrix eye(int);
Matrix givens(int, int, int, double);

std::ostream &operator<<(std::ostream &, const Matrix &);
std::ostream &operator<<(std::ostream &, Matrix &&);
Matrix operator*(double, const Matrix &);

bool is_zero(double n, double tol = 1e-14);
#endif
