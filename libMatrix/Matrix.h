#pragma once
#include <mkl.h>

class Matrix {
    const double tol = 1.0e-10;
    // matrix(size, size)
    double *mat = NULL;
    int    size = 0;

    inline int idx(int i, int j);
    void fill(const double val);
    void copy(const Matrix &matrix);
    void gen_rand(const double val_min, const double val_max, const MKL_INT seed);
    int lu_decomp(const double *a, int *ipivot, double *lu);
    int inverse(double *a, double *a_inv);
    
public:
    Matrix();
    Matrix(const int size);
    Matrix(const Matrix &matrix);
    virtual ~Matrix();
    Matrix& operator=(const Matrix &rhs);
    Matrix& operator=(const double val);
    double operator[] (const int i);
    Matrix operator+(const Matrix &rhs);
    Matrix operator+(const double rhs);
    Matrix operator-(const Matrix &rhs);
    Matrix operator-(const double rhs);
    Matrix operator*(const Matrix &rhs);
    Matrix operator*(const double rhs);
    Matrix operator/(const Matrix &rhs);
    Matrix operator/(const double rhs);
    bool operator==(const Matrix &rhs);
    bool operator==(const double rhs);
    bool operator!=(const Matrix &rhs);
    bool operator!=(const double rhs);
    Matrix& operator+=(const Matrix &rhs);
    Matrix& operator+=(const double rhs);
    Matrix& operator-=(const Matrix &rhs);
    Matrix& operator-=(const double rhs);
    Matrix& operator*=(const Matrix &rhs);
    Matrix& operator*=(const double rhs);
    Matrix& operator/=(const Matrix &rhs);
    Matrix& operator/=(const double rhs);
    void set(const double val_min, const double val_max, const int seed);
    void set(const int i, const int j, const double val);
    void show();
    Matrix& transpose();
    Matrix& inverse();
};
