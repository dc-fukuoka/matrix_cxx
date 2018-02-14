#include <iostream>
#include <iomanip>
#include <cstddef>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <mkl.h>

class Matrix {
    const double tol = 1.0e-10;
    // matrix(size, size)
    double *mat = NULL;
    int    size = 0;

    inline int idx(int i, int j) {return size*i + j;} // row major
    
    void fill(const double val) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (auto i=0; i<size*size; i++)
	    mat[i] = val;
    }
    
    void copy(const Matrix &matrix) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (auto i=0; i<size*size; i++) {
	    mat[i] = matrix.mat[i];
	}
    }
    
    void gen_rand(const double val_min, const double val_max, const MKL_INT seed) {
	VSLStreamStatePtr stream;
	int stat;
	MKL_INT brng = VSL_BRNG_MT19937;
	int method   = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE;
	
	stat = vslNewStream(&stream, brng, seed);
	assert(stat == VSL_ERROR_OK);
	stat = vdRngUniform(method, stream, size*size, mat, val_min, val_max);
	assert(stat == VSL_ERROR_OK);
	stat = vslDeleteStream(&stream);
	assert(stat == VSL_ERROR_OK);
    }
    
    int lu_decomp(const double *a, int *ipivot, double *lu) {
	int i, j, k;
	int ip, tmp_ip;
	double tmp, max0, w;
#ifdef _OPENMP	
#pragma omp parallel private(i)
	{
#pragma omp for
#endif
	    for (i=0; i<size*size; i++)
		lu[i] = a[i];
#ifdef _OPENMP
#pragma omp for
#endif
	    for (i=0; i<size; i++)
		ipivot[i] = i;
#ifdef _OPENMP
	} // omp parallel
#endif

	for (k=0; k<size-1; k++) {
	    max0 = std::abs(lu[idx(k, k)]);
	    ip   = k;
	    for (i=k+1; i<size; i++) {
		tmp = std::abs(lu[idx(i, k)]);
		if (tmp > max0) {
		    max0 = tmp;
		    ip   = i;
		}
	    }
	    if (max0 <= tol) {
		std::cout << "one of diagonal component is smaller than " << tol << '\n';
		return -1;
	    }
	    
	    if (ip != k) {
#ifdef _OPENMP
#pragma omp parallel for private(tmp)
#endif
		for (j=k; j<size; j++) {
		    tmp            = lu[idx(ip, j)];
		    lu[idx(ip, j)] = lu[idx(k , j)];
		    lu[idx(k , j)] = tmp;
		}
		tmp_ip     = ipivot[ip];
		ipivot[ip] = ipivot[k];
		ipivot[k]  = tmp_ip;
#ifdef _OPENMP
#pragma omp parallel for private(tmp)
#endif
		for (j=0; j<k; j++) {
		    tmp            = lu[idx(k,  j)];
		    lu[idx(k,  j)] = lu[idx(ip, j)];
		    lu[idx(ip, j)] = tmp;
		}
	    }
#ifdef _OPENMP
#pragma omp parallel for private(w)
#endif
	    for (i=k+1; i<size; i++) {
		w             = lu[idx(i, k)]/lu[idx(k, k)];
		lu[idx(i, k)] = w;
		for (j=k+1; j<size; j++)
		    lu[idx(i, j)] = lu[idx(i, j)] - w*lu[idx(k, j)];
	    }
	}
	return 0;
    }
    
    int inverse(double *a, double *a_inv) {
	double lu[size*size];
	int ipivot[size];
	int i, j, k;
	double unit_vec[size], y[size];
	double tmp;
	
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (i=0; i<size; i++) {
	    ipivot[i] = 0;
	    for (j=0; j<size; j++) {
		lu[idx(i, j)]     = 0.0;
		a_inv[idx(i, j)]  = 0.0;
	    }
	}
	int ret = lu_decomp(a, ipivot, lu);
	if (ret) {
	    std::cout << "LU decomposition failed.\n";
	    return ret;
	}
	for (k=0; k<size; k++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
	    for (i=0; i<size; i++)
		unit_vec[i] = 0.0;
	    unit_vec[k] = 1.0;
	    // forward substitution
	    y[0] = unit_vec[ipivot[0]];
	    for (i=1; i<size; i++) {
		tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
		for (j=0; j<i; j++)
		    tmp = tmp + lu[idx(i, j)]*y[j];
		y[i] = unit_vec[ipivot[i]] - tmp;
	    }
	    
	    // backward substitution
	    a_inv[idx(size-1, k)] = y[size-1]/lu[idx(size-1, size-1)];
	    for (i=size-2; i>=0; i--) {
		tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
		for (j=i+1; j<size; j++)
		    tmp = tmp + lu[idx(i, j)]*a_inv[idx(j, k)];
		a_inv[idx(i, k)] = (y[i] - tmp)/lu[idx(i, i)];
	    }
	}
	return 0;
    }
    
public:
    Matrix() : size(1) {
	mat = new double[size*size];
	fill(0.0);
    }
    
    Matrix(const int size) : size(size) {
	mat = new double[size*size];
	fill(0.0);
    }
    Matrix(const Matrix &matrix) : size(matrix.size) {
	mat = new double[size*size];
	copy(matrix);
    }
    
    virtual ~Matrix() {
	if (mat) delete[] mat;
    }
    
    Matrix& operator=(const Matrix &rhs) {
	if (this != &rhs) {
	    if (size != rhs.size) {
		delete[] mat;
		size = rhs.size;
		mat = new double[size*size];
	    }
	    copy(rhs);
	}
	return *this;
    }

    Matrix& operator=(const double val) {
	fill(val);
	return *this;
    }
    
    double operator[] (const int i) {
	return mat[i];
    }
    
    Matrix operator+(const Matrix &rhs) {
	Matrix matrix(*this);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (auto i=0; i<size*size; i++)
	    matrix.mat[i] += rhs.mat[i];
	return matrix;
    }

    Matrix operator+(const double rhs) {
	Matrix matrix(*this);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (auto i=0; i<size*size; i++)
	    matrix.mat[i] += rhs;
	return matrix;
    }
    
    Matrix operator-(const Matrix &rhs) {
	Matrix matrix(*this);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (auto i=0; i<size*size; i++)
	    matrix.mat[i] -= rhs.mat[i];
	return matrix;
    }

    Matrix operator-(const double rhs) {
	Matrix matrix(*this);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (auto i=0; i<size*size; i++)
	    matrix.mat[i] -= rhs;
	return matrix;
    }
    
    Matrix operator*(const Matrix &rhs) {
	int i, j, k;
	Matrix matrix(size);
#ifdef _OPENMP
#pragma omp parallel for private(i, j, k)
#endif
	for (i=0; i<size; i++) {
	    for (k=0; k<size; k++) {
		for (j=0; j<size; j++) {
		    matrix.mat[idx(i, j)] += mat[idx(i, k)]*rhs.mat[idx(k, j)];
		}
	    }
	}
	return matrix;
    }

    Matrix operator*(const double rhs) {
	Matrix matrix(*this);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (auto i=0; i<size*size; i++)
	    matrix.mat[i] *= rhs;
	return matrix;
    }

    Matrix operator/(const Matrix &rhs) {
	Matrix rhs_inv(size), matrix(*this);
	auto ret = inverse(rhs.mat, rhs_inv.mat);
	if (ret) {
	    std::cout << "LU decomposition or inverse failed.\n";
	    return Matrix(0);
	}
	matrix *= rhs_inv;
	return matrix;
    }

    Matrix operator/(const double rhs) {
	Matrix matrix(*this);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (auto i=0; i<size*size; i++)
	    matrix.mat[i] /= rhs;
	return matrix;
    }
    
    bool operator==(const Matrix &rhs) {
	double diff, max_err;
	max_err = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(diff) reduction(max:max_err)
#endif
	for (auto i=0; i<size*size; i++) {
	    diff    = std::abs(mat[i]-rhs.mat[i]);
	    max_err = std::max(max_err, diff);
	}
	if (max_err >= tol) {
	    std::cout << "maximum error is: " << max_err << '\n';
	    return false;
	}
	return true;
    }

    bool operator==(const double rhs) {
	double diff, max_err;
	max_err = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(diff) reduction(max:max_err)
#endif
	for (auto i=0; i<size*size; i++) {
	    diff    = std::abs(mat[i]-rhs);
	    max_err = std::max(max_err, diff);
	}
	if (max_err >= tol) {
	    std::cout << "maximum error is: " << max_err << '\n';
	    return false;
	}
	return true;
    }
    
    bool operator!=(const Matrix &rhs) {
	if (*this == rhs) {
	    return false;
	} else {
	    return true;
	}
    }

    bool operator!=(const double rhs) {
	if (*this == rhs) {
	    return false;
	} else {
	    return true;
	}
    }

    Matrix& operator+=(const Matrix &rhs) {
	*this = *this + rhs;
	return *this;
    }

    Matrix& operator+=(const double rhs) {
	*this = *this + rhs;
	return *this;
    }

    Matrix& operator-=(const Matrix &rhs) {
	*this = *this - rhs;
	return *this;
    }

    Matrix& operator-=(const double rhs) {
	*this = *this - rhs;
	return *this;
    }

    Matrix& operator*=(const Matrix &rhs) {
	*this = *this * rhs;
	return *this;
    }

    Matrix& operator*=(const double rhs) {
	*this = *this * rhs;
	return *this;
    }

    Matrix& operator/=(const Matrix &rhs) {
	*this = *this / rhs;
	return *this;
    }

    Matrix& operator/=(const double rhs) {
	*this = *this / rhs;
	return *this;
    }
    
    void setMatrix(const double val_min, const double val_max, const int seed) {
	gen_rand(val_min, val_max, seed);
    }
    
    void showMatrix() {
	for (auto i=0; i<size; i++) {
	    for (auto j=0; j<size; j++) {
		std::cout << std::setw(14) << std::setprecision(5) << std::scientific <<mat[idx(i, j)] << " ";
	    }
	    std::cout << std::endl;
	}
    }

    Matrix& transpose() {
	// not good...
	Matrix tmp = *this;
	int i, j;
#ifdef _OPENMP
#pragma omp parallel for private(i, j)
#endif
	for (i=0; i<size; i++) {
	    for (j=0; j<size; j++) {
		if (i != j)
		    mat[idx(i, j)] = tmp.mat[idx(j, i)];
	    }
	}
	return *this;
    }
    
    Matrix& inverse() {
	Matrix inv = *this;
	auto ret = inverse(mat, inv.mat);
	if (ret) {
	    std::cout << __func__ << " failed.\n";
	    fill(-1.0);
	}
	return *this;
    }
    void set_at(const int i, const int j, const double val) {
	mat[idx(i, j)] = val;
    }
};

int main(int argc, char **argv) {
    int size = (argc == 1) ? 4 : std::atoi(argv[1]);
    Matrix A(size), B(size), C(size);
    A.setMatrix(-1.0, 1.0, 5555);
    B.setMatrix(-2.0, 2.0, 7777);
    C = 0.0;
    // to have the same order as fortran
    A.transpose();
    B.transpose();

    // test for pivoting
    for (auto i=0; i<size; i++)
	B.set_at(i, i, 0.0);
    
    std::cout << "A:\n";
    A.showMatrix();
    std::cout << "B:\n";
    B.showMatrix();
    C += A;
    C -= A;
    std::cout << "C: \n";
    C.showMatrix();
    C = A*B;
    std::cout << "C = A*B\n";
    C.showMatrix();

    C /= B;
    std::cout << "C /= B(= A)\n";
    C.showMatrix();
    if (C == A)
	std::cout << "C == A\n";
    if (A != B)
	std::cout << "A != B\n";

    Matrix Ainv = A.inverse();
    std::cout << "Ainv = A^-1\n";
    Ainv.showMatrix();

    C = C - C;
    std::cout << "trivial test: C = C - C (=0)\n";
    C.showMatrix();

    if (C == 0.0)
	std::cout << "C == 0.0\n";

    Matrix D = 1.0;
    if (C != 2.0)
	std::cout << "C != 2.0\n";
    
    return 0;
}
