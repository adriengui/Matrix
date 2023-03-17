#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include "Matrix.hpp"

using namespace std;

Matrix::Matrix() : Matrix(1, 1, Zero) { }

Matrix::Matrix(int m, int n, Type t) {
	int i, j;

	if(m < 1)
		m = 1;

	if(n < 1)
		n = 1;

	M = m;
	N = n;

	mat = vec(M, vector<float> (N, 0));

	switch(t) {
		case Identity:
			if(M == N) {
				for(i = 0 ; i < M ; i++)
					(*this)(i, i) = 1;
			}

			break;

		case Random:
			for(i = 0 ; i < M ; i++)
				for(j = 0 ; j < N ; j++)
					(*this)(i, j) = pow(-1, (rand()&1)) * (rand()%1 + (rand()%100)/100.0);

			break;

		default:

			break;
	}
}

Matrix::Matrix(int m, int n, list<float> l) : Matrix(m, n, Zero) {
	int i = 0;

	if(l.size() != (unsigned int) M*N)
		return;

	for(float f : l) {
		(*this)(i/N, i%N) = f;

		i++;
	}
}

Matrix::Matrix(int m, int n, string s) : Matrix(m, n, Zero) {
	int i, j;
	string line;
	ifstream f(s);

	if(!f.is_open())
		return;

	for(i = 0 ; i < M ; i++) {
		getline(f, line);
		istringstream l(line);

		for(j = 0 ; j < N ; j++)
			l >> (*this)(i, j);
	}

	f.close();
}

Matrix Matrix::plus(Matrix m) {
	int i, j;
	Matrix s(M, N, Zero);

	if(M != m.getM() || N != m.getN())
		return Matrix(1, 1, Zero);

	for(i=0;i<M;i++)
		for(j = 0 ; j < N ; j++)
			s(i, j) = (*this)(i, j) + m(i, j);

	return s;
}

Matrix Matrix::minus(Matrix m) {
	int i, j;
	Matrix s(M, N, Zero);

	if(M != m.getM() || N != m.getN())
		return Matrix(1, 1, Zero);

	for(i = 0 ; i < M ; i++)
		for(j = 0 ; j < N ; j++)
			s(i, j) = (*this)(i, j) - m(i, j);

	return s;
}

Matrix Matrix::times(Matrix m) {
	int i, j, k;
	float x;
	Matrix s(M, m.getN(), Zero);

	if(N != m.getM())
		return Matrix(1, 1, Zero);

	for(i = 0 ; i < M ; i++) {
		for(j = 0 ; j < m.getN() ; j++) {
			x = 0;

			for(k = 0 ; k < N ; k++)
				x += (*this)(i, k) * m(k, j);

			s(i, j) = x;
		}
	}

	return s;
}

Matrix Matrix::operator+(const Matrix& m) {
	return(plus(m));
}

Matrix Matrix::operator-(const Matrix& m) {
	return(minus(m));
}

Matrix Matrix::operator*(const Matrix& m) {
	return(times(m));
}

string Matrix::toString() const {
	ostringstream s;
	string s1;
	const float n = LIM_TOSTRING;

	for(auto v : mat) {
		for(float f : v) {
			if(abs(f) < n)
				f = 0;
            s << setprecision(6) << f << "\t";
			// s << setfill(' ') << setw(12) << setprecision(6) << f << " ";
		}

		s << endl;
	}

	s1 = s.str();

	for(auto &i : s1)
        if(i == '.')
            i = ',';

	return s1;
}

void Matrix::print() {
	cout << toString();
}

ostream& operator<<(ostream& os, Matrix& m) {
	os << m.toString();

	return os;
}

void Matrix::exportToFile(string s) {
	ofstream f(s);

	if(!f.is_open())
		return;

	f << toString();

	f.close();
}

Matrix Matrix::getColumn(int n) {
	int i;
	Matrix v(M, 1, Zero);

	if(n >= N)
		return v;

	for(i = 0 ; i < M ; i++)
		v(i, 0) = (*this)(i, n);

	return v;
}

Matrix Matrix::getRow(int m) {
	Matrix v(1, N, Zero);

	if(m >= M)
		return v;

	v.mat[0] = mat[m];

	return v;
}

void Matrix::removeColumn(int n) {
	int i;

	if(n >= N || N == 1)
		return;

	for(i = 0 ; i < M ; i++)
		mat[i].erase(mat[i].begin()+n);

	N--;
}

void Matrix::removeRow(int m) {
	if(m >= M || M == 1)
		return;

	mat.erase(mat.begin()+m);

	M--;
}

float Matrix::max() {
	int i, j;
	float n = (*this)(0, 0);

	for(i = 0 ; i < M ; i++)
		for(j = 0 ; j < N ; j++)
			if((*this)(i, j) > n)
				n = (*this)(i, j);

	return n;
}

float Matrix::min() {
	int i, j;
	float n = (*this)(0, 0);

	for(i = 0 ; i < M ; i++)
		for(j = 0 ; j < N ; j++)
			if((*this)(i, j) < n)
				n = (*this)(i, j);

	return n;
}

bool Matrix::invertible() {
	const float n = LIM;

	return !(abs(det()) < n);
}

bool Matrix::orthogonal() {
	Matrix product = *this * transpose();

	if(!square())
		return 0;

	return (product == Matrix(M, M, Identity));
}

bool Matrix::diagDominant() {
	int i, j;
	float n;

	if(!square())
        return 0;

	for(i = 0 ; i < M ; i++) {
		n = 0;

		for(j = 0 ; j < N ; j++)
			if(i != j)
				n += abs((*this)(i, j));

		if(abs((*this)(i, i)) < n)
			return 0;
	}

	return 1;
}

bool Matrix::upperTri() {
	int i, j;
	const float n = LIM;

	if(!square())
        return 0;

	for(i = 1 ; i < M ; i++)
		for(j = 0 ; j < i ; j++)
			if(!(abs((*this)(i, j)) < n))
				return 0;

	return 1;
}

bool Matrix::lowerTri() {
	int i, j;
	const float n = LIM;

	if(!square())
        return 0;

	for(i = 0 ; i < M-1 ; i++)
        for(j = i+1 ; j < N ; j++)
			if(!(abs((*this)(i, j)) < n))
				return 0;

	return 1;
}

bool Matrix::square() {
    if(M != N)
		return 0;

	if(M == 1)
		return 0;

    return 1;
}

Matrix Matrix::toColumn() {
	return toRow().transpose();
}

Matrix Matrix::toRow() {
	int i, j;
	Matrix s(1, M*N, Zero);

	for(i = 0 ; i < M ; i++)
		for(j = 0 ; j < N ; j++)
			s(0, i*N+j) = (*this)(i, j);

	return s;
}

Matrix Matrix::transpose() {
	int i, j;
	Matrix s(N, M, Zero);

	for(i = 0 ; i < N ; i++)
		for(j = 0 ; j < M ; j++)
			s(i, j) = (*this)(j, i);

	return s;
}

float Matrix::det() {
	int i;
	float d = 0;
	const float n = LIM;

	if(!square())
		return 0;

	if(M == 2)
		return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);

	for(i = 0 ; i < M ; i++)
		if(!(abs((*this)(i, 0)) < n))
			d += pow(-1, i) * (*this)(i, 0) * Minor(i, 0);

	return d;
}

float Matrix::Minor(int i, int j) {
	return subMat(i, j).det();
}

Matrix Matrix::operator*(float n) {
	int i, j;
	Matrix s(M, N, Zero);

	for(i = 0 ; i < M ; i++)
		for(j = 0 ; j < N ; j++)
			 s(i, j) = (*this)(i, j) * n;

	return s;
}

Matrix Matrix::operator/(float n) {
	return (*this) * (1/n);
}

Matrix Matrix::operator*=(float n) {
	*this = *this * n;

	return *this;
}

Matrix Matrix::operator/=(float n) {
	*this = *this / n;

	return *this;
}

Matrix Matrix::inv() {
	// if(!invertible())
		// return Matrix(1, 1, Zero);
    float d = det();
    const float n = LIM;

    if(abs(d) < n)
        return Matrix(1, 1, Zero);

    return com().transpose() / d;
}

Matrix Matrix::invUpperTri() {
    int i, k, n;
    float sum;
    const float e = LIM;
	Matrix s(M, M, Zero);

    if(!square())
		return Matrix(1, 1, Zero);

    for(i = 0 ; i < M ; i++) {
        if(abs((*this)(i, i)) < e)
			return Matrix(1, 1, Zero);

		s(i, i) = 1 / (*this)(i, i);
    }

	for(n = 1 ; n < M ; n++) {
		for(i = 0 ; i < M-n ; i++) {
			sum = 0;

			for(k = 1 ; k <= n ; k++)
				sum -= (*this)(i, i+k) * s(i+k, i+n);

            sum /= (*this)(i, i);

			s(i, i+n) = sum;
		}
	}

	return s;
}

Matrix Matrix::invLowerTri() {
    int i, k, n;
    float sum;
    const float e = LIM;
	Matrix s(M, M, Zero);

    if(!square())
		return Matrix(1, 1, Zero);

    for(i = 0 ; i < M ; i++) {
        if(abs((*this)(i, i)) < e)
			return Matrix(1, 1, Zero);

		s(i, i) = 1 / (*this)(i, i);
    }

	for(n = 1 ; n < M ; n++) {
		for(i = n ; i < M ; i++) {
			sum = 0;

			for(k = 1 ; k <= n ; k++)
				sum -= (*this)(i, i-k) * s(i-k, i-n);

            sum /= (*this)(i, i);

			s(i, i-n) = sum;
		}
	}

	return s;
}

void Matrix::QR(Matrix &Q, Matrix &R) {
    int i, j, k;
	Matrix A = *this, QA = *this, Qp, u, v, vvt;
	vector<Matrix> Qv;

	Qv = vector<Matrix> (M-1, Matrix(M, M, Zero));

	for(i = 1 ; i < M-1 ; i++)
		for(j = 0 ; j < i ; j++)
			Qv[i](j, j) = 1;

	for(i = 0 ; i < M-1 ; i++) {
		u = A.getColumn(0);
		u(0, 0) -= u.norm();

		v = u / u.norm();
		vvt = v * v.transpose();

		Qp = Matrix(M-i, M-i, Identity);
		Qp = Qp - vvt * 2;

		for(j = i ; j < M ; j++)
			for(k = i ; k < M ; k++)
				Qv[i](j, k) = Qp(j-i, k-i);

		QA = Qv[i] * QA;

		A = QA;

		for(j = 0 ; j <= i ; j++) {
			A.removeColumn(0);
			A.removeRow(0);
		}
	}

	Q = Qv[0];

	for(i = 1 ; i < M-1 ; i++)
		Q = Q * Qv[i];

	R = QA;
}

Matrix Matrix::invQR() {
    Matrix Q, R;

    if(!square())
        return Matrix(1, 1, Zero);

    QR(Q, R);

    return R.invUpperTri() * Q.transpose();
}

Matrix Matrix::decompQR() {
    Matrix Q, R;

    if(!square())
        return Matrix(1, 1, Zero);

    QR(Q, R);

    return Q * R;
}

void Matrix::LU(Matrix &L, Matrix &U) {
    int i, j;
    const float n = LIM;
    Matrix A = *this, Linv;
    vector<Matrix> Lv(M-1, Matrix(M, M, Identity));

    for(i = 0 ; i < M-1 ; i++) {
        for(j = i+1 ; j < M ; j++) {
            if(abs(A(i, i)) < n) {
                L = Matrix(1, 1, Zero);
                U = Matrix(1, 1, Zero);

                return;
            }

            Lv[i](j, i) = - A(j, i) / A(i, i);
        }

        A = Lv[i] * A;
    }

    Linv = Lv[M-2];

    for(i = M-3 ; i >= 0 ; i--)
        Linv = Linv * Lv[i];

    U = A;
    L = Linv.invLowerTri();
}

Matrix Matrix::invLU() {
    Matrix L, U;

    if(!square())
        return Matrix(1, 1, Zero);

    LU(L, U);

    return U.invUpperTri() * L.invLowerTri();
}

Matrix Matrix::decompLU() {
    Matrix L, U;

    if(!square())
        return Matrix(1, 1, Zero);

    LU(L, U);

    return L * U;
}

double Matrix::detLU() {
    int i, j;
    const float n = LIM;
    Matrix A = *this, L;

    if(!square())
        return 0;

    for(i = 0 ; i < M-1 ; i++) {
        L = Matrix(M, M, Identity);

        for(j = i+1 ; j < M ; j++) {
            if(abs(A(i, i)) < n)
                return 0;

            L(j, i) = - A(j, i) / A(i, i);
        }

        A = L * A;
    }

    return A.prodDiag();
}

double Matrix::prodDiag() {
    int i;
    double f = (double) (*this)(0, 0);

    if(!square())
        return 0;

    for(i = 1 ; i < M; i++)
        f *= (double) (*this)(i, i);

    return f;
}

float Matrix::norm() {
	int i;
	float n = 0;

	if(N != 1)
		return 0;

	for(i = 0 ; i < M ; i++)
		n += pow((*this)(i, 0), 2);

	return sqrt(n);
}

Matrix Matrix::subMat(int l, int c) {
	int i, j, x=0, y;
	Matrix s(M-1, N-1, Zero);

	if(l >= M || c >= N || M == 1 || N == 1)
		return Matrix(1, 1, Zero);

	for(i = 0 ; i < M ; i++) {
		if(i != l) {
			y=0;

			for(j = 0 ; j < N ; j++)
				if(j != c)
					s(x, y++) = (*this)(i, j);

			x++;
		}
	}

	return s;
}

Matrix Matrix::com() {
	int i, j;
	Matrix s(M, N, Zero);

	if(M != N)
		return Matrix(1, 1, Zero);

	for(i = 0 ; i < M ; i++)
		for(j = 0 ; j < N ; j++)
			s(i, j) = pow(-1, i+j) * Minor(i, j);

	return s;
}

bool Matrix::operator==(Matrix m) {
	int i, j;
	const float n = LIM_EQUAL;

	if(M != m.getM() || N != m.getN())
		return 0;

	for(i = 0 ; i < M ; i++)
		for(j = 0 ; j < N ; j++)
			if(!(abs((*this)(i, j) - m(i, j)) < n))
				return 0;

	return 1;
}

float& Matrix::operator()(int i, int j) {
	return mat[i][j];
}

int Matrix::getM() const {
	return M;
}

int Matrix::getN() const {
	return N;
}
