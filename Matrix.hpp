#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <string>

using namespace std;

enum Type {Identity, Zero, Random};
typedef vector<vector<float>> vec;

#define LIM 1e-5
#define LIM_TOSTRING 1e-6
#define LIM_EQUAL 1e-5

class Matrix {
	private:
		int M;
		int N;
		vec mat;

	public:
		Matrix();
		Matrix(int m, int n, Type t);
		Matrix(int m, int n, list<float> l);
		Matrix(int m, int n, string s);

		Matrix plus(Matrix m);
		Matrix minus(Matrix m);
		Matrix times(Matrix m);

		Matrix operator+(const Matrix& m);
		Matrix operator-(const Matrix& m);
		Matrix operator*(const Matrix& m);

		string toString() const;
		void print();
		friend ostream& operator<<(ostream& os, Matrix& m);
		void exportToFile(string s);

		Matrix getColumn(int n);
		Matrix getRow(int m);
		void removeColumn(int n);
		void removeRow(int m);
		float max();
		float min();

		bool invertible();
		bool orthogonal();
		bool diagDominant();
		bool upperTri();
		bool lowerTri();
		bool square();

		Matrix toColumn();
		Matrix toRow();

		Matrix transpose();
		float det();
		float Minor(int i, int j);
		Matrix operator*(float n);
		Matrix operator/(float n);
		Matrix operator*=(float n);
		Matrix operator/=(float n);
		Matrix inv();
		Matrix invUpperTri();
		Matrix invLowerTri();
		void QR(Matrix &Q, Matrix &R);
		Matrix invQR();
		Matrix decompQR();
		void LU(Matrix &L, Matrix &U);
		Matrix invLU();
		Matrix decompLU();
		double detLU();
		double prodDiag();
		float norm();

		Matrix subMat(int l, int c);
		Matrix com();
		bool operator==(Matrix m);

		float& operator()(int i, int j);
		int getM() const;
		int getN() const;

};
