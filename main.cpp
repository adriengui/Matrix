#include <iostream>
#include <list>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "Matrix.hpp"

using namespace std;

int main() {

	srand(time(nullptr));

	/* Matrix m1(5, 5, Random), m1com, m2(5, 5, Random), m3, I(5, 5, Identity), m5, m5i, m5m5i, m50, m6, m7, m8, m9, m10, m11, m12;
	float a, b;
	list<float> l;

	cout << "m1 =" << endl << m1 << endl;
	cout << "m2 =" << endl << m2 << endl;
	cout << "I =" << endl << I << endl;

	m2.exportToFile("m2.txt");

	m3 = Matrix(5, 5, "m2.txt");

	cout << "m3 = m2 =" << endl << m3 << endl;
	m3 = m3.transpose();
	cout << "m3 <= transpose(m3) =" << endl << m3 << endl;

	I *= 10;
	cout << "I <= I * 10 =" << endl << I << endl;

	m3 = m3 + I;
	cout << "m3 <= m3 + I =" << endl << m3 << endl;
	cout << "diagDominant(m3) ? " << boolalpha << m3.diagDominant() << endl << endl;

	m5 = m1 * m3;
	cout << "m5 = m1 * m3 =" << endl << m5 << endl;
	cout << "det(m5) = " << m5.det() << endl << endl;
	cout << "invertible(m5) ? " << m5.invertible() << endl << endl;

	m5i = m5.inv();
	cout << "inv(m5) =" << endl << m5i << endl;

	m5m5i = m5 * m5i;
	cout << "m5 * inv(m5) =" << endl << m5m5i << endl;

	I /= 10;
	cout << "I =" << endl << I << endl;
	cout << "m5 * inv(m5) = I ? " << (m5m5i == I) << endl << endl;

	m50 = m5 - m1;
	m50 = m50 / 2;
	m50(1, 1) = m50(1, 0) * m50(0, 1) + m50(0, 0);
	m50.removeColumn(3);
	m50.removeRow(4);

	cout << "m50 =" << endl << m50 << endl;
	cout << "min(m50) = " << m50.min() << endl;
	cout << "max(m50) = " << m50.max() << endl << endl;

	m6 = m50.getColumn(3);
	m7 = m50.getRow(2);
	cout << "column(m50)(3) =" << endl << m6 << endl;
	cout << "row(m50)(2) =" << endl << m7 << endl;

	a = (rand()%201)/100.0-1;
	b = sqrt(1-pow(a, 2));

	l = {
		a, b,
		-b, a
	};

	m8 = Matrix(2, 2, l);
	cout << "m8 =" << endl << m8 << endl;
	cout << "orthogonal(m8) ? " << m8.orthogonal() << endl << endl;

	m9 = m8.toColumn();
	m10 = m8.toRow();
	cout << "toColumn(m8) =" << endl << m9 << endl;
	cout << "toRow(m8) =" << endl << m10 << endl;

	m11 = m8;
	m11(0, 1) = 0;
	cout << "m11 =" << endl << m11 << endl;
	cout << "upperTriangular(m11) ? " << m11.upperTri() << endl;
	cout << "lowerTriangular(m11) ? " << m11.lowerTri() << endl;
	cout << "orthogonal(m11) ? " << m11.orthogonal() << endl;
	cout << "diagDominant(m11) ? " << m11.diagDominant() << endl << endl;

	cout << "m1 =" << endl << m1 << endl;

	m1com = m1.com();
	m12 = m1com.subMat(3, 3);
	cout << "com(m1) =" << endl << m1com << endl;
	cout << "minor(m1)(0, 1) = " << m1.Minor(0, 1) << endl << endl;
	cout << "subMat(com(m1))(3, 3) = " << endl << m12; */

	/* list<float> l = {
		12, -51, 4,
		6, 167, -68,
		-4, 24, -41
	};

	list<float> l2 = {
		-4, 20, 35, 5,
		-4, -30, -15, 55,
		-8, 40, -80, -65,
		23, -15, 30, 15
	};

	list<float> l3 = {
		0, 1, 2,
		4, 4.5, -7,
		-1, -1, 2
	}; */


	int M = 100, c = 5;//, t;//, n = 0, t;
	for(int i=0;i<c;i++) {
        Matrix A(M, M, Random);//, AdecompQR = A.decompQR();//, AdecompLU = A.decompLU();
        /* t = (A == Ainv);

        if(!t)
            n++; */
        // t = (A == AdecompLU);
        cout << (A == A.decompQR()) << endl;// " - " << t << endl;

        /* if(!t)
            cout << A << endl << endl << AdecompLU << endl << endl;
        */
        // cout << t << endl;
	}

	// cout << endl << (float) (n*100.0/c) << endl;
	// A.exportToFile("mat_A");
	// Ainv.exportToFile("mat_Ainv");
	// AAinv.exportToFile("mat_AAinv");

	// cout << A << endl << Ainv << endl << AAinv << endl;

	// A.invLU();

    // cout << AAinv << endl;
    // cout << A << endl;
    // A.exportToFile("testDetLU");
    // cout << A.det() << endl;
    // cout << A << endl << Ainv << endl << (A == Ainv) << endl;
	// cout << (AAinv == Matrix(M, M, Identity)) << endl << Ainv.getM() << endl;// << endl;// << endl << AAinv;
    // }
	/* Matrix A(8, 8, Random), Ainv, AAinv;

	for(int i=0;i<7;i++)
        for(int j=i+1;j<8;j++)
            A(i, j) = 0;

    Ainv = A.invLowerTri();
    AAinv = A * Ainv;

    cout << A << endl << Ainv << endl << AAinv; */

	return EXIT_SUCCESS;
}
