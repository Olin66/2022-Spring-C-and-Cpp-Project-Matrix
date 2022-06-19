#include "matrix.h"
#include <iostream>
#include <complex>
#define BASIC_MATRIX_INT BasicMatrix<int>
#define SPARSE_MATRIX_INT SparseMatrix<int>
#define BASIC_MATRIX_COMPLEX BasicMatrix<complex<int>>
#define SPARSE_MATRIX_COMPLEX SparseMatrix<complex<int>>
using namespace std;
using namespace mat;

void test() {
    cout << "test successfully!" << endl;
}

void bmat_multi_test1() {
    cout << "Basic matrix<int> multiplication test1 begin:" << endl;
    BASIC_MATRIX_INT *bm1 = new BASIC_MATRIX_INT(4, 4, 2);
    cout << "Matrix 1:" << endl;
    bm1->show();
    Matrix<int> * m = Matrix<int>::eye(4, 4, BASIC_MATRIX);
    cout << "Matrix 2:" << endl;
    m->show();
    ((BasicMatrix<int>*)m)->dotProduct(*bm1);
    cout << "Matrix 2 dot product matrix 1(new matrix 2):" << endl;
    m->show();
    cout << "Matrix 1 cross product matrix 2(new matrix 1):" << endl;
    bm1->crossProduct(*((BASIC_MATRIX_INT*)m));
    bm1->show();
    delete bm1;
    delete m;
    cout << "test end" << endl << endl;
}

void bmat_multi_test2() {
    cout << "Basic matrix<int> multiplication test2 begin:" << endl;
    int *data1 = new int[6]{1,2,3,4,5,6};
    int *data2 = new int[12]{1,2,3,4,5,6,6,5,4,3,2,1};
    BASIC_MATRIX_INT *bm1 = new BASIC_MATRIX_INT(2, 3, data1);
    BASIC_MATRIX_INT *bm2 = new BASIC_MATRIX_INT(3, 4, data2);
    cout << "Matrix 1:" << endl;
    bm1->show();
    cout << "Matrix 2:" << endl;
    bm2->show();
    cout << "Matrix 1 cross product matrix 2(new matrix 1):" << endl;
    bm1->crossProduct(*bm2);
    bm1->show();
    delete bm1;
    delete bm2;
    cout << "test end" << endl << endl;
}

void bmat_multi_test3() {
    cout << "Basic matrix<complex> multiplication test3 begin:" << endl;
    complex<int> *c_data1 = new complex<int>[2]{
        complex<int>(0, 1),
        complex<int>(0, -1)
    };
    BASIC_MATRIX_COMPLEX *bm1 = new BASIC_MATRIX_COMPLEX(2, 1, c_data1);
    complex<int> *c_data2 = new complex<int>[4]{
        complex<int>(1, 1),
        complex<int>(1, -1),
        complex<int>(2, 1),
        complex<int>(5, -1)
    };
    BASIC_MATRIX_COMPLEX *bm2 = new BASIC_MATRIX_COMPLEX(2, 2, c_data2);
    cout << "Matrix 1:" << endl;
    bm1->show();
    cout << "Matrix 2:" << endl;
    bm2->show();
    cout << "Matrix 1 dot product matrix 2(new matrix 1):" << endl;
    bm1->dotProduct(*bm2);
    bm1->show();
    delete bm1;
    delete bm2;
    cout << "test end" << endl << endl;
}

void smat_multi_test1() {
    // cout << "Sparse matrix<int> multiplication test1 begin:" << endl;
    // SPARSE_MATRIX_INT *bm1 = new SPARSE_MATRIX_INT(4, 4, 2);
    // cout << "Matrix 1:" << endl;
    // bm1->show();
    // Matrix<int> * m = Matrix<int>::eye(4, 4, BASIC_MATRIX);
    // cout << "Matrix 2:" << endl;
    // m->show();
    // ((BasicMatrix<int>*)m)->dotProduct(*bm1);
    // cout << "Matrix 2 dot product matrix 1(new matrix 2):" << endl;
    // m->show();
    // cout << "Matrix 1 cross product matrix 2(new matrix 1):" << endl;
    // bm1->crossProduct(*((SPARSE_MATRIX_INT*)m));
    // bm1->show();
    // delete bm1;
    // delete m;
    // cout << "test end" << endl << endl;
}

void smat_multi_test2() {
    cout << "Sparse matrix<int> multiplication test2 begin:" << endl;
    int *data1 = new int[6]{1,2,3,4,5,6};
    int *data2 = new int[12]{1,2,3,4,5,6,6,5,4,3,2,1};
    SPARSE_MATRIX_INT *bm1 = new SPARSE_MATRIX_INT(2, 3, data1);
    SPARSE_MATRIX_INT *bm2 = new SPARSE_MATRIX_INT(3, 4, data2);
    cout << "Matrix 1:" << endl;
    bm1->show();
    cout << "Matrix 2:" << endl;
    bm2->show();
    cout << "Matrix 1 cross product matrix 2(new matrix 1):" << endl;
    bm1->crossProduct(*bm2);
    bm1->show();
    delete bm1;
    delete bm2;
    cout << "test end" << endl << endl;
}

void smat_multi_test3() {
    cout << "Sparse matrix<complex> multiplication test3 begin:" << endl;
    complex<int> *c_data1 = new complex<int>[2]{
        complex<int>(0, 1),
        complex<int>(0, -1)
    };
    SPARSE_MATRIX_COMPLEX *bm1 = new SPARSE_MATRIX_COMPLEX(2, 1, c_data1);
    complex<int> *c_data2 = new complex<int>[4]{
        complex<int>(1, 1),
        complex<int>(1, -1),
        complex<int>(2, 1),
        complex<int>(5, -1)
    };
    SPARSE_MATRIX_COMPLEX *bm2 = new SPARSE_MATRIX_COMPLEX(2, 2, c_data2);
    cout << "Matrix 1:" << endl;
    bm1->show();
    cout << "Matrix 2:" << endl;
    bm2->show();
    cout << "Matrix 1 dot product matrix 2(new matrix 1):" << endl;
    bm1->dotProduct(*bm2);
    bm1->show();
    delete bm1;
    delete bm2;
    cout << "test end" << endl << endl;
}