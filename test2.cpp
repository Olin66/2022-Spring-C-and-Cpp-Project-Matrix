#include "matrix.h"
#include <iostream>
#include <complex>
#define BASIC_MATRIX_INT BasicMatrix<int>
#define SPARSE_MATRIX_INT SparseMatrix<int>
#define BASIC_MATRIX_COMPLEX BasicMatrix<complex<int>>
#define SPARSE_MATRIX_COMPLEX SparseMatrix<complex<int>>
using namespace std;
using namespace mat;

void bmat_multi_test1() {
    cout << "Basic matrix<int> multiplication test1 begin:" << endl;
    BASIC_MATRIX_INT *bm1 = new BASIC_MATRIX_INT(4, 1, 2);
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
    cout << "Matrix 2:" << endl;
    bm2->show();
    cout << "Matrix 1 cross product matrix 2(new matrix 1:"  << endl;
    bm1->crossProduct(*bm2);
    bm1->show();
    delete bm1;
    delete bm2;
    cout << "test end" << endl << endl;
}

void bmat_multi_test4() {
    cout << "Basic matrix<int> multiplication test4 begin:" << endl;
    BASIC_MATRIX_INT *bm1 = new BASIC_MATRIX_INT(4, 1, 2);
    cout << "Matrix 1:" << endl;
    bm1->show();
    BASIC_MATRIX_INT *bm2 = (BASIC_MATRIX_INT*)Matrix<int>::eye(4, 4, BASIC_MATRIX);
    cout << "Matrix 2:" << endl;
    bm2->show();
    bm1->dotProduct(*bm2);
    cout << "Matrix 1 dot product matrix 2(new matrix 1):" << endl;
    bm1->show();
    cout << "Matrix 1 cross product matrix 1(new matrix 1):" << endl;
    bm1->crossProduct(*bm1);
    bm1->show();
    delete bm1;
    delete bm2;
    cout << "test end" << endl << endl;
}

void smat_multi_test1() {
    cout << "Sparse matrix<int> multiplication test1 begin:" << endl;
    vector<Triple<int>> tri1;
    tri1.push_back(Triple<int>(1, 2, 1));
    tri1.push_back(Triple<int>(4, 4, 1));
    vector<Triple<int>> tri2;
    tri2.push_back(Triple<int>(1, 1, 1));
    tri2.push_back(Triple<int>(2, 1, 1));
    tri2.push_back(Triple<int>(3, 1, 1));
    SPARSE_MATRIX_INT *bm1 = new SPARSE_MATRIX_INT(5, 5, tri1);
    cout << "Matrix 1:" << endl;
    bm1->show();
    SPARSE_MATRIX_INT *bm2 = new SPARSE_MATRIX_INT(5, 2, tri2);
    cout << "Matrix 2:" << endl;
    bm2->show();
    // cout << "Matrix 2:" << endl;
    // bm2->show();
    // ((BasicMatrix<int>*)m)->dotProduct(*bm1);
    // cout << "Matrix 2 dot product matrix 1(new matrix 2):" << endl;
    // bm2->show();
    cout << "Matrix 1 cross product matrix 2(new matrix 1):" << endl;
    bm1->crossProduct(*bm2);
    bm1->show();
    delete bm1;
    delete bm2;
    cout << "test end" << endl << endl;
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
    cout << "Matrix 2:" << endl;
    bm2->show();
    cout << "Matrix 1 cross product matrix 2(new matrix 1:"  << endl;
    bm1->crossProduct(*bm2);
    bm1->show();
    delete bm1;
    delete bm2;
    cout << "test end" << endl << endl;
}

void bmat_conv_test1() {
    cout << "Basic matrix convolution test begin:" << endl;
    const int TEST_SIZE = 16;
    int *_data1 = new int[TEST_SIZE]{
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16};
    BASIC_MATRIX_INT bm1(4, 4, _data1);

    int *_data2 = new int[9]{
        -1, -2, -1,
        0, 0, 0,
        1, 2, 1};
    BASIC_MATRIX_INT bm2(3, 3, _data2);
    bm1.show();
    bm2.show();
    BASIC_MATRIX_INT* mat_ans = bm1.convolve(bm2, 1, 1);
    cout << " Answer of matrix convolution" << endl;
    mat_ans->show();
    cout << "Test end." << endl << endl; 
}

void smat_conv_test1() {
    cout << "Sparse matrix convolution begin:" << endl;
    const int TEST_SIZE = 16;
    vector<Triple<int>> tri1;
    tri1.push_back(Triple<int>(1, 2, 1));
    tri1.push_back(Triple<int>(4, 4, 1));
    vector<Triple<int>> tri2;
    tri2.push_back(Triple<int>(0, 0, 1));
    tri2.push_back(Triple<int>(1, 2, -1));
    SPARSE_MATRIX_INT bm1(5, 5, tri1);
    SPARSE_MATRIX_INT bm2(3, 3, tri2);
    cout << "Matrix 1:" << endl;
    bm1.show();
    cout << "Matrix 2:" << endl;
    bm2.show();
    SPARSE_MATRIX_INT* mat_ans = bm1.convolve(bm2, 1, 0);
    cout << " Answer of matrix convolution(stride=1, padding=0)" << endl;
    mat_ans->show();
    
    cout << "Test end." << endl;
}