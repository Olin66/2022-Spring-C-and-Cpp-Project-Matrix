#include "matrix-ex.h"
#include "matrix.h"
#include "sparse-matrix.h"
#include "basic-matrix.h"
#include <opencv2/opencv.hpp>
#include <cstring>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <complex>

using namespace std;
using namespace mat;
using namespace cv;

#define BASIC_MATRIX_INT BasicMatrix<int>
#define BASIC_MATRIX_DOUBLE BasicMatrix<double>
#define BASIC_MATRIX_COMPLEX_INT BasicMatrix<complex<int>>
#define SPARSE_MATRIX_INT SparseMatrix<int>
#define SPARSE_MATRIX_DOUBLE SparseMatrix<double>
#define SPARSE_MATRIX_COMPLEX_INT SparseMatrix<complex<int>>


void bmat_constructor_test(){
    cout<<"Basic matrix constructor test begin: "<<endl;
    
    cout<<"(1) ";
    BASIC_MATRIX_INT bm1(3, 3);
    bm1.show();
    cout<<endl;

    cout<<"(2) ";
    BASIC_MATRIX_DOUBLE bm2(2,2,1.3);
    bm2.show();
    cout<<endl;

    cout<<"(3) ";
    complex<int>* arr1 = new complex<int>[15];
    for (int i = 0; i < 15; i++)
        arr1[i] = complex<int>(1, i);
    BASIC_MATRIX_COMPLEX_INT bm3(5, 3, arr1);
    bm3.show();
    cout<<endl;

    cout<<"(4) ";
    BASIC_MATRIX_DOUBLE bm4(bm2);
    bm4.show();
    cout<<endl;

    cout<<"(5) ";
    vector<vector<double>> v1(2, vector<double>(4));
    v1[0][1] = 5.7;
    v1[0][3] = 6.2;
    v1[1][1] = 0.1;
    v1[0][0] = 3.2;
    v1[1][2] = 10.55;
    BASIC_MATRIX_DOUBLE bm5(v1);
    bm5.show();
    cout<<endl;

    cout<<"(6) ";
    bm4 = bm5;
    bm4.show();
    cout<<endl;

    cout<<"test end"<<endl<<endl;
}

void bmat_arithmetic_test(){
    cout<<"Basic matrix arithmetic test begin: "<<endl;
    cout<<"(ori 1) ";
    BASIC_MATRIX_DOUBLE bm1(2, 4);
    bm1.show();
    cout<<"(new 1) ";
    vector<vector<double>> v1(2, vector<double>(4));
    v1[0][1] = 5.7;
    v1[0][3] = 6.2;
    v1[1][1] = 0.1;
    v1[0][0] = 3.2;
    v1[1][2] = 10.55;
    BASIC_MATRIX_DOUBLE bm2(v1);
    cout<<"add ";
    bm2.show();
    bm1.add(bm2);
    bm1.show();
    cout<<endl;

    cout<<"(ori 2) ";
    BASIC_MATRIX_COMPLEX_INT bm3(2, 2, complex<int>(1, 2));
    bm3.show();
    cout<<"(new 2) ";
    vector<vector<complex<int>>> v2(2, vector<complex<int>>(2));
    v2[0][1] = complex<int>(2, 3);
    v2[1][1] = complex<int>(7, 7);
    BASIC_MATRIX_COMPLEX_INT bm4(v2);
    cout<<"add ";
    bm4.show();
    bm3.add(bm4);
    bm3.show();
    cout<<endl;

    cout<<"(ori 3) ";
    bm2.show();
    vector<vector<double>> v3(2, vector<double>(4));
    v3[0][1] = 5.7;
    v3[0][3] = -3.2;
    v3[1][1] = 0.08;
    v3[0][0] = 3.22;
    v3[1][2] = -1.55;
    v3[1][3] = 3.4;
    BASIC_MATRIX_DOUBLE bm5(v3);
    cout<<"subtract ";
    bm5.show();
    bm2.subtract(bm5);
    cout<<"(new 3) ";
    bm2.show();
    cout<<endl;

    cout<<"(ori 4) "; 
    bm4.show();
    bm4.scalarMultiply(3);
    cout<<"multiply 3"<<endl;
    bm4.show();
    cout<<endl;

    cout<<"(ori 5) "<<endl;
    bm5.show();
    bm5.scalarDivide(2.1);
    cout<<"divide 2.1"<<endl;
    bm5.show();
    cout<<endl;

    cout<<"(ori 6-1) ";
    bm1.show();
    cout<<"(ori 6-2) ";
    bm2.show();
    cout<<"(new 6) ";
    BASIC_MATRIX_DOUBLE bm6 = bm1 + bm2;
    bm6.show();
    cout<<endl;

    cout<<"test end"<<endl<<endl;
}

void bmat_reshape_slice_test(){
    cout<<"Basic matrix reshape and slice test begin: "<<endl;
    cout<<"(ori 1) ";
    vector<vector<double>> v1(2, vector<double>(4));
    v1[0][1] = 5.7;
    v1[0][3] = 6.2;
    v1[1][1] = 0.1;
    v1[0][0] = 3.2;
    v1[1][2] = 10.55;
    BASIC_MATRIX_DOUBLE bm1(v1);
    bm1.show();
    cout<<"reshape to row=4 col=2"<<endl;
    bm1.reshape(4, 2);
    cout<<"(new 1) ";
    bm1.show();
    cout<<endl;

    cout<<"(ori 2) ";
    bm1.show();
    cout<<"slice row to (1,2)"<<endl;
    bm1.sliceRow(1, 2);
    cout<<"(new 2) ";
    bm1.show();
    cout<<endl;

    cout<<"(ori 3) ";
    bm1.show();
    cout<<"slice col to (1,1)"<<endl;
    bm1.sliceCol(1, 1);
    cout<<"(new 3) ";
    bm1.show();
    cout<<endl;

    cout<<"(ori 4) ";
    vector<vector<complex<int>>> v2(4, vector<complex<int>>(3));
    v2[0][1] = complex<int>(2, 3);
    v2[1][1] = complex<int>(7, 7);
    v2[2][2] = complex<int>(1, 2);
    BASIC_MATRIX_COMPLEX_INT bm2(v2);
    bm2.show();
    cout<<"slice row to (1, 3) col to (2, 2)"<<endl;
    bm2.slice(1,3,2,2);
    cout<<"(new 4) ";
    bm2.show();
    cout<<endl;

    cout<<"test end"<<endl<<endl;
}

void mat_opencv_test(){
    cout<<"Basic matrix opencv test begin: "<<endl;
    Mat img = imread("./img/btn1.png");
    Mat grey;
    cvtColor(img, grey, CV_BGR2GRAY);
    BASIC_MATRIX_INT bm1(grey);
    cout<<"(1) ";
    bm1.show();
    BASIC_MATRIX_INT bm2(50, 50, 240);
    Mat* mat1 = bm2.getCvMat();
    imwrite("./img/test1.png", *mat1);
    Matrix<int>* bm3 = Matrix<int>::eye(30, 30, true);
    for (size_t i = 0; i < bm3->getRow(); i++)
        bm3->setByIndex(4, i, 200);
    Mat* mat2 = bm3->getCvMat();
    imwrite("./img/test2.png", *mat2);
    imwrite("./img/test3.png", *bm1.getCvMat());
    SparseMatrix<int> bm4(100, 100);
    Mat* mat3 = bm4.getCvMat();
    imwrite("./img/test4.png", *mat3);

    cout<<"test end"<<endl<<endl;
}

void smat_constructor_test(){
    cout<<"Sparse matrix constructor test begin: "<<endl;
    SPARSE_MATRIX_INT sm1(2, 2);
    cout<<"(1) ";
    sm1.show();
    cout<<endl;

    cout<<"(2) ";
    complex<int>* arr1 = new complex<int>[15];
    for (int i = 0; i < 15; i++)
        arr1[i] = complex<int>(1, i);
    SPARSE_MATRIX_COMPLEX_INT sm2(5, 3, arr1);
    sm2.show();
    cout<<endl;

    vector<Triple<int>> tri1;
    tri1.push_back(Triple<int>(1, 2, 1));
    tri1.push_back(Triple<int>(4, 4, 1));
    tri1.push_back(Triple<int>(1, 1, 1));
    tri1.push_back(Triple<int>(2, 1, 1));
    tri1.push_back(Triple<int>(3, 1, 1));
    SPARSE_MATRIX_INT sm3(5, 5, tri1);
    cout<<"(3) ";
    sm3.show();
    cout<<endl;
}

void smat_bigdata_test(){
    cout << "Big sparse matrix test begin:" << endl;
    vector<Triple<int>> tri1;
    tri1.push_back(Triple<int>(1, 2, 1));
    tri1.push_back(Triple<int>(4, 4, 4));
    vector<Triple<int>> tri2;
    tri2.push_back(Triple<int>(1, 1, 1));
    tri2.push_back(Triple<int>(2, 1, 1));
    tri2.push_back(Triple<int>(3, 1, 1));
    tri2.push_back(Triple<int>(4, 6, 9));
    SparseMatrix<int> *bm1 = new SparseMatrix<int>(50, 5000, tri1);
    SparseMatrix<int> *bm2 = new SparseMatrix<int>(5000, 50, tri2);
    bm2->crossProduct(*bm1);
    bm2->show();
    cout << "Test end." << endl << endl;
}