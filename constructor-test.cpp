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

void bmat_constructor_test(){
    cout<<"Basic matrix constructor test begin: "<<endl;
    
    cout<<"(1) ";
    BASIC_MATRIX_INT bm1(3, 3);
    bm1.show();

    cout<<"(2) ";
    BASIC_MATRIX_DOUBLE bm2(2,2,1.3);
    bm2.show();

    cout<<"(3) ";
    int* arr1 = new int[15];
    for (int i = 0; i < 15; i++)
        arr1[i] = i * 2;
    BASIC_MATRIX_INT bm3(5, 3, arr1);
    bm3.show();

    cout<<"(4) ";
    BASIC_MATRIX_DOUBLE bm4(bm2);
    bm4.show();

    cout<<"(5) ";
    vector<vector<double>> v1(2, vector<double>(4));
    v1[0][1] = 5.7;
    v1[0][3] = 6.2;
    v1[1][1] = 0.1;
    v1[0][0] = 3.2;
    v1[1][2] = 10.55;
    BASIC_MATRIX_DOUBLE bm5(v1);
    bm5.show();

    cout<<"(6) ";
    bm4 = bm5;
    bm4.show();

    cout<<"test end"<<endl<<endl;
}

void bmat_add_test(){
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
    bm1.add(bm2);
    bm1.show();

    
}

int main(){
    bmat_constructor_test();
    bmat_add_test();
}