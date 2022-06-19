#include <iostream>
#include <vector>
#include "matrix.h"
#include "basic-matrix.h"
#include <opencv2/opencv.hpp>
#define BASIC_MATRIX_INT SparseMatrix<int>//BasicMatrix<int>
#define SPARSE_MATRIX_INT SparseMatrix<int>
using namespace std;
using namespace mat;
using namespace cv;
void matrix_dot_product();
void matrix_cross_product();
void matrix_convolution();

extern void bmat_constructor_test();
extern void bmat_add_test();
extern void bmat_multi_test1();
extern void bmat_multi_test2();
extern void bmat_multi_test3();
extern void smat_multi_test1();
extern void smat_multi_test2();
extern void smat_multi_test3();
extern void bmat_conv_test1();

int main() {
    // bmat_constructor_test();

    return 0;
    // matrix_dot_product();
    // matrix_cross_product();
    // //matrix_convolution();
    // return 0;

    Matrix<int> *m;
    BasicMatrix<int> bm1(3, 2);
    BasicMatrix<int> bm2(3, 2);
    BasicMatrix<int> rm(1, 2);

    rm.show();

    m = &bm1;
    cout<<bm1.getCol()<<endl;
    cout<<m->getCol()<<endl;

//    cout << bm.getData()[1] << endl;
//    cout << bm.getData()[2] << endl;
//    cout << bm.getData()[3] << endl;
//    cout << bm.getData()[4] << endl;
//    cout << bm.getData()[5] << endl;

    vector<vector<int>> v(3, vector<int>(2));
    v[0][1] = 5;
    v[2][1] = 6;
    v[1][1] = 1;
    v[0][0] = 3;
    BasicMatrix<int> qm(v);
    cout<<qm.getSize()<<endl;

    try {
        bm1.add(rm);
    }catch (ex::MismatchedSizeException& e){
        cout<<e.what()<<endl;
    }

    bm1.add(qm);

    BasicMatrix<int> bm3 = bm1 + qm;

    bm3.show();

    BasicMatrix<int> mm(bm1);

    // BasicMatrix<int> a = bm * qm;

    // for (int i = 0; i < bm1.getSize(); i++)
    // {
    //     cout<<mm.getData()[i]<<endl;;
    // }
    
    // sm1.show();

    cout<<endl;
    cout<<endl;
    cout<<endl;
    cout<<endl;
    BasicMatrix<int> m1(v);
    m1.show();
    BasicMatrix<int> m2(v);
    BasicMatrix<int> m3 = m1 + m2;
    m3.show();

    m3 = m3 * 2;

    m3.show();

    m3 = 2 * m3;
    m3.show();

    BasicMatrix<int> m4(2, 2, 100);
    m4.show();

    vector<Triple<int>> v3;
    v3.push_back(Triple<int>(1, 1, 4));
    v3.push_back(Triple<int>(1, 0, 1));
    // v3.push_back(Triple<int>(1, 1, 2));
    // v3.push_back(Triple<int>(2, 2, 4));
    SparseMatrix<int> sm2(2, 2, v3);
    sm2.show();
    SparseMatrix<int> sm3(sm2);
    // cout<<sm3.getTriples()[3]->_row<<endl;
    // cout<<sm3.getTriples()[3]->_col<<endl;
    // cout<<sm3.getTriples()[3]->val<<endl;
    sm3.show();

    SparseMatrix<int> sm4(2, 2);
    sm4 = sm3;
    sm4.show();

    Mat mat1(5, 5, CV_8UC1);
    BasicMatrix<int> mmm1(mat1);
    Mat mat2 = imread("./img/conan1.png");
    Mat grey;
    cvtColor(mat2,grey,CV_BGR2GRAY);
    // imshow("grey",grey);
	BasicMatrix<int> mm2(grey);
    // mm2.show();
    SparseMatrix<int> mm3(grey);
    mm3.show();
}

void matrix_dot_product() {
    cout << "matrix dot product test begins:" << endl;
    const int TEST_SIZE = 6;
    int *_data = new int[TEST_SIZE]{1,2,3,4,5,6};
    BASIC_MATRIX_INT bm1(3, 2, _data);
    BASIC_MATRIX_INT bm2(3, 2, _data);
    bm1.show();
    bm1.add(bm2);
    bm1.show();
    bm1.dotProduct(bm2);
    bm1.show();
    int *_data1 = new int[TEST_SIZE/2]{1,2,3};
    BASIC_MATRIX_INT bm3(3, 1, _data1);
    bm3.show();
    bm2.show();
    bm3.dotProduct(bm2);
    bm3.show();
    cout << endl;
}

void matrix_cross_product() {
    cout << "matrix cross product test begins:" << endl;
    const int TEST_SIZE = 6;
    int *_data = new int[TEST_SIZE]{1,2,3,4,5,6};
    BASIC_MATRIX_INT bm1(3, 2, _data);
    BASIC_MATRIX_INT bm2(2, 3, _data);
    bm1.show();
    bm2.show();
    bm1.crossProduct(bm2); 
    bm1.show();
    cout << endl;
}

void matrix_convolution() {
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
}