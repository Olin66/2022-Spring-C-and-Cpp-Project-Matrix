#include <iostream>
#include <vector>
#include "matrix.h"
#include "basic-matrix.h"
#define BASIC_MATRIX_INT BasicMatrix<int>
#define SPARSE_MATRIX_INT SparseMatrix<int>
using namespace std;
using namespace mat;

void matrix_dot_product() {
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
}

void matrix_cross_product() {
    const int TEST_SIZE = 6;
    int *_data = new int[TEST_SIZE]{1,2,3,4,5,6};
    BASIC_MATRIX_INT bm1(3, 2, _data);
    BASIC_MATRIX_INT bm2(2, 3, _data);
    bm1.show();
    bm2.show();
    bm1.crossProduct(bm2);
    bm1.show();
}

int main() {
    // matrix_cross_product();
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
    

    SparseMatrix<int> sm1(v);
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
}