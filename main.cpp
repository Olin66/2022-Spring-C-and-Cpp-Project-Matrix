#include <iostream>
#include <vector>
#include "matrix.h"
#include "basic-matrix.h"

using namespace std;
using namespace mat;

int main() {
    Matrix<int> *m;
    BasicMatrix<int> bm(3, 2);
    BasicMatrix<int> rm(1, 2);
    m = &bm;
    cout<<bm.getCol()<<endl;
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
        bm.add(rm);
    }catch (ex::MismatchedSizeException& e){
        cout<<e.what()<<endl;
    }
}