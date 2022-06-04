#include <iostream>
#include <vector>
#include "matrix.h"
#include "basic-matrix.h"

using namespace std;

int main(){
    Matrix<int>* m;
    BasicMatrix<int> bm(1, 2);
    BasicMatrix<int> rm(1,2);
    m = &bm;
    cout<<bm.getCol()<<endl;
    cout<<m->getCol()<<endl;

    vector<vector<int>> v(3, vector<int>(2));
    v[0][1] = 5;
    v[2][1] = 6;
    v[1][1] = 1;
    v[0][0] = 3;
    BasicMatrix<int> qm(v);
}