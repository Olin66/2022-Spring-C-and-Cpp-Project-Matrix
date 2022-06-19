#include"matrix.h"
#include"matrix-ex.h"
#include"basic-matrix.h"
#include"sparse-matrix.h"

using namespace mat;
// void bmat_max_test();
// void bmat_min_test();
// void bmat_average_test();
// void bmat_eigenvalue_vector_test();
// void bmat_transposition_text();
// void bmat_Jacobi_eigen_test();//使用雅克比算法求实对称矩阵的特征值和特征向量。
// void bmat_trace_test();
// void bmat_det_test();
// void bmat_inverse_test();
// void smat_max_test();
// void smat_min_test();
// void smat_average_test();
// void smat_eigenvalue_vector_test();
// void smat_transposition_text();
// void smat_Jacobi_eigen_test();//使用雅克比算法求实对称矩阵的特征值和特征向量。
// void smat_trace_test();
// void smat_det_test();
// void smat_inverse_test();
// void exception_test1();
// void exception_test2();


void bmat_transpose_text() {
    cout << "Basic matrix transposition begin" << endl;
    BasicMatrix<int> test(3, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            test.setByIndex(i, j, i * j);
        }

    }
    cout << "Matrix :" << endl;
    test.show();
    test.transpose();
    cout << "Transpose matrix is: " << endl;
    test.show();
    cout << "test end" << endl << endl;
}

void bmat_max_test() {
    cout << "Basic matrix max test begin:" << endl;
    BasicMatrix<int> test(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            test.setByIndex(i, j, i + j);
        }

    }
    cout << "Matrix:" << endl;
    test.show();
    cout << "The max value is " << test.getMax() << endl;
    cout << "The max value of first row is " << test.getRowMax(0) << endl;
    cout << "The max value of first col is " << test.getColMax(0) << endl;
    cout << "test end" << endl << endl;

}

void bmat_min_test() {
    cout << "Basic matrix min test begin:" << endl;
    BasicMatrix<int> test(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            test.setByIndex(i, j, i + j);
        }

    }
    cout << "Matrix:" << endl;
    test.show();
    cout << "The min value is " << test.getMin() << endl;
    cout << "The min value of first row is " << test.getRowMin(0) << endl;
    cout << "The min value of first col is " << test.getColMin(0) << endl;
    cout << "test end" << endl << endl;
}

void bmat_average_test() {
    cout << "Basic matrix average test begin:" << endl;
    BasicMatrix<int> test(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            test.setByIndex(i, j, i + j);
        }

    }
    cout << "Matrix:" << endl;
    test.show();
    cout << "The average is " << test.getAvg() << endl;
    cout << "The average of first row is " << test.getRowAvg(0) << endl;
    cout << "The average of first col is " << test.getColAvg(0) << endl;
    cout << "test end" << endl << endl;
}

void bmat_eigenvalue_vector_test() {
    cout << "Basic matrix eigenvalue test begin:" << endl;
    BasicMatrix<double> test(3, 3);
    test.setByIndex(0, 0, 5.0);
    test.setByIndex(0, 1, -3.0);
    test.setByIndex(0, 2, 2.0);
    test.setByIndex(1, 0, 6.0);
    test.setByIndex(1, 1, -4.0);
    test.setByIndex(1, 2, 4.0);
    test.setByIndex(2, 0, 4.0);
    test.setByIndex(2, 1, -4.0);
    test.setByIndex(2, 2, 5.0);
    cout << "Matrix:" << endl;
    test.show();
    BasicMatrix<double> eigenvalue(3, 3);
    test.getEigenvalue(10, eigenvalue);
    cout << "The eigenvalue is in the diagonal" << endl;
    eigenvalue.show();
    cout << "Basic matirx eigenvecor test begin:" << endl;
    BasicMatrix<double> eigenvector(3, 1);
    for (int i = 0; i < test.getCol(); i++) {
        test.getEigenvector(eigenvector, eigenvalue.getByIndex(i, i));//将所求的特征值依次带入求得特征向量。
        cout << "The No." << i << " eigenvector is:" << endl;
        eigenvector.show();
    }
    cout << "test end" << endl << endl;

}

void bmat_Jacobi_eigen_test() {
    cout << "Basic matrix Jacobi test begin:" << endl;
    BasicMatrix<double> test(3, 3);
    test.setByIndex(0, 0, 2.0);
    test.setByIndex(0, 1, -1.0);
    test.setByIndex(0, 2, 0.0);
    test.setByIndex(1, 0, -1.0);
    test.setByIndex(1, 1, 2.0);
    test.setByIndex(1, 2, -1.0);
    test.setByIndex(2, 0, 0.0);
    test.setByIndex(2, 1, -1.0);
    test.setByIndex(2, 2, 2.0);
    cout << "Matrix:" << endl;
    test.show();
    BasicMatrix<double> eigenvalue(3, 3);
    BasicMatrix<double> eigenvector(3, 3);
    test.getEigen(eigenvector, eigenvalue, 0.00001, 20);
    cout << "The eigenvalues are on diagonal:" << endl;
    eigenvalue.show();
    cout << "Basic matirx eigenvecor test begin:" << endl;
    eigenvector.show();
    cout << "test end" << endl << endl;
}

void bmat_trace_test() {
    cout << "Basic matrix trace test begin:" << endl;
    BasicMatrix<int> test(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            test.setByIndex(i, j, i + j);
        }

    }
    cout << "Matrix:" << endl;
    test.show();
    cout << "The trace is " << test.getTrace() << endl;
    cout << "test end" << endl << endl;
}

void bmat_det_test() {

    cout << "Basic matrix determinant test begin:" << endl;
    BasicMatrix<double> test(3, 3);
    test.setByIndex(0, 0, 5.0);
    test.setByIndex(0, 1, -3.0);
    test.setByIndex(0, 2, 2.0);
    test.setByIndex(1, 0, 6.0);
    test.setByIndex(1, 1, -4.0);
    test.setByIndex(1, 2, 4.0);
    test.setByIndex(2, 0, 4.0);
    test.setByIndex(2, 1, -4.0);
    test.setByIndex(2, 2, 5.0);
    cout << "Matrix:" << endl;
    test.show();
    cout << "The determinant of Matrix is: " << test.getDeterminant() << endl;
    cout << "test end" << endl << endl;
}

void bmat_inverse_test() {
    cout << "Basic matrix inverse test begin:" << endl;
    BasicMatrix<double> test(3, 3);
    test.setByIndex(0, 0, 5.0);
    test.setByIndex(0, 1, -3.0);
    test.setByIndex(0, 2, 2.0);
    test.setByIndex(1, 0, 6.0);
    test.setByIndex(1, 1, -4.0);
    test.setByIndex(1, 2, 4.0);
    test.setByIndex(2, 0, 4.0);
    test.setByIndex(2, 1, -4.0);
    test.setByIndex(2, 2, 5.0);
    cout << "Matrix:" << endl;
    test.show();
    BasicMatrix<double> inver(test);
    inver.inverse();
    cout << "The inverse of Matrix is: " << endl;
    inver.show();
    cout << "The result of A*A^-1 should be I" << endl;
    test.crossProduct(inver);
    test.show();
    cout << "test end" << endl << endl;

}

void exception_test1() {
    cout << " matrix exception test1 begin:" << endl;
    BasicMatrix<int> test(3, 2);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            test.setByIndex(i, j, i + j);
        }

    }
    cout << "Matrix: " << endl;
    test.getDeterminant();
}

void exception_test2() {
    cout << " matrix exception test2 begin:" << endl;
    BasicMatrix<int> test(2, 2);
    test.setByIndex(0, 0, 1);
    test.setByIndex(0, 1, 0);
    test.setByIndex(1, 0, 1);
    test.setByIndex(1, 1, 0);
    cout << "Matrix: " << endl;
    test.show();
    test.inverse();

}

void smat_max_test() {
    cout << "Sparse matrix max test begin:" << endl;
    SparseMatrix<int> test(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            test.setByIndex(i, j, i + j);
        }

    }
    cout << "Matrix:" << endl;
    test.show();
    cout << "The max value is " << test.getMax() << endl;
    cout << "The max value of first row is " << test.getRowMax(0) << endl;
    cout << "The max value of first col is " << test.getColMax(0) << endl;
    cout << "test end" << endl << endl;
}

void smat_min_test() {
    cout << "Sparse matrix min test begin:" << endl;
    SparseMatrix<int> test(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            test.setByIndex(i, j, i + j);
        }

    }
    cout << "Matrix:" << endl;
    test.show();
    cout << "The min value is " << test.getMin() << endl;
    cout << "The min value of first row is " << test.getRowMin(0) << endl;
    cout << "The min value of first col is " << test.getColMin(0) << endl;
    cout << "test end" << endl << endl;

}

void smat_average_test() {
    cout << "Sparse matrix average test begin:" << endl;
    SparseMatrix<int> test(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            test.setByIndex(i, j, i + j);
        }

    }
    cout << "Matrix:" << endl;
    test.show();
    cout << "The average is " << test.getAvg() << endl;
    cout << "The average of first row is " << test.getRowAvg(0) << endl;
    cout << "The average of first col is " << test.getColAvg(0) << endl;
    cout << "test end" << endl << endl;
}

void smat_eigenvalue_vector_test() {
    cout << "Sparse matrix eigenvalue test begin:" << endl;
    SparseMatrix<double> test(3, 3);
    test.setByIndex(0, 0, 5.0);
    test.setByIndex(0, 1, -3.0);
    test.setByIndex(0, 2, 2.0);
    test.setByIndex(1, 0, 6.0);
    test.setByIndex(1, 1, -4.0);
    test.setByIndex(1, 2, 4.0);
    test.setByIndex(2, 0, 4.0);
    test.setByIndex(2, 1, -4.0);
    test.setByIndex(2, 2, 5.0);
    cout << "Matrix:" << endl;
    test.show();
    SparseMatrix<double> eigenvalue(3, 3);
    test.getEigenvalue(10, eigenvalue);
    cout << "The eigenvalue is in the diagonal" << endl;
    eigenvalue.show();
    cout << "Basic matirx eigenvecor test begin:" << endl;
    SparseMatrix<double> eigenvector(3, 1);
    for (int i = 0; i < test.getCol(); i++) {
        test.getEigenvector(eigenvector, eigenvalue.getByIndex(i, i));//将所求的特征值依次带入求得特征向量。
        cout << "The No." << i << " eigenvector is:" << endl;
        eigenvector.show();
    }
    cout << "test end" << endl << endl;


}

void smat_transposition_text() {
    cout << "Sparse matrix transposition begin" << endl;
    SparseMatrix<int> test(3, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            test.setByIndex(i, j, i * j);
        }

    }
    cout << "Matrix :" << endl;
    test.show();
    test.transpose();
    cout << "Transpose matrix is: " << endl;
    test.show();
    cout << "test end" << endl << endl;

}

void smat_Jacobi_eigen_test() {
    cout << "Sparse matrix Jacobi test begin:" << endl;
    SparseMatrix<double> test(3, 3);
    test.setByIndex(0, 0, 2.0);
    test.setByIndex(0, 1, -1.0);
    test.setByIndex(0, 2, 0.0);
    test.setByIndex(1, 0, -1.0);
    test.setByIndex(1, 1, 2.0);
    test.setByIndex(1, 2, -1.0);
    test.setByIndex(2, 0, 0.0);
    test.setByIndex(2, 1, -1.0);
    test.setByIndex(2, 2, 2.0);
    cout << "Matrix:" << endl;
    test.show();
    SparseMatrix<double> eigenvalue(3, 3);
    SparseMatrix<double> eigenvector(3, 3);
    test.getEigen(eigenvector, eigenvalue, 0.00001, 20);
    cout << "The eigenvalues are on diagonal:" << endl;
    eigenvalue.show();
    cout << "Sparse matirx eigenvecor test begin:" << endl;
    eigenvector.show();
    cout << "test end" << endl << endl;
}

void smat_trace_test() {
    cout << "Sparse matrix trace test begin:" << endl;
    SparseMatrix<int> test(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            test.setByIndex(i, j, i + j);
        }

    }
    cout << "Matrix:" << endl;
    test.show();
    cout << "The trace is " << test.getTrace() << endl;
    cout << "test end" << endl << endl;
}

void smat_det_test() {
    cout << "Sparse matrix determinant test begin:" << endl;
    SparseMatrix<double> test(3, 3);
    test.setByIndex(0, 0, 5.0);
    test.setByIndex(0, 1, -3.0);
    test.setByIndex(0, 2, 2.0);
    test.setByIndex(1, 0, 6.0);
    test.setByIndex(1, 1, -4.0);
    test.setByIndex(1, 2, 4.0);
    test.setByIndex(2, 0, 4.0);
    test.setByIndex(2, 1, -4.0);
    test.setByIndex(2, 2, 5.0);
    cout << "Matrix:" << endl;
    test.show();
    cout << "The determinant of Matrix is: " << test.getDeterminant() << endl;
    cout << "test end" << endl << endl;
}

void smat_inverse_test() {
    cout << "Sparse matrix inverse test begin:" << endl;
    SparseMatrix<double> test(3, 3);
    test.setByIndex(0, 0, 5.0);
    test.setByIndex(0, 1, -3.0);
    test.setByIndex(0, 2, 2.0);
    test.setByIndex(1, 0, 6.0);
    test.setByIndex(1, 1, -4.0);
    test.setByIndex(1, 2, 4.0);
    test.setByIndex(2, 0, 4.0);
    test.setByIndex(2, 1, -4.0);
    test.setByIndex(2, 2, 5.0);
    cout << "Matrix:" << endl;
    test.show();
    SparseMatrix<double> inver(test);
    inver.inverse();
    cout << "The inverse of Matrix is: " << endl;
    inver.show();
    cout << "The result of A*A^-1 should be I" << endl;
    test.crossProduct(inver);
    test.show();
    cout << "test end" << endl << endl;
}