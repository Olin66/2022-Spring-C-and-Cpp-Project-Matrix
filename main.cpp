#include <iostream>
#include <vector>
#include "matrix.h"
#include "basic-matrix.h"
#include <opencv2/opencv.hpp>
using namespace std;
using namespace mat;
using namespace cv;
extern void bmat_constructor_test();
extern void bmat_arithmetic_test();
extern void bmat_reshape_slice_test();
extern void smat_constructor_test();
extern void smat_bigdata_test();
extern void mat_opencv_test();
extern void bmat_multi_test1();
extern void bmat_multi_test2();
extern void bmat_multi_test3();
extern void bmat_multi_test4();
extern void smat_multi_test1();
extern void smat_multi_test2();
extern void smat_multi_test3();
extern void bmat_conv_test1();
extern void smat_conv_test1();
extern void bmat_max_test();
extern void bmat_min_test();
extern void bmat_average_test();
extern void bmat_eigenvalue_vector_test();
extern void bmat_transpose_text();
extern void bmat_Jacobi_eigen_test();
extern void bmat_trace_test();
extern void bmat_det_test();
extern void bmat_inverse_test();
extern void smat_max_test();
extern void smat_min_test();
extern void smat_average_test();
extern void smat_eigenvalue_vector_test();
extern void smat_transposition_text();
extern void smat_Jacobi_eigen_test();
extern void smat_trace_test();
extern void smat_det_test();
extern void smat_inverse_test();

int main() {
    smat_bigdata_test();
}
