#ifndef MATRIX_H
#define MATRIX_H

#include <ccomplex>
#include <cstring>
#include <opencv2/opencv.hpp>

#include "basic-matrix.h"
#include "sparse-matrix.h"

template<class T>
class Matrix {
private:
    int size;
    int row;
    int col;
public:

    Matrix(int row, int col) {
        this->row = row;
        this->col = col;
        this->size = this->row * this->col;
    }

    Matrix(const cv::Mat &mat) {
        this->row = mat.rows;
        this->col = mat.cols;
        this->size = this->row * this->col;
    }

    void setSize(const int size) {
        this->size = size;
    }

    void setRow(const int row) {
        this->row = row;
    }

    void setCol(const int col) {
        this->col = col;
    }

    int getSize() {
        return this->size;
    }

    int getRow() {
        return this->row;
    }

    int getCol() {
        return this->col;
    }

    virtual Matrix<T> &add(const BasicMatrix<T> &) = 0;

    virtual Matrix<T> &add(const SparseMatrix<T> &) = 0;

    virtual Matrix<T> &subtract(const BasicMatrix<T> &) = 0;

    virtual Matrix<T> &subtract(const SparseMatrix<T> &) = 0;

    virtual Matrix<T> &scalarMultiply(short) = 0;

    virtual Matrix<T> &scalarMultiply(int) = 0;

    virtual Matrix<T> &scalarMultiply(long) = 0;

    virtual Matrix<T> &scalarMultiply(long long) = 0;

    virtual Matrix<T> &scalarMultiply(float) = 0;

    virtual Matrix<T> &scalarMultiply(double) = 0;

    virtual Matrix<T> &scalarMultiply(long double) = 0;

    virtual Matrix<T> &scalarMultiply(std::complex<short>) = 0;

    virtual Matrix<T> &scalarMultiply(std::complex<int>) = 0;

    virtual Matrix<T> &scalarMultiply(std::complex<long>) = 0;

    virtual Matrix<T> &scalarMultiply(std::complex<long long>) = 0;

    virtual Matrix<T> &scalarMultiply(std::complex<float>) = 0;

    virtual Matrix<T> &scalarMultiply(std::complex<double>) = 0;

    virtual Matrix<T> &scalarMultiply(std::complex<long double>) = 0;

    virtual Matrix<T> &scalarDivide(short) = 0;

    virtual Matrix<T> &scalarDivide(int) = 0;

    virtual Matrix<T> &scalarDivide(long) = 0;

    virtual Matrix<T> &scalarDivide(long long) = 0;

    virtual Matrix<T> &scalarDivide(float) = 0;

    virtual Matrix<T> &scalarDivide(double) = 0;

    virtual Matrix<T> &scalarDivide(long double) = 0;

    virtual Matrix<T> &scalarDivide(std::complex<short>) = 0;

    virtual Matrix<T> &scalarDivide(std::complex<int>) = 0;

    virtual Matrix<T> &scalarDivide(std::complex<long>) = 0;

    virtual Matrix<T> &scalarDivide(std::complex<long long>) = 0;

    virtual Matrix<T> &scalarDivide(std::complex<float>) = 0;

    virtual Matrix<T> &scalarDivide(std::complex<double>) = 0;

    virtual Matrix<T> &scalarDivide(std::complex<long double>) = 0;

    // virtual Matrix &dotProduct(const BasicMatrix &) = 0;

    // virtual Matrix &dotProduct(const SparseMatrix &) = 0;

    // virtual Matrix &crossProduct(const BasicMatrix &) = 0;

    // virtual Matrix &crossProduct(const SparseMatrix &) = 0;

//
//    virtual Matrix &transpose() = 0;
//
//    virtual Matrix &inverse() = 0;
//
//    virtual Matrix &conjugate() = 0;
//
//    virtual T getMax() = 0;
//
//    virtual T getMin() = 0;
//
//    virtual T getSum() = 0;
//
//    virtual T getAvg() = 0;
//
//    virtual T getMaxRow(int) = 0;
//
//    virtual T getMaxCol(int) = 0;
//
//    virtual T getMinRow(int) = 0;
//
//    virtual T getMinCol(int) = 0;
//
//    virtual T getSumRow(int) = 0;
//
//    virtual T getSumCol(int) = 0;
//
//    virtual T getAvgRow(int) = 0;
//
//    virtual T getAvgCol(int) = 0;
//
//    virtual T getEigenvalue() = 0;
//
//    virtual Matrix& getEigenvector() = 0;
//
//    virtual T getTrace() = 0;
//
//    virtual T getDeterminant() = 0;
//
//    virtual void reshape(int row, int col) = 0;
//
//    virtual void sliceRow(int row1, int row2) = 0;
//
//    virtual void sliceCol(int col1, int col2) = 0;
//
//    virtual void slice(int row1, int row2, int col1, int col2) = 0;
//
//    virtual Matrix& convolve(Matrix&) = 0;

};


#endif