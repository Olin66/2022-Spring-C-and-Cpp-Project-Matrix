#ifndef MATRIX_H
#define MATRIX_H

#include <ccomplex>
#include <cstring>
#include <opencv2/opencv.hpp>

template<class T>
class Matrix {
private:
    int row;
    int col;
public:

    Matrix(int row, int col) {
        this->row = row;
        this->col = col;
    }

    Matrix(cv::Mat mat){
        this->row = mat.rows;
        this->col = mat.cols;
    }

    void setRow(int row) {
        this->row = row;
    }

    void setCol(int col) {
        this->col = col;
    }

    int getRow() {
        this->row;
    }

    int getCol() {
        this->col;
    }

    virtual Matrix &add(const Matrix &) = 0;

    virtual Matrix &subtract(const Matrix &) = 0;

    virtual Matrix &scalarMultiply(short) = 0;

    virtual Matrix &scalarMultiply(int) = 0;

    virtual Matrix &scalarMultiply(long) = 0;

    virtual Matrix &scalarMultiply(long long) = 0;

    virtual Matrix &scalarMultiply(float) = 0;

    virtual Matrix &scalarMultiply(double) = 0;

    virtual Matrix &scalarMultiply(long double) = 0;

    virtual Matrix &scalarMultiply(std::complex<short>) = 0;

    virtual Matrix &scalarMultiply(std::complex<int>) = 0;

    virtual Matrix &scalarMultiply(std::complex<long>) = 0;

    virtual Matrix &scalarMultiply(std::complex<long long>) = 0;

    virtual Matrix &scalarMultiply(std::complex<float>) = 0;

    virtual Matrix &scalarMultiply(std::complex<double>) = 0;

    virtual Matrix &scalarMultiply(std::complex<long double>) = 0;

    virtual Matrix &scalarDivide(short) = 0;

    virtual Matrix &scalarDivide(int) = 0;

    virtual Matrix &scalarDivide(long) = 0;

    virtual Matrix &scalarDivide(long long) = 0;

    virtual Matrix &scalarDivide(float) = 0;

    virtual Matrix &scalarDivide(double) = 0;

    virtual Matrix &scalarDivide(long double) = 0;

    virtual Matrix &scalarDivide(std::complex<short>) = 0;

    virtual Matrix &scalarDivide(std::complex<int>) = 0;

    virtual Matrix &scalarDivide(std::complex<long>) = 0;

    virtual Matrix &scalarDivide(std::complex<long long>) = 0;

    virtual Matrix &scalarDivide(std::complex<float>) = 0;

    virtual Matrix &scalarDivide(std::complex<double>) = 0;

    virtual Matrix &scalarDivide(std::complex<long double>) = 0;

    virtual Matrix &dotProduct(const Matrix &) = 0;

    virtual Matrix &crossProduct(const Matrix &) = 0;

    virtual Matrix &transpose() = 0;

    virtual Matrix &inverse() = 0;

    virtual Matrix &conjugate() = 0;

    virtual T getMax() = 0;

    virtual T getMin() = 0;

    virtual T getSum() = 0;

    virtual T getAvg() = 0;

    virtual T getMaxRow(int) = 0;

    virtual T getMaxCol(int) = 0;

    virtual T getMinRow(int) = 0;

    virtual T getMinCol(int) = 0;

    virtual T getSumRow(int) = 0;

    virtual T getSumCol(int) = 0;

    virtual T getAvgRow(int) = 0;

    virtual T getAvgCol(int) = 0;

    virtual T getEigenvalue() = 0;

    virtual Matrix& getEigenvector() = 0;

    virtual T getTrace() = 0;

    virtual T getDeterminant() = 0;

    virtual void reshape(int row, int col) = 0;

    virtual void slice(int row1, int row2) = 0;

    virtual void slice(int col1, int col2) = 0;

    virtual void slice(int row1, int row2, int col1, int col2) = 0;

    virtual Matrix& convolve(Matrix&) = 0;

};


#endif