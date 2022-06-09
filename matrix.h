#ifndef MATRIX_H
#define MATRIX_H

#include <ccomplex>
#include <cstring>
// #include <opencv2/opencv.hpp>

#include "basic-matrix.h"
#include "sparse-matrix.h"

namespace mat {

    template<class T>
    class Matrix {
    private:
        long size;
        int row;
        int col;
    public:

        Matrix(int row, int col) {
            this->row = row;
            this->col = col;
            this->size = this->row * this->col;
        }

        // Matrix(const cv::Mat &mat) {
        //     this->row = mat.rows;
        //     this->col = mat.cols;
        //     this->size = this->row * this->col;
        // }

        void setSize(const long size) {
            this->size = size;
        }

        void setRow(const int row) {
            this->row = row;
        }

        void setCol(const int col) {
            this->col = col;
        }

        int getSize() const{
            return this->size;
        }

        int getRow() const{
            return this->row;
        }

        int getCol() const{
            return this->col;
        }

        virtual void add(const BasicMatrix<T> &) = 0;

        virtual void add(const SparseMatrix<T> &) = 0;

        virtual void subtract(const BasicMatrix<T> &) = 0;

        virtual void subtract(const SparseMatrix<T> &) = 0;

        virtual void scalarMultiply(short) = 0;

        virtual void scalarMultiply(int) = 0;

        virtual void scalarMultiply(long) = 0;

        virtual void scalarMultiply(long long) = 0;

        virtual void scalarMultiply(float) = 0;

        virtual void scalarMultiply(double) = 0;

        virtual void scalarMultiply(long double) = 0;

        virtual void scalarMultiply(std::complex<short>) = 0;

        virtual void scalarMultiply(std::complex<int>) = 0;

        virtual void scalarMultiply(std::complex<long>) = 0;

        virtual void scalarMultiply(std::complex<long long>) = 0;

        virtual void scalarMultiply(std::complex<float>) = 0;

        virtual void scalarMultiply(std::complex<double>) = 0;

        virtual void scalarMultiply(std::complex<long double>) = 0;

        virtual void scalarDivide(short) = 0;

        virtual void scalarDivide(int) = 0;

        virtual void scalarDivide(long) = 0;

        virtual void scalarDivide(long long) = 0;

        virtual void scalarDivide(float) = 0;

        virtual void scalarDivide(double) = 0;

        virtual void scalarDivide(long double) = 0;

        virtual void scalarDivide(std::complex<short>) = 0;

        virtual void scalarDivide(std::complex<int>) = 0;

        virtual void scalarDivide(std::complex<long>) = 0;

        virtual void scalarDivide(std::complex<long long>) = 0;

        virtual void scalarDivide(std::complex<float>) = 0;

        virtual void scalarDivide(std::complex<double>) = 0;

        virtual void scalarDivide(std::complex<long double>) = 0;

        virtual void dotProduct(const BasicMatrix<T> &) = 0;

        virtual void dotProduct(const SparseMatrix<T> &) = 0;

        virtual void crossProduct(const BasicMatrix<T> &) = 0;

        virtual void crossProduct(const SparseMatrix<T> &) = 0;

        virtual void transpose() = 0;

        virtual void inverse() = 0;

        virtual void conjugate() = 0;

        virtual T getMax() = 0;

        virtual T getMin() = 0;

        virtual T getSum() = 0;

        virtual T getAvg() = 0;

        virtual T getRowMax(int) = 0;

        virtual T getColMax(int) = 0;

        virtual T getRowMin(int) = 0;

        virtual T getColMin(int) = 0;

        virtual T getRowSum(int) = 0;

        virtual T getColSum(int) = 0;

        virtual T getRowAvg(int) = 0;

        virtual T getColAvg(int) = 0;

        virtual T getEigenvalue() = 0;

        virtual Matrix<T> &getEigenvector() = 0;

        virtual T getTrace() = 0;

        virtual T getDeterminant() = 0;

        virtual void reshape(int row, int col) = 0;

        virtual void sliceRow(int row1, int row2) = 0;

        virtual void sliceCol(int col1, int col2) = 0;

        virtual void slice(int row1, int row2, int col1, int col2) = 0;

        virtual Matrix<T> &convolve(BasicMatrix<T> &) = 0;

        virtual Matrix<T> &convolve(SparseMatrix<T> &) = 0;

    };

}

#endif