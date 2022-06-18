#ifndef MATRIX_H
#define MATRIX_H

#include <ccomplex>
#include <cstring>
// #include <opencv2/opencv.hpp>

#include "basic-matrix.h"
#include "sparse-matrix.h"

namespace mat {

    #define MATRIX_TYPE bool
    #define BASIC_MATRIX true
    #define SPARSE_MATRIX false

 namespace ex {
        class MatrixException;
        class MismatchedSizeException;
        class DuplicatedTripleException;
        class NotSquareException;
        class NoInverseException;
        class InvalidSizeException;
    }

    template<class T>
    class Matrix {
    private:
        long size;
    protected:
        int row;
        int col;
    public:

        Matrix(int row, int col) {
            if (row <= 0 || col <= 0) {
                throw ex::InvalidSizeException("creating Matrix", 1, row, col);
            }
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

        long getSize() const{
            return this->size;
        }

        int getRow() const{
            return this->row;
        }

        int getCol() const{
            return this->col;
        }

        virtual T getByIndex(int _row, int _col) const = 0;

        virtual void setByIndex(int _row, int _col, T val) = 0;

        virtual void add(const BasicMatrix<T> &) = 0;

        virtual void add(const SparseMatrix<T> &) = 0;

        virtual void subtract(const BasicMatrix<T> &) = 0;

        virtual void subtract(const SparseMatrix<T> &) = 0;

        virtual void scalarMultiply(T) = 0;

        virtual void scalarDivide(T) = 0;

        virtual void dotProduct(const BasicMatrix<T> &) = 0;

        virtual void dotProduct(const SparseMatrix<T> &) = 0;

        virtual void crossProduct(const BasicMatrix<T> &) = 0;

        virtual void crossProduct(const SparseMatrix<T> &) = 0;

        virtual void transpose() = 0;

        virtual void inverse() = 0;

        virtual void reverse() = 0;

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

        virtual T getTrace() = 0;

        virtual T getDeterminant() = 0;

        virtual void reshape(int row, int col) = 0;

        virtual void sliceRow(int row1, int row2) = 0;

        virtual void sliceCol(int col1, int col2) = 0;

        virtual void slice(int row1, int row2, int col1, int col2) = 0;

        virtual Matrix<T> &convolve(BasicMatrix<T> &, int stride = 1, int padding = 0) = 0;

        virtual Matrix<T> &convolve(SparseMatrix<T> &, int stride = 1, int padding = 0) = 0;

        virtual void exponent(int) = 0;

        virtual void show() {
            std::cout << "Base class Matrix" << std::endl;
        }

        static Matrix<T> * eye(int row, int col, MATRIX_TYPE type){
            if (row <= 0 || col <= 0) {
                throw ex::InvalidSizeException("creating Matrix", 1, row, col);
            }
            if (type){
                BasicMatrix<T> mat(row, col);
                for (size_t i = 0; i < row; i++)
                    mat.setByIndex(i, i, 1);
                return &mat;
            }
        }
    };

}

#endif