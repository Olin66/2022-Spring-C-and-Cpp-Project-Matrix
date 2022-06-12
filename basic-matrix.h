#ifndef BASIC_MATRIX_H
#define BASIC_MATRIX_H

#include <string>
#include <cstring>
#include <vector>
#include "matrix.h"
#include "sparse-matrix.h"
#include "matrix-ex.h"

namespace mat {

    template<typename>
    class Matrix;

    template<typename>
    class SparseMatrix;

    template<class T>
    class BasicMatrix : public Matrix<T> {
    private:
        T *m_data;
    public:
        BasicMatrix(int, int);

        BasicMatrix(int, int, T *);

        // BasicMatrix(const cv::Mat &mat);

        BasicMatrix(std::vector<std::vector<T>>);

        BasicMatrix(const BasicMatrix<T> &);

        BasicMatrix<T>& operator=(const BasicMatrix<T>&);

        void add(const BasicMatrix<T> &);

        void add(const SparseMatrix<T> &);

        void subtract(const BasicMatrix<T> &);

        void subtract(const SparseMatrix<T> &);

        void scalarMultiply(T);

        void scalarDivide(T);

        void dotProduct(const BasicMatrix<T> &);

        void dotProduct(const SparseMatrix<T> &);

        void crossProduct(const BasicMatrix<T> &);

        void crossProduct(const SparseMatrix<T> &);

        void transpose();

        void inverse();

        void conjugate();

        T getByIndex(int, int) const;

        void setByIndex(int, int, T);

        T getMax();

        T getMin();

        T getSum();

        T getAvg();

        T getRowMax(int);

        T getColMax(int);

        T getRowMin(int);

        T getColMin(int);

        T getRowSum(int);

        T getColSum(int);

        T getRowAvg(int);

        T getColAvg(int);

        T getEigenvalue();

        Matrix<T> &getEigenvector();

        T getTrace();

        T getDeterminant();

        void reshape(int row, int col);

        void sliceRow(int row1, int row2);

        void sliceCol(int col1, int col2);

        void slice(int row1, int row2, int col1, int col2);

        Matrix<T> &convolve(BasicMatrix<T> &);

        Matrix<T> &convolve(SparseMatrix<T> &);

        void exponent(int exp);

        void show();

        BasicMatrix<T> operator+(const BasicMatrix<T> &);

        BasicMatrix<T> operator+(const SparseMatrix<T> &);

        BasicMatrix<T> operator-(const BasicMatrix<T> &);

        BasicMatrix<T> operator-(const SparseMatrix<T> &);

        template<typename P>
        BasicMatrix<T> operator*(P);

        T* getData(){
            return this->m_data;
        }

    };

    template<class T>
    BasicMatrix<T>::BasicMatrix(int row, int col):Matrix<T>(row, col) {
        this->m_data = new T[this->getSize()];
        std::memset(m_data, 0, sizeof(T) * this->getSize());
    }

    template<class T>
    BasicMatrix<T>::BasicMatrix(int row, int col, T *_data):Matrix<T>(row, col) {
        this->m_data = new T[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++)
            this->m_data[i] = _data[i];
    }

    // template<class T>
    // BasicMatrix<T>::BasicMatrix(const cv::Mat &mat):Matrix<T>(mat) {

    // }

    template<class T>
    BasicMatrix<T>::BasicMatrix(std::vector<std::vector<T>> mat): Matrix<T>(mat.size(), mat[0].size()) {
        this->m_data = new T[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++)
            this->m_data[i] = mat[i / this->getCol()][i % this->getCol()];
    }

    template<class T>
    BasicMatrix<T>::BasicMatrix(const BasicMatrix<T> & right): Matrix<T>(right.getRow(), right.getCol()) {
        this->m_data = new T[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++)
            this->m_data[i] = right.m_data[i];
    }

    template<class T>
    BasicMatrix<T>& BasicMatrix<T>::operator=(const BasicMatrix<T>& right){
        if (this == &right) return (*this);
        this->setRow(right.getRow());
        this->setCol(right.getCol());
        this->setSize(right.getSize());
        delete [] m_data;
        m_data = new T[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++)
            this->m_data[i] = right.m_data[i];
        return (*this);
    }

    template<class T>
    T BasicMatrix<T>::getByIndex(int _row, int _col) const{
        return m_data[_row*this->col+_col];  
    }

    template<class T>
    void BasicMatrix<T>::setByIndex(int _row, int _col, T val) {
        m_data[_row*this->col+_col] = val;
    }

    template<class T>
    void BasicMatrix<T>::add(const BasicMatrix<T> &right) {
        if (this->getRow() != right.getRow() || this->getCol() != right.getCol())
            throw ex::MismatchedSizeException(this->getRow(), this->getCol(), right.getRow(), right.getCol(),
                                              "matrix addition");
        for (size_t i = 0; i < this->getSize(); i++) {
            this->m_data[i] = this->m_data[i] + right.m_data[i];
        }
    }

    template<class T>
    void BasicMatrix<T>::add(const SparseMatrix<T> &) {
    }

    template<class T>
    void BasicMatrix<T>::subtract(const BasicMatrix<T> &) {
    }

    template<class T>
    void BasicMatrix<T>::subtract(const SparseMatrix<T> &) {
    }

    template<class T>
    void BasicMatrix<T>::scalarMultiply(T) {
    }

    template<class T>
    void BasicMatrix<T>::scalarDivide(T) {
    }

    template<class T>
    void BasicMatrix<T>::dotProduct(const BasicMatrix<T> &right) {
        if (this->row != right.row) {
            throw ex::MismatchedSizeException(*this, right, "matrix dot product");
        } else if (this->col != 1 && this->col != right.col) {
            throw ex::MismatchedSizeException(*this, right, "matrix dot product");
        }

        int r = right.row;
        int c = right.col;
        if (this->col == 1) {
            BasicMatrix<T> mat(r, c);
            for (int i = 0; i < r; i++) {
                for (int j = 0; j < c; j++) {
                    mat.setByIndex(i, j, this->getByIndex(i,0) * (right.getByIndex(i,j)));
                }
            }
            *this = mat;
        } else {
            this->col = c;
            for (int i = 0; i < r; i++) {
                for (int j = 0; j < c; j++) {
                    this->setByIndex(i, j, this->getByIndex(i,j) * (right.getByIndex(i,j)));
                }
            }
        }
    }

    template<class T>
    void BasicMatrix<T>::dotProduct(const SparseMatrix<T> &) {

    }

    template<class T>
    void BasicMatrix<T>::crossProduct(const BasicMatrix<T> &right) {
        if (this->col != right.row) {
            throw ex::MismatchedSizeException(*this, right, "matrix cross product");
        }
        int r = this->row;
        int c = right.col;
        BasicMatrix<T> mat(r,c);
        for (int k = 0; k < right.row; k++) {
            for (int i = 0; i < r; i++) {
                T temp = this->getByIndex(i, k);
                for (int j = 0; j < c; j++) {
                    mat.setByIndex(i, j, mat.getByIndex(i, j) + temp * right.getByIndex(k, j));
                }
            }
        }
        *this = mat;
    }

    template<class T>
    void BasicMatrix<T>::crossProduct(const SparseMatrix<T> &) {

    }

    template<class T>
    void BasicMatrix<T>::transpose() {

    }

    template<class T>
    void BasicMatrix<T>::inverse() {

    }

    template<class T>
    void BasicMatrix<T>::conjugate() {

    }

    template<class T>
    T BasicMatrix<T>::getMax() {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getMin() {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getSum() {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getAvg() {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getEigenvalue() {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getRowMin(int) {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getColMin(int) {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getRowMax(int) {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getColMax(int) {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getRowSum(int) {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getColSum(int) {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getRowAvg(int) {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getColAvg(int) {
        return m_data[0];
    }

    template<class T>
    Matrix<T> &BasicMatrix<T>::getEigenvector() {
    }

    template<class T>
    T BasicMatrix<T>::getTrace() {
        return m_data[0];
    }

    template<class T>
    T BasicMatrix<T>::getDeterminant() {
        return m_data[0];
    }

    template<class T>
    void BasicMatrix<T>::reshape(int row, int col) {
        long _size = row * col;
        if (this->getSize() == _size) {
            this->setRow(row);
            this->setCol(col);
        } else throw ex::MismatchedSizeException(this->getRow(), this->getCol(), row, col, "matrix reshaping");
    }

    template<class T>
    void BasicMatrix<T>::sliceRow(int row1, int row2) {

    }

    template<class T>
    void BasicMatrix<T>::sliceCol(int col1, int col2) {

    }

    template<class T>
    void BasicMatrix<T>::slice(int row1, int row2, int col1, int col2) {

    }

    template<class T>
    Matrix<T> &BasicMatrix<T>::convolve(BasicMatrix<T> &) {
    }

    template<class T>
    Matrix<T> &BasicMatrix<T>::convolve(SparseMatrix<T> &) {
    }

    template<class T>
    void BasicMatrix<T>::exponent(int exp) {
        
    }

    template<class T>
    void BasicMatrix<T>::show() {
        using namespace std;
        cout << "Basic Matrix:" << endl;
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++)
            cout << getByIndex(i, j) << " ";
            cout << endl;
        }
    }

    template<class T>
    BasicMatrix<T> BasicMatrix<T>::operator+(const BasicMatrix<T> &) {
        return BasicMatrix<T>(0, 0);
    }

    template<class T>
    BasicMatrix<T> BasicMatrix<T>::operator+(const SparseMatrix<T> &) {
        return BasicMatrix<T>(0, 0);
    }

    template<class T>
    BasicMatrix<T> BasicMatrix<T>::operator-(const BasicMatrix<T> &) {
        return BasicMatrix<T>(0, 0);
    }

    template<class T>
    BasicMatrix<T> BasicMatrix<T>::operator-(const SparseMatrix<T> &) {
        return BasicMatrix<T>(0, 0);
    }

    template<class T>
    template<typename P>
    BasicMatrix<T> BasicMatrix<T>::operator*(P){
    }
}

#endif