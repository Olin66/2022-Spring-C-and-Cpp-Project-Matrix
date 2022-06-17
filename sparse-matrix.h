#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "matrix-ex.h"
#include "matrix.h"
#include "basic-matrix.h"
#include <cstring>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <set>

namespace mat {

    template<typename T>
    struct Triple {
        long _row;
        long _col;
        T val;
    };

    template<typename>
    class Matrix;

    template<typename>
    class BasicMatrix;

    template<class T>
    class SparseMatrix : public Matrix<T> {
    private:
        std::set<Triple<T>> triples;
    public:
        SparseMatrix(int, int);

        // SparseMatrix(const cv::Mat &mat);

        SparseMatrix(std::vector<std::vector<T>>);

        SparseMatrix(int, int, T*);

        SparseMatrix(int, int, std::vector<Triple<T>>);

        SparseMatrix(const SparseMatrix<T> &);

        SparseMatrix<T>& operator=(const SparseMatrix<T> &);

        void add(const BasicMatrix<T> &);

        void add(const SparseMatrix<T> &);

        void subtract(const BasicMatrix<T> &);

        void subtract(const SparseMatrix<T> &);

        void scalarMultiply(T);

        void scalarDivide(T);

        void scalarDivide(std::complex<long double>);

        void dotProduct(const BasicMatrix<T> &);

        void dotProduct(const SparseMatrix<T> &);

        void crossProduct(const BasicMatrix<T> &);

        void crossProduct(const SparseMatrix<T> &);

        void transpose();

        void inverse();

        void conjugate();

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

        T getByIndex(int _row, int _col) const;

        void setByIndex(int _row, int _col, T val);

        T getTrace();

        T getDeterminant();

        void reshape(int row, int col);

        void sliceRow(int row1, int row2);

        void sliceCol(int col1, int col2);

        void slice(int row1, int row2, int col1, int col2);

        Matrix<T> &convolve(BasicMatrix<T> &);

        Matrix<T> &convolve(SparseMatrix<T> &);

        void exponent(int exp);
    };

    template<class T>
    SparseMatrix<T>::SparseMatrix(int row, int col): Matrix<T>(row, col) {}

    // template<class T>
    // SparseMatrix<T>::SparseMatrix(const cv::Mat &mat) {

    // }

    template<class T>
    SparseMatrix<T>::SparseMatrix(std::vector<std::vector<T>> mat): Matrix<T>(mat.size(), mat[0].size()) {
        for (size_t i = 0; i < mat.size(); i++)
        {
            for (size_t j = 0; j < mat[i].size(); j++)
            {
                if (mat[i][j] != 0) {
                    Triple<T> triple(i, j, mat[i][j]);
                    triples.insert(triple);
                }
            }
            
        }
        
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(int row, int col, T* _data): Matrix<T>(row, col){
        for (size_t i = 0; i < this->getSize(); i++)
        {
            if (_data[i] != 0) {
                Triple<T> triple(i/this->getCol(), i%this->getCol(), _data[i]);
                triples.insert(triple);
            }
        }
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(int row, int col, std::vector<Triple<T>> mat): Matrix<T>(row, col){
        
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(const SparseMatrix<T> & right): Matrix<T>(right.getRow(), right.getCol()){
        this->triples(right.triples);
    }

    template<class T>
    SparseMatrix<T>& SparseMatrix<T>::operator=(const SparseMatrix<T> & right){
        this->setSize(right.getSize());
        this->setRow(right.getRow());
        this->setCol(right.getCol());
        this->triples(right.triples);
    }
    template <class T>
T SparseMatrix<T>::getByIndex(int _row, int _col) const {
    T target=0;
    for (int i = 0; i < this->triples.size(); i++)
    {
        if (this->triples[i]._col==_col &&this->triples[i]._row==_row)
        {
            return target=this->triples[i].val;
        }
        
    }
    
    
    
    return target;
}

template <class T>
void SparseMatrix<T>::setByIndex(int _row, int _col, T val) {
    for (int i = 0; i < this->triples.size(); i++)
    {
        if (this->triples[i]._row==_row&&this->triples[i]._col==_col)
        {
            this->triples[i].val=val;
            return;
        }
    }
    Triple<T>triple(_row,_col,val);
    this->triples.push_back(triple);
    return;
}
    template<class T>
    void SparseMatrix<T>::add(const BasicMatrix<T> &) {
    }

    template<class T>
    void SparseMatrix<T>::add(const SparseMatrix<T> &) {
    }

    template<class T>
    void SparseMatrix<T>::subtract(const BasicMatrix<T> &) {
    }

    template<class T>
    void SparseMatrix<T>::subtract(const SparseMatrix<T> &) {
    }

    template<class T>
    void SparseMatrix<T>::scalarMultiply(T) {
    }

    template<class T>
    void SparseMatrix<T>::scalarDivide(T) {
    }

    template<class T>
    void SparseMatrix<T>::scalarDivide(std::complex<long double>) {
    }

    template<class T>
    void SparseMatrix<T>::dotProduct(const BasicMatrix<T> &) {

    }

    template<class T>
    void SparseMatrix<T>::dotProduct(const SparseMatrix<T> &) {

    }

    template<class T>
    void SparseMatrix<T>::crossProduct(const BasicMatrix<T> &) {

    }

    template<class T>
    void SparseMatrix<T>::crossProduct(const SparseMatrix<T> &) {

    }

    template<class T>
    void SparseMatrix<T>::transpose() {

    }

    template<class T>
    void SparseMatrix<T>::inverse() {

    }

    template<class T>
    void SparseMatrix<T>::conjugate() {

    }

    template<class T>
    T SparseMatrix<T>::getMax() {
        T max=0;
        for (int i = 0;i< this->triples.size(); i++)
        {
            if (max<this->triples[i].val)
            {
                max=this->triples[i].val;
            }
            
        }
        return max;
        
        
        
    }

    template<class T>
    T SparseMatrix<T>::getMin() {
        T min=0;
        for (int i = 0; i < this->triples.size(); i++)
        {
            if (min>this->triples[i].val)
            {
                min=this->triples[i].val;
            }
            
        }
        return min;
        
        
    }

    template<class T>
    T SparseMatrix<T>::getSum() {
        T sum=0;
        for (int i = 0; i < this->triples.size(); i++)
        {
            sum+=this->triples[i].val;
        }
        return sum;
        
    }

    template<class T>
    T SparseMatrix<T>::getAvg() {
        return this->getSum/(this->col*this->row);
    }

    template<class T>
    T SparseMatrix<T>::getRowMax(int row) {
        T max=0;
        for (int i = 0; i < this->triples.size(); i++)
        {
            if (max<this->triples[i].val&&this->triples[i]._row==row)
            {
                max=this->triples[i].val;
            }
            
        }
        return max;
    }

    template<class T>
    T SparseMatrix<T>::getColMax(int col) {
        T max=0;
        for (int i = 0; i < this->triples.size(); i++)
        {
            if (max<this->triples[i].val&&this->triples[i]._col==col)
            {
                max=this->triples[i].val;
            }
            
        }
        return max;
    }

    template<class T>
    T SparseMatrix<T>::getRowMin(int row) {
        T min=0;
        for (int i = 0; i < this->triples.size(); i++)
        {
            if (min>this->triples[i].val&&this->triples[i]._row==row)
            {
                min=this->triples[i].val;
            }
            
        }
        return min;

    }

    template<class T>
    T SparseMatrix<T>::getColMin(int col) {
         T min=0;
        for (int i = 0; i < this->triples.size(); i++)
        {
            if (min>this->triples[i].val&&this->triples[i]._col==col)
            {
                min=this->triples[i].val;
            }
            
        }
        return min;
    }

    template<class T>
    T SparseMatrix<T>::getRowSum(int row) {
        T sum=0;
        for (int i = 0; i < this->triples.size(); i++)
        {
            if (this->triples[i]._row==row)
            {
                sum+=this->triples[i].val;
            }
            
        }
        return sum;
    }

    template<class T>
    T SparseMatrix<T>::getColSum(int col) {
        T sum=0;
        for (int i = 0; i < this->triples.size(); i++)
        {
            if (this->triples[i]._col==col)
            {
                sum+=this->triples[i].val;
            }
            
        }
        return sum;
    }

    template<class T>
    T SparseMatrix<T>::getRowAvg(int row) {
        return this->getRowSum(row)/this->col;
    }

    template<class T>
    T SparseMatrix<T>::getColAvg(int col) {
        return this->getColSum(col)/this->row;
    }

    template<class T>
    T SparseMatrix<T>::getEigenvalue() {
    }

    template<class T>
    Matrix<T> &SparseMatrix<T>::getEigenvector() {
    }

    template<class T>
    T SparseMatrix<T>::getTrace() {
    }

    template<class T>
    T SparseMatrix<T>::getDeterminant() {
    }

    template<class T>
    void SparseMatrix<T>::reshape(int row, int col) {

    }

    template<class T>
    void SparseMatrix<T>::sliceRow(int row1, int row2) {

    }

    template<class T>
    void SparseMatrix<T>::sliceCol(int col1, int col2) {

    }

    template<class T>
    void SparseMatrix<T>::slice(int row1, int row2, int col1, int col2) {

    }

    template<class T>
    Matrix<T> &SparseMatrix<T>::convolve(BasicMatrix<T> &) {
    }

    template<class T>
    Matrix<T> &SparseMatrix<T>::convolve(SparseMatrix<T> &) {
    }

    template<class T>
    void SparseMatrix<T>::exponent(int exp){
        
    }

}

#endif