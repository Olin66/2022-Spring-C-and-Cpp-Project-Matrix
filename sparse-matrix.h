#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "matrix.h"
#include "basic-matrix.h"
#include "matrix-ex.h"
#include <cstring>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <set>
#include <map>

namespace mat {
    namespace ex {
        class MatrixException;
        class MismatchedSizeException;
        class DuplicatedTripleException;
        class NotSquareException;
        class NoInverseException;
        class InvalidSizeException;
        class InvalidTripleException;
    }
    template<typename T>
    struct Triple {
        int _row;
        int _col;
        T val;

        Triple(int _row, int _col, T val): _row(_row), _col(_col), val(val) {}

         bool operator <(const Triple<T> & tri) const
            {
                if (this->_row == tri._row && this->_col == tri._col) return false;
                else return true;
            }

        Triple(const Triple<T> & right){
            this->_row = right._row;
            this->_col = right._col;
            this->val = right.val;
        }
    };

    template<typename>
    class Matrix;

    template<typename>
    class BasicMatrix;

    template<class T>
    class SparseMatrix : public Matrix<T> {
    private:
        long size;
        std::map<int, Triple<T>*> tri_map;
    public:
        SparseMatrix(int, int);

        // SparseMatrix(const cv::Mat &mat);

        SparseMatrix(std::vector<std::vector<T>>);

        SparseMatrix(int, int, T*);

        SparseMatrix(int, int, std::vector<Triple<T>>);

        SparseMatrix(int, int, std::map<int, Triple<T>*>);

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

        void reverse();

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

        Matrix<T> &convolve(BasicMatrix<T> &, int stride = 1, int padding = 0);

        Matrix<T> &convolve(SparseMatrix<T> &, int stride = 1, int padding = 0);

        void exponent(int exp);

        std::map<int , Triple<T>*> getTriples() const {
            return this->tri_map;
        }

        void show();
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
                    setByIndex(i, j, mat[i][j]);
                }
            }
        }
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(int row, int col, T* _data): Matrix<T>(row, col){
        for (size_t i = 0; i < this->getSize(); i++)
        {
            if (_data[i] != 0) {
                setByIndex(i/this->getCol(), i%this->getCol(), _data[i]);
            }
        }
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(int row, int col, std::vector<Triple<T>> mat): Matrix<T>(row, col){
        std::set<Triple<T>> temp(mat.begin(), mat.end());
        if (mat.size() != temp.size()) {
            throw ex::DuplicatedTripleException();
        }
        for (size_t i = 0;i < mat.size();i++)
        {
            Triple<T>* t = new Triple<T>(mat[i]);
            long index = t->_row * this->getCol() + t->_col;
            if (index >= this->getSize() || t->_row < 0 || t->_row >= row || t->_col < 0 || t->_col >= col)
                {
                    for (auto it = this->tri_map.begin();it != this->tri_map.end();it++)
                {
                    if (it->second != nullptr) delete it->second;
                }
                    throw ex::InvalidTripleException(t->_row, t->_col);}
            tri_map[index] = t;
        }
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(int row, int col, std::map<int , Triple<T>*> _map): Matrix<T>(row, col){
        for (auto it = _map.begin();it != _map.end();it++)
        {
            Triple<T> temp = it;
            Triple<T>* t = new Triple<T>(temp->_row, temp->_col, temp->val);
            long index = t->_row * this->getCol() + t->_col;
            if (index >= this->getSize() || t->_row < 0 || t->_row >= row || t->_col < 0 || t->_col >= col)
               { 
                for (auto it = this->tri_map.begin();it != this->tri_map.end();it++)
                {
                    if (it->second != nullptr) delete it->second;
                }
                
                throw ex::InvalidTripleException(t->_row, t->_col);
            }
        }
        
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(const SparseMatrix<T> & right): Matrix<T>(right.getRow(), right.getCol()){
        for (auto it = right.tri_map.begin();it != right.tri_map.end();it++)
        {
            auto tri = it->second;
            this->setByIndex(tri->_row, tri->_col, tri->val);
        }
    }

    template<class T>
    SparseMatrix<T>& SparseMatrix<T>::operator=(const SparseMatrix<T> & right){
        this->setSize(right.getSize());
        this->setRow(right.getRow());
        this->setCol(right.getCol());
        for (auto it = this->tri_map.begin(); it != tri_map.end();it++)
        {
            if (it->second != nullptr) delete it->second;
        }
        this->tri_map.clear();
        for (auto it = right.tri_map.begin(); it != right.tri_map.end();it++)
        {
            auto tri = it->second;
            this->setByIndex(tri->_row, tri->_col, tri->val);
        }
    }

    template <class T>
    T SparseMatrix<T>::getByIndex(int _row, int _col) const {
        int index = _row * this->col + _col;
        if (tri_map.find(index)->second == nullptr)
            return 0;
        else 
            return tri_map.find(index)->second->val;
    }

    template <class T>
    void SparseMatrix<T>::setByIndex(int _row, int _col, T val) {
        long index = _row * this->col + _col;
        if (tri_map[index] == nullptr) {
            Triple<T>* t = new Triple<T>(_row, _col, val);
            tri_map[index] = t;
        } else {
            tri_map[index]->val = val;
        }
    }

    template<class T>
    void SparseMatrix<T>::add(const BasicMatrix<T> &) {
    }

    template<class T>
    void SparseMatrix<T>::add(const SparseMatrix<T> &right) {
        SparseMatrix<T> mat(this->row, this->col);
        for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
            for (auto j = right.tri_map.begin(); j != right.tri_map.end(); j++) {
                Triple<T>* LP = i->second;
                Triple<T>* RP = j->second;
                if (LP->_row == RP->_row && LP->_row == RP->_row) {
                    mat.setByIndex(LP->_row, LP->_col, LP->val + RP->val);
                }
            }
        }
        *this = mat;
    }

    template<class T>
    void SparseMatrix<T>::subtract(const BasicMatrix<T> &) {
    }

    template<class T>
    void SparseMatrix<T>::subtract(const SparseMatrix<T> &right) {
        SparseMatrix<T> mat(this->row, this->col);
        for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
            for (auto j = right.tri_map.begin(); j != right.tri_map.end(); j++) {
                Triple<T>* LP = i->second;
                Triple<T>* RP = j->second;
                if (LP->_row == RP->_row && LP->_row == RP->_row) {
                    mat.setByIndex(LP->_row, LP->_col, LP->val - RP->val);
                }
            }
        }
        *this = mat;
    }

    template<class T>
    void SparseMatrix<T>::scalarMultiply(T x) {
        for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
            Triple<T> *tri = i->second;
            this->setByIndex(tri->_row, tri->_col, tri->val * x);
        }
    }

    template<class T>
    void SparseMatrix<T>::scalarDivide(T x) {
        for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
            Triple<T> *tri = i->second;
            this->setByIndex(tri->_row, tri->_col, tri->val / x);
        }
    }

    template<class T>
    void SparseMatrix<T>::dotProduct(const BasicMatrix<T> &) {

    }

    template<class T>
    void SparseMatrix<T>::dotProduct(const SparseMatrix<T> &right) {
        if (this->row != right.row) {
            throw ex::MismatchedSizeException(*this, right, "matrix dot product");
        } else if (this->col != 1 && this->col != right.col) {
            throw ex::MismatchedSizeException(*this, right, "matrix dot product");
        }

        int r = right.row;
        int c = right.col;
        if (this->col == 1) {
            SparseMatrix<T> mat(r, 1);
            for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
                for (auto j = right.tri_map.begin(); j != right.tri_map.end(); j++) {
                    Triple<T>* LP = i->second;
                    Triple<T>* RP = j->second;
                    if (LP->_row == RP->_row && LP->_row == RP->_row) {
                        mat.setByIndex(RP->_row, RP->_col, LP->val * RP->val);
                    }
                }
            }
            *this = mat;
        } else {
            SparseMatrix<T> mat(r, c);
            for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
                for (auto j = right.tri_map.begin(); j != right.tri_map.end(); j++) {
                    Triple<T>* LP = i->second;
                    Triple<T>* RP = j->second;
                    if (LP->_row == RP->_row &&
                        LP->_col == RP->_col) {
                        mat.setByIndex(RP->_row, RP->_col, LP->val * RP->val);
                    }
                }
            }
            *this = mat;
        }
    }

    template<class T>
    void SparseMatrix<T>::crossProduct(const BasicMatrix<T> &) {

    }

    template<class T>
    void SparseMatrix<T>::crossProduct(const SparseMatrix<T> &right) {
        if (this->col != right.col) {
            throw ex::MismatchedSizeException(*this, right, "matrix cross product");
        }
        int r = this->row;
        int c = right.col;
        SparseMatrix<T> mat(r, c);
        for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
            Triple<T> *l_point = i->second;
            for (auto j = right.tri_map.begin(); j != right.tri_map.end(); j++) {
                Triple<T> *r_point = j->second;
                if (l_point->_col == r_point->_row) {
                    mat.setByIndex(l_point->_row, r_point->_col, mat.getByIndex(l_point->_row, r_point->_col) + l_point->val * r_point->val);
                }
            }
        }
        *this = mat;
    }

    template<class T>
    void SparseMatrix<T>::transpose() {

    }

    template<class T>
    void SparseMatrix<T>::inverse() {

    }

    template<class T>
    void SparseMatrix<T>::reverse() {
        
    }

    template<class T>
    void SparseMatrix<T>::conjugate() {

    }

    template<class T>
    T SparseMatrix<T>::getMax() {
        T max=0;
        // for (int i = 0;i< this->tri_map.size(); i++)
        // {
        //     if (max<this->tri_map[i].val)
        //     {
        //         max=this->tri_map[i].val;
        //     }
            
        // }
        return max;
    }

    template<class T>
    T SparseMatrix<T>::getMin() {
        T min=0;
        // for (int i = 0; i < this->tri_map.size(); i++)
        // {
        //     if (min>this->tri_map[i].val)
        //     {
        //         min=this->tri_map[i].val;
        //     }
            
        // }
        return min;
    }

    template<class T>
    T SparseMatrix<T>::getSum() {
        T sum=0;
        // for (int i = 0; i < this->tri_map.size(); i++)
        // {
        //     sum+=this->tri_map[i].val;
        // }
        return sum;
        
    }

    template<class T>
    T SparseMatrix<T>::getAvg() {
        return this->getSum()/(this->col*this->row);
    }

    template<class T>
    T SparseMatrix<T>::getRowMax(int row) {
        T max=0;
        // for (int i = 0; i < this->tri_map.size(); i++)
        // {
        //     if (max<this->tri_map[i].val&&this->tri_map[i]._row==row)
        //     {
        //         max=this->tri_map[i].val;
        //     }
            
        // }
        return max;
    }

    template<class T>
    T SparseMatrix<T>::getColMax(int col) {
        T max=0;
        // for (int i = 0; i < this->tri_map.size(); i++)
        // {
        //     if (max<this->tri_map[i].val&&this->tri_map[i]._col==col)
        //     {
        //         max=this->tri_map[i].val;
        //     }
        // }
        return max;
    }

    template<class T>
    T SparseMatrix<T>::getRowMin(int row) {
        T min=0;
        // for (int i = 0; i < this->tri_map.size(); i++)
        // {
        //     if (min>this->tri_map[i].val&&this->tri_map[i]._row==row)
        //     {
        //         min=this->tri_map[i].val;
        //     }
            
        // }
        return min;

    }

    template<class T>
    T SparseMatrix<T>::getColMin(int col) {
         T min=0;
        // for (int i = 0; i < this->tri_map.size(); i++)
        // {
        //     if (min>this->tri_map[i].val&&this->tri_map[i]._col==col)
        //     {
        //         min=this->tri_map[i].val;
        //     }
            
        // }
        return min;
    }

    template<class T>
    T SparseMatrix<T>::getRowSum(int row) {
        T sum=0;
        // for (int i = 0; i < this->tri_map.size(); i++)
        // {
        //     if (this->tri_map[i]._row==row)
        //     {
        //         sum+=this->tri_map[i].val;
        //     }
            
        // }
        return sum;
    }

    template<class T>
    T SparseMatrix<T>::getColSum(int col) {
        T sum=0;
        // for (int i = 0; i < this->tri_map.size(); i++)
        // {
        //     if (this->tri_map[i]._col==col)
        //     {
        //         sum+=this->tri_map[i].val;
        //     }
            
        // }
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
        long _size = row * col;
    if (this->getSize() == _size) {
        this->setRow(col);
        this->setCol(row);
        for (auto i = tri_map.begin(); i != tri_map.end(); i++)
        {
            Triple<T> *tri = i->second;
            int temp = tri->_col;
            tri->_col = tri->_row;
            tri->_row = temp;
        }
    } else
        throw ex::MismatchedSizeException(this->getRow(), this->getCol(), row, col, "matrix reshaping");
        
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
    Matrix<T> &SparseMatrix<T>::convolve(BasicMatrix<T> &, int stride, int padding) {
    }

    template<class T>
    Matrix<T> &SparseMatrix<T>::convolve(SparseMatrix<T> &, int stride, int padding) {
    }

    template<class T>
    void SparseMatrix<T>::exponent(int exp){
        
    }

    template<class T>
    void SparseMatrix<T>::show(){
        using namespace std;
        T mat[this->getRow()][this->getCol()];
        memset(mat, 0, sizeof(mat));
        // cout<<tri_map[3]->_row<<endl;
        // cout<<tri_map[3]->_col<<endl;
        // cout<<tri_map[3]->val<<endl;
        for (auto i = tri_map.begin(); i != tri_map.end(); i++)
        {
            auto tri = i->second;
            mat[tri->_row][tri->_col] = tri->val;
        }
        for (size_t i = 0; i < this->getRow(); i++)
        {
            for (size_t j = 0;j < this->getCol(); j++)
            {
                cout<<mat[i][j]<<" ";
            }
            cout<<endl;
        }
    }
}

#endif