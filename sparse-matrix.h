#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "matrix.h"
#include "basic-matrix.h"
#include "matrix-ex.h"
#include <opencv2/opencv.hpp>
#include <cstring>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <set>
#include <map>
using namespace std;
using namespace cv;
namespace mat {
    namespace ex {
        class MatrixException;

        class MismatchedSizeException;

        class DuplicatedTripleException;

        class NotSquareException;

        class NoInverseException;

        class InvalidSizeException;

        class InvalidTripleException;

        class InvalidChannelDepth;
    }
    template<typename T>
    struct Triple {
        int _row;
        int _col;
        T val;

        Triple(int _row, int _col, T val) : _row(_row), _col(_col), val(val) {}

        bool operator<(const Triple<T> &tri) const {
            if (this->_row == tri._row && this->_col == tri._col) return false;
            else return true;
        }

        Triple(const Triple<T> &right) {
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
        map<int, Triple<T> *> tri_map;
    public:
        SparseMatrix(int, int);

        SparseMatrix(const cv::Mat &mat);

        SparseMatrix(vector<vector<T>>);

        SparseMatrix(int, int, T *);

        SparseMatrix(int, int, vector<Triple<T>>);

        SparseMatrix(int, int, map<int, Triple<T> *>);

        SparseMatrix(const SparseMatrix<T> &);

        SparseMatrix<T> &operator=(const SparseMatrix<T> &);

        ~SparseMatrix<T>();

        void add(const SparseMatrix<T> &);

        void subtract(const SparseMatrix<T> &);

        void scalarMultiply(T);

        void scalarDivide(T);

        void dotProduct(const SparseMatrix<T> &);

        void crossProduct(const SparseMatrix<T> &);

        bool getEigen(SparseMatrix<T> &eigenvector, SparseMatrix<T> &eigenvalue, double error, double iterator);//雅克比法计算特征值和特征向量

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

        T getEigenvalue(int LoopNumber,  SparseMatrix<T> &result);

        void Gaussian_Eliminate(SparseMatrix<T>&ans,SparseMatrix<T>&eigenmatirx);

        SparseMatrix<T> &getEigenvector(SparseMatrix<T> &eigenvector,const T lamda);

        T getByIndex(int _row, int _col) const;

        void setByIndex(int _row, int _col, T val);

        T getTrace();

        T getDeterminant();

        void reshape(int row, int col);

        void sliceRow(int row1, int row2);

        void sliceCol(int col1, int col2);

        void slice(int row1, int row2, int col1, int col2);

        SparseMatrix<T> *convolve(SparseMatrix<T> &, int stride = 1, int padding = 0);

        void exponent(int exp);

        map<int, Triple<T> *> getTriples() const {
            return this->tri_map;
        }

        void QR(SparseMatrix<T>&Q,SparseMatrix<T>&R);

        void show();

        cv::Mat* getCvMat();
    };

    template<class T>
    Mat* SparseMatrix<T>::getCvMat(){
        bool flags[this->getRow()][this->getCol()];
        memset(flags, false, sizeof(flags));
        Mat* mat = new Mat(this->getRow(), this->getCol(), CV_8UC1);
        for (auto it = tri_map.begin(); it != tri_map.end(); it++)
        {
            auto tri = it->second;
            double re = real(getByIndex(tri->_row, tri->_col));
            mat->at<uchar>(tri->_row, tri->_col) = re;
            flags[tri->_row][tri->_col];
        }
        for (size_t i = 0; i < this->getRow(); i++)
        {
            for (size_t j = 0;j < this->getCol();j++)
            {
                if (!flags[i][j]) mat->at<uchar>(i, j) = 0;
            }
        }
        return mat;
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(int row, int col): Matrix<T>(row, col) {}

    template<class T>
    SparseMatrix<T>::SparseMatrix(const cv::Mat &mat): Matrix<T>(mat) {
        if (mat.channels() != 1)
            throw ex::InvalidChannelDepth(mat.channels());
        for (size_t i = 0; i < this->getRow(); i++)
        {
            for (size_t j = 0; j < this->getCol(); j++)
            {
                if ((T)mat.at<uchar>(i,j) != 0)
                    setByIndex(i, j, (T)mat.at<uchar>(i,j));
            }
        }
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(vector<vector<T>> mat): Matrix<T>(mat.size(), mat[0].size()) {
        for (size_t i = 0; i < mat.size(); i++) {
            for (size_t j = 0; j < mat[i].size(); j++) {
                if (mat[i][j] != 0) {
                    setByIndex(i, j, mat[i][j]);
                }
            }
        }
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(int row, int col, T *_data): Matrix<T>(row, col) {
        for (size_t i = 0; i < this->getSize(); i++) {
            if (_data[i] != 0) {
                setByIndex(i / this->getCol(), i % this->getCol(), _data[i]);
            }
        }
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(int row, int col, vector<Triple<T>> mat): Matrix<T>(row, col) {
        set<Triple<T>> temp(mat.begin(), mat.end());
        if (mat.size() != temp.size()) {
            throw ex::DuplicatedTripleException();
        }
        for (size_t i = 0; i < mat.size(); i++) {
            Triple<T> *t = new Triple<T>(mat[i]);
            long index = t->_row * this->getCol() + t->_col;
            if (index >= this->getSize() || t->_row < 0 || t->_row >= row || t->_col < 0 || t->_col >= col) {
                for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
                    if (it->second != nullptr) delete it->second;
                }
                throw ex::InvalidTripleException(t->_row, t->_col);
            }
            tri_map[index] = t;
        }
    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(int row, int col, map<int, Triple<T> *> _map): Matrix<T>(row, col) {
        for (auto it = _map.begin(); it != _map.end(); it++) {
            Triple<T> temp = it;
            Triple<T> *t = new Triple<T>(temp->_row, temp->_col, temp->val);
            long index = t->_row * this->getCol() + t->_col;
            if (index >= this->getSize() || t->_row < 0 || t->_row >= row || t->_col < 0 || t->_col >= col) {
                for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
                    if (it->second != nullptr) delete it->second;
                }

                throw ex::InvalidTripleException(t->_row, t->_col);
            }
        }

    }

    template<class T>
    SparseMatrix<T>::SparseMatrix(const SparseMatrix<T> &right): Matrix<T>(right.getRow(), right.getCol()) {
        for (auto it = right.tri_map.begin(); it != right.tri_map.end(); it++) {
            auto tri = it->second;
            this->setByIndex(tri->_row, tri->_col, tri->val);
        }
    }

    template<class T>
    SparseMatrix<T> &SparseMatrix<T>::operator=(const SparseMatrix<T> &right) {
        this->setRow(right.getRow());
        this->setCol(right.getCol());
        this->setSize(right.getSize());
        for (auto it = this->tri_map.begin(); it != tri_map.end(); it++) {
            if (it->second != nullptr) delete it->second;
        }
        this->tri_map.clear();
        for (auto it = right.tri_map.begin(); it != right.tri_map.end(); it++) {
            auto tri = it->second;
            this->setByIndex(tri->_row, tri->_col, tri->val);
        }
    }

    template<class T>
    SparseMatrix<T>::~SparseMatrix<T>() {
        for (auto it = this->tri_map.begin(); it != tri_map.end(); it++) {
            if (it->second != nullptr) delete it->second;
        }
        this->tri_map.clear();
    }

    
    template <class T>
    bool SparseMatrix<T>::getEigen(SparseMatrix<T> &eigenvector, SparseMatrix<T> &eigenvalue, double error, double iterator) {
        for (int i = 0; i < this->col; i++) {
            eigenvector.setByIndex(i, i, 1.0);
            for (int j = 0; j < this->col; j++) {
                if (i != j) {
                    eigenvector.setByIndex(i, j, 0.0);
                }
            }
        }
        int count = 0;
        while (true) {
            T max = this->getByIndex(0, 1);
            int row = 0, col = 1;
            for (int i = 0; i < this->row; i++) //找到非对角线最大元素
            {
                for (int j = 0; j < this->col; j++) {
                    T d = fabs(this->getByIndex(i, j));
                    if (i != j && d > max) {
                        max = d;
                        row = i;
                        col = j;
                    }
                }
            }
            if (max < error || count > iterator) {
                break;
            }
            count++;
            T pp = this->getByIndex(row, row);
            T pq = this->getByIndex(row, col);
            T qq = this->getByIndex(col, col);
            //设置旋转角度
            T Angle = 0.5 * atan2(-2 * pq, qq - pp);
            T Sin = sin(Angle);
            T Cos = cos(Angle);
            T Sin2 = sin(2 * Angle);
            T Cos2 = cos(2 * Angle);
            this->setByIndex(row, row, pp * Cos * Cos + qq * Sin * Sin + 2 * pq * Cos * Sin);
            this->setByIndex(col, col, pp * Sin * Sin + qq * Cos * Cos - 2 * pq * Cos * Sin);
            this->setByIndex(row, col, (qq - pp) * Sin2 / 2 + pq * Cos2);
            this->setByIndex(col, row, this->getByIndex(row, col));
            for (int i = 0; i < this->col; i++) {
                if (i != col && i != row) {
                    max = this->getByIndex(i, row);
                    this->setByIndex(i, row, this->getByIndex(i, col) * Sin + max * Cos);
                    this->setByIndex(i, col, this->getByIndex(i, col) * Cos - max * Sin);
                }
            }
            for (int j = 0; j < this->col; j++) {
                if (j != col && j != row) {
                    max = this->getByIndex(row, j);
                    this->setByIndex(row, j, this->getByIndex(col, j) * Sin + max * Cos);
                    this->setByIndex(col, j, this->getByIndex(col, j) * Cos - max * Sin);
                }
            }
            //计算特征向量
            for (int i = 0; i < this->col; i++) {
                max = eigenvector.getByIndex(i, row);
                eigenvector.setByIndex(i, row, eigenvector.getByIndex(i, col) * Sin + max * Cos);
                eigenvector.setByIndex(i, col, eigenvector.getByIndex(i, col) * Cos - max * Sin);
            }
        }
        map<T, int> mapEigen;
        for (int i = 0; i < this->col; i++) {
            eigenvalue.setByIndex(i, i, this->getByIndex(i, i));
            mapEigen.insert(make_pair(eigenvalue.getByIndex(i, i), i));
        }

        SparseMatrix<T> temp(this->col, this->col);
        typename map<T, int>::reverse_iterator it = mapEigen.rbegin();
        for (int j = 0; it != mapEigen.rend(), j < this->col; it++, j++) {
            for (int i = 0; i < this->col; i++) {
                temp.setByIndex(i, j, eigenvector.getByIndex(i, it->second));
            }
            eigenvalue.setByIndex(j, j, it->first);
        }

        for (int i = 0; i < this->col; i++) {
            T sum = 0;
            for (int j = 0; j < this->col; j++) {
                sum += temp.getByIndex(j, i);
            }
            if (sum < 0) {
                for (int j = 0; j < this->col; j++) {
                    temp.setByIndex(j, i, -temp.getByIndex(j, i));
                }
            }
        }
        eigenvector = temp;
        return true;
    }

    template<class T>
    T SparseMatrix<T>::getByIndex(int _row, int _col) const {
        int index = _row * this->col + _col;
        if (tri_map.find(index) == tri_map.end())
            return 0;
        else
            return tri_map.find(index)->second->val;
    }

    template<class T>
    void SparseMatrix<T>::setByIndex(int _row, int _col, T val) {
        int index = _row * this->col + _col;
        if (tri_map[index] == nullptr) {
            tri_map[index] = new Triple<T>(_row, _col, val);
        } else {
            tri_map[index]->val = val;
        }
    }

    template<class T>
    void SparseMatrix<T>::add(const SparseMatrix<T> &right) {
        SparseMatrix<T> mat(this->row, this->col);
        for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
            for (auto j = right.tri_map.begin(); j != right.tri_map.end(); j++) {
                Triple<T> *LP = i->second;
                Triple<T> *RP = j->second;
                if (LP->_row == RP->_row && LP->_row == RP->_row) {
                    mat.setByIndex(LP->_row, LP->_col, LP->val + RP->val);
                }
            }
        }
        *this = mat;
    }

    template<class T>
    void SparseMatrix<T>::subtract(const SparseMatrix<T> &right) {
        SparseMatrix<T> mat(this->row, this->col);
        for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
            for (auto j = right.tri_map.begin(); j != right.tri_map.end(); j++) {
                Triple<T> *LP = i->second;
                Triple<T> *RP = j->second;
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
    void SparseMatrix<T>::dotProduct(const SparseMatrix<T> &right) {
        if (this->row != right.row) {
            throw ex::MismatchedSizeException(*this, right, "matrix dot product");
        } else if (this->col != 1 && this->col != right.col) {
            throw ex::MismatchedSizeException(*this, right, "matrix dot product");
        }

        int r = right.getRow();
        int c = right.getCol();
        if (this->col == 1) {
            SparseMatrix<T> mat(r, 1);
            for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
                for (auto j = right.tri_map.begin(); j != right.tri_map.end(); j++) {
                    Triple<T> *LP = i->second;
                    Triple<T> *RP = j->second;
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
                    Triple<T> *LP = i->second;
                    Triple<T> *RP = j->second;
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
    void SparseMatrix<T>::crossProduct(const SparseMatrix<T> &right) {
        if (this->col != right.row) {
            throw ex::MismatchedSizeException(*this, right, "matrix cross product");
        }
        int r = this->row;
        int c = right.col;
        SparseMatrix<T> mat(r, c);
        for (auto i = this->tri_map.begin(); i != this->tri_map.end(); i++) {
            Triple<T> *tri = i->second;
            for (auto j = right.tri_map.begin(); j != right.tri_map.end(); j++) {
                Triple<T> *trj = j->second;
                if (tri->_col == trj->_row) {
                    mat.setByIndex(tri->_row, trj->_col, mat.getByIndex(tri->_row, trj->_col) + tri->val * trj->val);
                }
            }
        }
        *this = mat;
    }

    template <class T>
    void SparseMatrix<T>::transpose() {
        SparseMatrix<T> mat(this->col, this->row);
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                mat.setByIndex(j, i, this->getByIndex(i, j));
            }
        }
        *this = mat;
    }

    template<class T>
    void SparseMatrix<T>::inverse() {
        T det = this->getDeterminant();
        if (this->col != this->col || det == 0) //行列式为0的矩阵不可逆
        {
            throw ex::NoInverseException(*this, "matrix inverse");
        }
        SparseMatrix<T> temp(this->row - 1, this->col - 1);
        SparseMatrix<T> adjoint(this->row, this->col);
        if (this->col == 1) {
            return; //一阶矩阵逆是本身
        }
        for (int i = 0; i < this->col; i++) {
            for (int j = 0; j < this->col; j++) {
                for (int k = 0; k < this->col - 1; k++) {
                    for (int l = 0; l < this->col - 1; l++) {
                        temp.setByIndex(k, l, this->getByIndex(k >= i ? k + 1 : k, l >= j ? l + 1 : l));
                    }
                }
                adjoint.setByIndex(j, i, temp.getDeterminant()); //进行了一次转置。
                if ((i + j) % 2 == 1) {
                    adjoint.setByIndex(j, i, -adjoint.getByIndex(j, i)); //将每一个元素的代数余子式加上正负号生成伴随矩阵
                }
            }
        }
        for (int i = 0; i < this->col; i++) {
            for (int j = 0; j < this->row; j++) {
                this->setByIndex(i, j, adjoint.getByIndex(i, j) / det);
            }
        }
    }


    template<class T>
    void SparseMatrix<T>::reverse() {
        int _size = this->getSize();
        map<int, Triple<T> *> _map;
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            int index = it->first;
            Triple<T> *tri = it->second;
            _map[_size - index - 1] = new Triple<T>(this->row - tri->_row - 1, this->col - tri->_col - 1, tri->val);
        }
        this->tri_map = _map;
    }

    template<class T>
    void SparseMatrix<T>::conjugate() {
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                if (imag(this->getByIndex(i, j)) != 0) {
                    T temp = real(this->getByIndex(i, j)) - imag(this->getByIndex(i, j));
                    this->setByIndex(i, j, temp);
                }
            }
        }
    }

    template<class T>
    T SparseMatrix<T>::getMax() {
        T max = 0;
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            if (max < it->second->val) {
                max = it->second->val;
            }

        }
        return max;
    }

    template<class T>
    T SparseMatrix<T>::getMin() {
        T min = 0;
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            if (min > it->second->val) {
                min = it->second->val;
            }

        }
        return min;
    }

    template<class T>
    T SparseMatrix<T>::getSum() {
        T sum = 0;
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            sum += (it->second->val);
        }
        return sum;

    }

    template<class T>
    T SparseMatrix<T>::getAvg() {
        return this->getSum() / (this->col * this->row);
    }

    template<class T>
    T SparseMatrix<T>::getRowMax(int row) {
        T max = 0;
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            if (it->second->_row == row && max < it->second->val) {
                max = it->second->val;
            }

        }
        return max;
    }

    template<class T>
    T SparseMatrix<T>::getColMax(int col) {
        T max = 0;
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            if (it->second->_col == col && max < it->second->val) {
                max = it->second->val;
            }

        }
        return max;
    }

    template<class T>
    T SparseMatrix<T>::getRowMin(int row) {
        T min = 0;
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            if (it->second->_row == row && min > it->second->val) {
                min = it->second->val;
            }

        }
        return min;

    }

    template<class T>
    T SparseMatrix<T>::getColMin(int col) {
        T min = 0;
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            if (it->second->_col == col && min > it->second->val) {
                min = it->second->val;
            }

        }
        return min;
    }

    template<class T>
    T SparseMatrix<T>::getRowSum(int row) {
        T sum = 0;
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            if (it->second->_row == row) {
                sum += it->second->val;
            }

        }
        return sum;
    }

    template<class T>
    T SparseMatrix<T>::getColSum(int col) {
        T sum = 0;
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            if (it->second->_col == col) {
                sum += it->second->val;
            }

        }
        return sum;
    }

    template<class T>
    T SparseMatrix<T>::getRowAvg(int row) {
        return this->getRowSum(row) / this->col;
    }

    template<class T>
    T SparseMatrix<T>::getColAvg(int col) {
        return this->getColSum(col) / this->row;
    }

    template <class T>
    T SparseMatrix<T>::getEigenvalue(int LoopNumber,  SparseMatrix<T> &result) {
        if (this->col != this->row) {
            throw ex::NotSquareException(this->row, this->col, "eigen value");
        }
        SparseMatrix<T> tempA (*this);//这是一个临时的矩阵，用来保存每一次被QR分解迭代的对象
        SparseMatrix<T> tempR(this->row,this->col);
        SparseMatrix<T>  tempQ(this->row,this->col);

        for (int i = 0; i < LoopNumber; i++) {
            tempA.QR(tempQ, tempR);
            SparseMatrix<T>tempRR(tempR);
            tempRR.crossProduct(tempQ);
            tempA = tempRR ;//下一次迭代的矩阵由RQ = Q'AQ给出

        }
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                if (i == j)
                    continue;
                tempA.setByIndex(i, j, 0);
            }
        }//将对角阵的非对角元写成0
        //totalQ = totalQ.transpose();
        result = tempA;
        return true;
    }

    template <class T>
    SparseMatrix<T> &SparseMatrix<T>::getEigenvector(SparseMatrix<T> &eigenvector,const T lamda) {
        SparseMatrix<T>eigenM(this->col,this->row);
        for (int i = 0; i <this->row; i++)
        {
            for (int j = 0; j <this->col; j++)
            {
                eigenM.setByIndex(i,j,this->getByIndex(i,j));
                if (i == j)
                    eigenM.setByIndex(i,j,this->getByIndex(i,j)-lamda);

            }
        }
        Gaussian_Eliminate(eigenvector, eigenM);
        return eigenvector;
    }

    template<class T>
    T SparseMatrix<T>::getTrace() {
        if (this->col != this->row) {
            throw ex::NotSquareException(*this, "matrix trace");
        }
        T trace = 0;
        for (int i = 0; i < this->col; i++) {
            trace += getByIndex(i, i);
        }
        return trace;
    }

    template<class T>
    T SparseMatrix<T>::getDeterminant() {
        if (this->col != this->row) {
            throw ex::NotSquareException(*this, "matrix determinant");
        }
        T det = 0;
        if (this->col == 1) {
            det = getByIndex(0, 0);
        } else if (this->col == 2) {
            det = this->getByIndex(0, 0) * getByIndex(1, 1) - getByIndex(0, 1) * getByIndex(1, 0); //一阶二阶直接计算
        } else {
            for (int k = 0; k < this->col; k++) {
                SparseMatrix<T> M(this->row - 1, this->col - 1); //为代数余子式申请内存
                for (int i = 0; i < this->col - 1; i++) {
                    for (int j = 0; j < this->col - 1; j++) {
                        M.setByIndex(i, j, this->getByIndex(i + 1, j < k ? j : j + 1)); //为代数余子式赋值
                    }
                }
                if (this->getByIndex(0, k) != 0) //如果是零可以直接不继续算
                {
                    det += this->getByIndex(0, k) * (T)(M.getDeterminant() * (T)(((2 + k) % 2) == 1 ? -1 : 1)); //从一第行展开，采用递归算法计算行列式
                }
            }
        }
        return det;
    }

    template<class T>
    void SparseMatrix<T>::reshape(int row, int col) {
        long _size = row * col;
        if (this->getSize() == _size) {
            this->setRow(col);
            this->setCol(row);
            for (auto it = tri_map.begin(); it != tri_map.end(); it++) {
                Triple<T> *tri = it->second;
                int temp = tri->_col;
                tri->_col = tri->_row;
                tri->_row = temp;
            }
        } else
            throw ex::MismatchedSizeException(this->getRow(), this->getCol(), row, col, "matrix reshaping");

    }

    template<class T>
    void SparseMatrix<T>::sliceRow(int row1, int row2) {
        for (auto it = tri_map.begin(); it != tri_map.end(); it++) {
            Triple<T> *tri = it->second;
            if (tri->_row < row1 || tri->_row > row2) {
                delete tri;
                tri_map.erase(it++);
            }
        }
    }

    template<class T>
    void SparseMatrix<T>::sliceCol(int col1, int col2) {
        for (auto it = tri_map.begin(); it != tri_map.end(); it++) {
            Triple<T> *tri = it->second;
            if (tri->_col < col1 || tri->_col > col2) {
                delete tri;
                tri_map.erase(it++);
            }
        }
    }

    template<class T>
    void SparseMatrix<T>::slice(int row1, int row2, int col1, int col2) {
        for (auto it = tri_map.begin(); it != tri_map.end(); it++) {
            Triple<T> *tri = it->second;
            if (tri->_row < row1 || tri->_row > row2 || tri->_col < col1 || tri->_col > col2) {
                delete tri;
                tri_map.erase(it++);
            }
        }
    }

    template<class T>
    SparseMatrix<T> *SparseMatrix<T>::convolve(SparseMatrix<T> &right, int stride, int padding) {
        if (right.row != right.col) {
            throw ex::NotSquareException(right, "doing matrix convolution");
        }
        int r = (this->row - right.row + 2 * padding) / stride + 1;
        int c = (this->col - right.col + 2 * padding) / stride + 1;
        SparseMatrix<T> *mat = new SparseMatrix<T>(r, c);
        SparseMatrix<T> rev(right);
        rev.reverse();

        SparseMatrix<T> ext(this->row + 2 * padding, this->col + 2 * padding);
        for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++) {
            Triple<T> *tri = it->second;
            ext.setByIndex(tri->_row + padding, tri->_col + padding, this->getByIndex(tri->_row, tri->_col));
        }

        for (auto i = rev.tri_map.begin(); i != rev.tri_map.end(); i++) {
            for (auto j = ext.tri_map.begin(); j != ext.tri_map.end(); j++) {
                Triple<T> *tri = i->second;
                Triple<T> *trj = j->second;
                int _r = trj->_row - tri->_row;
                int _c = trj->_col - tri->_col;
                if (0 <= _r && _r < r && 0 <= _c && _c < c) {
                    mat->setByIndex(_r, _c, mat->getByIndex(_r, _c) + tri->val * trj->val);
                }
            }
        }
        return mat;
    }


    template<class T>
    void SparseMatrix<T>::exponent(int exp) {
        if (this->getRow() != this->getCol())
            throw ex::NotSquareException(*this, "doing matrix exponential");
        SparseMatrix<T> temp(*this);
        for (int i = 0; i < exp; i++)
            this->crossProduct(temp);
    }

    template <class T>
    void SparseMatrix<T>::QR(SparseMatrix<T>&Q,SparseMatrix<T>&R){
        int i,j,k,r,m;
        T temp,sum,dr,cr,hr;
        SparseMatrix<T>ur(this->row*this->row,1);
        SparseMatrix<T>pr(this->row*this->row,1);
        SparseMatrix<T>wr(this->row*this->row,1);

        SparseMatrix<T>q1(this->row,this->row);
        SparseMatrix<T>emp(this->row,this->row);

        for(i=0;i<this->row;i++)//将a放入temp中

            for(j=0;j<this->col;j++)
            {
                emp.setByIndex(i,j,this->getByIndex(i,j));
            };
        for(i=0;i<this->row;i++)//定义单位矩阵

            for(j=0;j<this->col;j++)
            {
                if(i==j)Q.setByIndex(i,j,1);

                else Q.setByIndex(i,j,0);
            };

        for(r=0;r<this->col;r++)
        {
            temp=0;

            for(k=r;k<this->col;k++)

                temp+=fabs(this->getByIndex(k,r));

            if(temp>=0.0)
            {
                sum=0;

                for(k=r;k<this->row;k++)
                    sum+=this->getByIndex(k,r)*this->getByIndex(k,r);
                dr=sqrt(sum);
                if(this->getByIndex(r,r)>0.0)m=-1;

                else m=1;
                cr=m*dr;
                hr=cr*(cr-this->getByIndex(r,r));

                for(i=0;i<this->col;i++)//定义ur
                {
                    if(i<r)ur.setByIndex(i,0,0);

                    if(i==r)ur.setByIndex(i,0,this->getByIndex(r,r)-cr);

                    if(i>r)ur.setByIndex(i,0,this->getByIndex(i,r));
                };

                for(i=0;i<this->row;i++)//定义wr
                {
                    sum=0;

                    for(j=0;j<this->row;j++)

                        sum+=Q.getByIndex(i,j)*ur.getByIndex(j,0);
                    wr.setByIndex(i,0,sum);
                };

                for(i=0;i<this->row;i++)//定义qr
                    for(j=0;j<this->row;j++)
                    {
                        q1.setByIndex(i,j,Q.getByIndex(i,j)-wr.getByIndex(i,0)*ur.getByIndex(j,0)/hr);

                    };
                for(i=0;i<this->row;i++)//定义qr+1
                    for(j=0;j<this->row;j++)
                    {
                        Q.setByIndex(i,j,q1.getByIndex(i,j));
                    };
                for(i=0;i<this->col;i++)//定义pr
                {
                    sum=0;
                    for(j=0;j<this->col;j++)
                        sum+=this->getByIndex(j,i)*ur.getByIndex(j,0);
                    pr.setByIndex(i,0,sum/hr);
                };

                for(i=0;i<this->row;i++)
                    for(j=0;j<this->col;j++)
                    {
                        this->setByIndex(i,j,this->getByIndex(i,j)-ur.getByIndex(i,0)*pr.getByIndex(j,0));

                    };
            };
        };
        for(i=0;i<this->row;i++)
            for(j=0;j<this->col;j++)
            {
                if(fabs(this->getByIndex(i,j))<0.0)this->setByIndex(i,j,0);
            };

        for(i=0;i<this->row;i++)
            for(j=0;j<this->col;j++)
            {
                R.setByIndex(i,j,this->getByIndex(i,j));

            };

        for(i=0;i<this->row;i++)//将a取出
            for(j=0;j<this->col;j++)
            {
                this->setByIndex(i,j,emp.getByIndex(i,j));
            }
    }
    template <class T>
    void SparseMatrix<T>::Gaussian_Eliminate(SparseMatrix<T>&ans,SparseMatrix<T>&eigenmatirx){
        T max;
        short row;//每列的最大值及其行数
        T temp[eigenmatirx.col];
        for (int j = 0; j < eigenmatirx.col - 1; j++)//j为参考列
        {
            //找出列最大值及其行数
            max = abs(eigenmatirx.getByIndex(j,j));
            row = j;
            for (int i = j+1; i < eigenmatirx.col; i++)
            {
                if (abs(eigenmatirx.getByIndex(i,j)) > max)
                {
                    max = abs(eigenmatirx.getByIndex(i,j));
                    row = i;
                }
            }
            //将最大值行与第一行交换
            if (row != j)
            {
                for (int i = j; i < eigenmatirx.col; i++)
                    temp[i] = eigenmatirx.getByIndex(row,i);
                for (int i = j; i < eigenmatirx.col; i++)
                    eigenmatirx.setByIndex(row,i,eigenmatirx.getByIndex(j,i));
                for (int i = j; i < eigenmatirx.col; i++)
                    eigenmatirx.setByIndex(j,i,temp[i]);
            }
            //开始按列消元,即每次将除最大值行以外的一列清零
            for (int i = j + 1; i < eigenmatirx.col; i++)
                for (int k = j + 1; k <eigenmatirx.col; k++)
                    eigenmatirx.setByIndex(i,k,eigenmatirx.getByIndex(i,k)-eigenmatirx.getByIndex(i,j)/eigenmatirx.getByIndex(j,j)*eigenmatirx.getByIndex(j,k));

        }
        ans.setByIndex(eigenmatirx.col-1,0,1);
        for (int i = eigenmatirx.col - 2; i >= 0; i--)
        {
            ans.setByIndex(i,0,0);
            for (int j = eigenmatirx.row - 1; j > i; j--)
                ans.setByIndex(i,0,ans.getByIndex(i,0)-eigenmatirx.getByIndex(i,j)*ans.getByIndex(0,j));
            ans.setByIndex(i,0,ans.getByIndex(i,0)/eigenmatirx.getByIndex(i,i));

        }
    }

    template<class T>
    void SparseMatrix<T>::show() {
        if(this->getSize() < 2e5)
        {cout << "Sparse Matrix:" << endl;
        T mat[this->getRow()][this->getCol()];
        memset(mat, 0, sizeof(mat));
        for (auto i = tri_map.begin(); i != tri_map.end(); i++) {
            auto tri = i->second;
            mat[tri->_row][tri->_col] = tri->val;
        }
        for (size_t i = 0; i < this->getRow(); i++) {
            for (size_t j = 0; j < this->getCol(); j++) {
                cout << mat[i][j] << " ";
            }
            cout << endl;
        }}else{
            for (auto it = this->tri_map.begin(); it != this->tri_map.end(); it++)
            {
                Triple<T>* tri = it->second;
                cout<<"row="<<tri->_row<<" col="<<tri->_col<<" val="<<tri->val<<endl;
            }
        }
    }
}

#endif