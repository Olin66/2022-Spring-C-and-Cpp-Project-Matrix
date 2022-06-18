#ifndef BASIC_MATRIX_H
#define BASIC_MATRIX_H
#include "matrix-ex.h"
#include "matrix.h"
#include "sparse-matrix.h"
// #include <opencv2/opencv.hpp>
#include <cstring>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include<complex>
using namespace std;
// using namespace cv;
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

    template <typename>
    class Matrix;

    template <typename>
    class SparseMatrix;

    template <typename>
    struct Triple;

    template <class T>
    class BasicMatrix : public Matrix<T> {
    private:
        T *m_data;

    public:
        BasicMatrix(int, int, T val=0);

        BasicMatrix(int, int, T *);

        // BasicMatrix(const cv::Mat &mat);

        BasicMatrix(std::vector<std::vector<T>>);

        BasicMatrix(const BasicMatrix<T> &);

        BasicMatrix<T> &operator=(const BasicMatrix<T> &);

        ~BasicMatrix<T>();

        void add(const BasicMatrix<T> &);

        void subtract(const BasicMatrix<T> &);

        void scalarMultiply(T);

        void scalarDivide(T);

        void dotProduct(const BasicMatrix<T> &);

        void crossProduct(const BasicMatrix<T> &);

        void transpose();

        void inverse();

        void reverse();

        void conjugate();

        void Hessenberg();

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

        bool getEigenvalue(int LoopNumber, double error, BasicMatrix<T> &result);

        BasicMatrix<T> &getEigenvector(BasicMatrix<T>&eigenvector);

        bool getEigen(BasicMatrix<T> &eigenvector, BasicMatrix<T> &eigenvalue, double error, double iterator);//雅克比法计算特征值和特征向量

        T getTrace();

        T getDeterminant();

        void reshape(int row, int col);

        void loop(T[]);

        void sliceRow(int row1, int row2);

        void sliceCol(int col1, int col2);

        void slice(int row1, int row2, int col1, int col2);

        BasicMatrix<T> *convolve(BasicMatrix<T> &, int stride = 1, int padding = 0);

        void exponent(int exp);

        void show();

        BasicMatrix<T> operator+(const BasicMatrix<T> &);

        BasicMatrix<T> operator-(const BasicMatrix<T> &);

        template <typename P>
        BasicMatrix<T> operator*(P);

        template <typename P>
        friend BasicMatrix<T> operator*(P val, BasicMatrix<T>& right){
            return right * val;
        }

        T *getData() {
            return this->m_data;
        }
    };

    template <class T>
    BasicMatrix<T>::BasicMatrix(int row, int col, T val) : Matrix<T>(row, col) {
        this->m_data = new T[this->getSize()];
        if (val == 0) std::memset(m_data, 0, sizeof(T) * this->getSize());
        else {
            for (size_t i = 0; i < this->getSize(); i++)
                m_data[i] = val;
        }
    }

    template <class T>
    BasicMatrix<T>::BasicMatrix(int row, int col, T *_data) : Matrix<T>(row, col) {
        this->m_data = new T[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++)
            this->m_data[i] = _data[i];
    }

    // template<class T>
    // BasicMatrix<T>::BasicMatrix(const cv::Mat &mat): Matrix<T>(mat) {
    //     if (mat.type() != CV_8UC1)
    //         throw ex::InvalidChannelDepth(mat.type);
    //     for (size_t i = 0; i < this->getSize(); i++)
    //         this->m_data[i] = mat.data[i];

    // }

    template <class T>
    BasicMatrix<T>::BasicMatrix(std::vector<std::vector<T>> mat) : Matrix<T>(mat.size(), mat[0].size()) {
        this->m_data = new T[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++)
            this->m_data[i] = mat[i / this->getCol()][i % this->getCol()];
    }

    template <class T>
    BasicMatrix<T>::BasicMatrix(const BasicMatrix<T> &right) : Matrix<T>(right.getRow(), right.getCol()) {
        this->m_data = new T[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++)
            this->m_data[i] = right.m_data[i];
    }

    template <class T>
    BasicMatrix<T> &BasicMatrix<T>::operator=(const BasicMatrix<T> &right) {
        if (this == &right) return (*this);
        this->setRow(right.getRow());
        this->setCol(right.getCol());
        this->setSize(right.getSize());
        delete[] m_data;
        m_data = new T[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++)
            this->m_data[i] = right.m_data[i];
        return (*this);
    }

    template<class T>
    BasicMatrix<T>::~BasicMatrix<T>(){
        delete[] m_data;
    }

    template <class T>
    T BasicMatrix<T>::getByIndex(int _row, int _col) const {
        return m_data[_row * this->col + _col];
    }

    template <class T>
    void BasicMatrix<T>::setByIndex(int _row, int _col, T val) {
        m_data[_row * this->col + _col] = val;
    }

    template <class T>
    void BasicMatrix<T>::add(const BasicMatrix<T> &right) {
        if (this->getRow() != right.getRow() || this->getCol() != right.getCol())
            throw ex::MismatchedSizeException(*this, right,
                                              "matrix addition");
        for (size_t i = 0; i < this->getSize(); i++) {
            this->m_data[i] = this->m_data[i] + right.m_data[i];
        }
    }

    template <class T>
    void BasicMatrix<T>::subtract(const BasicMatrix<T> &right) {
        if (this->getRow() != right.getRow() || this->getCol() != right.getCol()) {
            throw ex::MismatchedSizeException(*this, right,
                                              "matrix subtraction");
        }
        for (size_t i = 0; i < this->getSize(); i++) {
            this->m_data[i] = this->m_data[i] - right.m_data[i];
        }
    }

    template <class T>
    void BasicMatrix<T>::scalarMultiply(T val) {
        for (size_t i = 0; i < this->getSize(); i++)
        {
            this->m_data[i] = this->m_data[i] * val;
        }
    }

    template <class T>
    void BasicMatrix<T>::scalarDivide(T val) {
        for (size_t i = 0; i < this->getSize(); i++)
        {
            this->m_data[i] = this->m_data[i] / val;
        }
    }

    template <class T>
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
                    mat.setByIndex(i, j, this->getByIndex(i, 0) * (right.getByIndex(i, j)));
                }
            }
            *this = mat;
        } else {
            this->col = c;
            for (int i = 0; i < r; i++) {
                for (int j = 0; j < c; j++) {
                    this->setByIndex(i, j, this->getByIndex(i, j) * (right.getByIndex(i, j)));
                }
            }
        }
    }

    template <class T>
    void BasicMatrix<T>::crossProduct(const BasicMatrix<T> &right) {
        if (this->col != right.row) {
            throw ex::MismatchedSizeException(*this, right, "matrix cross product");
        }
        int r = this->row;
        int c = right.col;
        BasicMatrix<T> mat(r, c);
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

    template <class T>
    void BasicMatrix<T>::transpose() {
        BasicMatrix<T> mat(this->col, this->row);
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                mat.setByIndex(j, i, this->getByIndex(i, j));
            }
        }
        *this = mat;
    }

    template <class T>
    void BasicMatrix<T>::inverse() { //用伴随矩阵求逆
        T det = this->getDeterminant();
        if (this->col != this->col || det == 0) //行列式为0的矩阵不可逆
        {
            throw ex::NoInverseException(*this, "matrix inverse");
        }
        BasicMatrix<T> temp(this->row - 1, this->col - 1);
        BasicMatrix<T> adjoint(this->row, this->col);
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
    void BasicMatrix<T>::reverse() {
        T *r_data = new T[this->getSize()];
        for (int i = 0; i < this->getSize(); i++)
            r_data[i] = m_data[this->getSize() - i - 1];
        delete [] m_data;
        m_data = r_data;
    }

    template <class T>
    void BasicMatrix<T>::conjugate() {
        for (int i = 0; i < this->getRow(); i++)
        {
            for (int j = 0; j < this->getCol(); j++)
            {
                if(imag(this->getByIndex(i,j))!=0){
                    T temp=real(this->getByIndex(i,j))-imag(this->getByIndex(i,j));
                    this->setByIndex(i,j,temp);
                }
            }

        }

    }


    template <class T>
    void BasicMatrix<T>::Hessenberg() { //求海森堡矩阵，上三角化
        BasicMatrix<T> A(this->row, this->col);
        int k;
        int i;
        int j;
        int Max;
        double temp;
        for (i = 0; i < this->col; i++) {

            for (j = 0; j < this->col; j++) {
                A.setByIndex(i, j, this->getByIndex(i, j));
            }
        }
        for (k = 1; k < this->col - 1; k++) {
            i = k - 1;
            Max = k;
            temp = abs(A.getByIndex(k, i));
            for (j = k + 1; j < this->col; j++) {
                if (abs(A.getByIndex(j, i)) > temp) {
                    temp = abs(A.getByIndex(j, i));
                    Max = j;
                }
            }

            this->setByIndex(0, 0, A.getByIndex(Max, i));
            i = Max;
            if (this->getByIndex(0, 0) != 0) {
                if (i != k) {
                    for (j = k - 1; j < this->col; j++) {
                        temp = A.getByIndex(i, j);
                        A.setByIndex(i, j, A.getByIndex(k, j));
                        A.setByIndex(k, j, temp);
                    }
                    for (j = 0; j < this->col; j++) {
                        temp = A.getByIndex(j, i);
                        A.setByIndex(j, i, A.getByIndex(j, k));
                        A.setByIndex(j, k, temp);
                    }
                }
                for (i = k + 1; i < this->col; i++) {
                    temp = A.getByIndex(i, k - 1) / this->getByIndex(0, 0);
                    A.setByIndex(i, k - 1, 0);
                    for (j = k; j < this->col; j++) {
                        A.setByIndex(i, j, -temp * A.getByIndex(k, j));
                    }
                    for (j = 0; j < this->col; j++) {
                        A.setByIndex(j, k, A.getByIndex(j, k) + temp * A.getByIndex(j, i));
                    }
                }
            }
        }
        for (i = 0; i < this->col; i++) {

            for (j = 0; j < this->col; j++) {
                this->setByIndex(i, j, A.getByIndex(i, j));
            }
        }
    }
    template <class T>
    T BasicMatrix<T>::getMax() {
        T max = getByIndex(0, 0);
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                T temp = getByIndex(i, j);
                if (max < temp) {
                    max = temp;
                }
            }
        }
        return max;
    }

    template <class T>
    T BasicMatrix<T>::getMin() {
        T min = getByIndex(0, 0);
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                T temp = getByIndex(i, j);
                if (min > temp) {
                    min = temp;
                }
            }
        }
        return min;
    }

    template <class T>
    T BasicMatrix<T>::getSum() {
        T sum = 0;
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                sum += getByIndex(i, j);
            }
        }
        return sum;
    }

    template <class T>
    T BasicMatrix<T>::getAvg() {
        int sum = getSum();
        return sum / this->getSize();
    }

    template <class T>
    bool BasicMatrix<T>::getEigenvalue(int LoopNumber, double error, BasicMatrix<T> &result) {
        if (this->col != this->row) {
            throw ex::NotSquareException(this->row, this->col, "getting the eigenvalue");
        }
        int i;
        int j;
        int k;
        int t;
        int m;
        int loop;
        T b;
        T c;
        T d;
        T g;
        T xy;
        T p;
        T q;
        T r;
        T x;
        T s;
        T e;
        T f;
        T z;
        T y;
        T temp;
        this->Hessenberg();
        BasicMatrix<T> A(*this);
        BasicMatrix<T> B(this->row, 2);
        result = B;
        m = this->getCol();
        loop = LoopNumber;
        while (m) {
            t = m - 1;
            while (t > 0) {
                temp = abs(A.getByIndex(t - 1, t - 1));
                temp += temp * error;
                if (abs(A.getByIndex(t, t - 1)) > temp) {
                    t--;
                } else {
                    break;
                }
            }
            if (t == m - 1) {
                result.setByIndex(m - 1, 0, A.getByIndex(m - 1, m - 1));
                result.setByIndex(m - 1, 1, 0);
                m -= 1;
                loop = LoopNumber;
            } else if (t == m - 2) {
                b = -A.getByIndex(m - 1, m - 1) - A.getByIndex(m - 2, m - 2);
                c = A.getByIndex(m - 1, m - 1) * A.getByIndex(m - 2, m - 2) - A.getByIndex(m - 1, m - 2) * A.getByIndex(m - 2, m - 1);
                d = b * b - 4 * c;
                y = sqrt(abs(d));
                if (d > 0) {
                    xy = 1;
                    if (b < 0) {
                        xy = -1;
                    }
                    result.setByIndex(m - 1, 0, -(b + xy * y) / 2);
                    result.setByIndex(m - 1, 1, 0);
                    result.setByIndex(m - 2, 0, c / result.getByIndex(m - 1, 0));
                    result.setByIndex(m - 2, 1, 0);

                } else {
                    result.setByIndex(m - 1, 0, -b / 2);
                    result.setByIndex(m - 2, 0, -b / 2);
                    result.setByIndex(m - 1, 1, y / 2);
                    result.setByIndex(m - 2, 1, -y / 2);
                }
                m -= 2;
                loop = LoopNumber;

            } else {
                if (loop < 1) {
                    cout << "no eigenvalue" << endl;
                    return false;
                }
                loop--;
                j = t + 2;
                while (j < m) {
                    A.setByIndex(j, j - 2, 0);
                    j++;
                }
                j = t + 3;
                while (j < m) {
                    A.setByIndex(j, j - 3, 0);
                    j++;
                }
                k = t;
                while (k < m - 1) {
                    if (k != t) {
                        p = A.getByIndex(k, k - 1);
                        q = A.getByIndex(k + 1, k - 1);
                        if (k != m - 2) {
                            r = A.getByIndex(k + 2, k - 1);
                        } else {
                            r = 0;
                        }

                    } else {
                        b = A.getByIndex(m - 1, m - 1);
                        c = A.getByIndex(m - 2, m - 2);
                        x = b + c;
                        y = b * c - A.getByIndex(m - 2, m - 1) * A.getByIndex(m - 1, m - 2);
                        p = A.getByIndex(t, t) * (A.getByIndex(t, t) - x) + A.getByIndex(t, t + 1) * A.getByIndex(t + 1, t) + y;
                        q = A.getByIndex(t + 1, t) * (A.getByIndex(t, t) + A.getByIndex(t + 1, t + 1) - x);
                        r = A.getByIndex(t + 1, t) * A.getByIndex(t + 2, t + 1);
                    }
                    if (p != 0 || q != 0 || r != 0) {
                        if (q < 0) {
                            xy = -1;
                        } else {
                            xy = 1;
                        }
                        s = xy * sqrt(p * p + q * q + r * r);
                        if (k != t) {
                            A.setByIndex(k, k - 1, -s);
                        }
                        e = -q / s;
                        f = -r / s;
                        x = -p / s;
                        y = -x - f * r / (p + s);
                        g = e * r / (p + s);
                        z = -x - e * q / (p + s);
                        for (j = k; j < m; j++) {
                            b = A.getByIndex(k, j);
                            c = A.getByIndex(k + 1, j);
                            p = x * b + e * c;
                            q = e * b + y * c;
                            r = f * b + g * c;
                            if (k != m - 2) {
                                b = A.getByIndex(k + 2, j);
                                p += f * b;
                                q += g * b;
                                r += z * b;
                                A.setByIndex(k + 2, j, r);
                            }
                            A.setByIndex(k + 1, j, q);
                            A.setByIndex(k, j, p);
                        }
                        j = k + 3;
                        if (j > m - 2) {
                            j = m - 1;
                        }
                        for (i = t; i < j + 1; i++) {
                            b = A.getByIndex(i, k);
                            c = A.getByIndex(i, k + 1);
                            p = x * b + e * c;
                            q = e * b + y * c;
                            r = f * b + g * c;
                            if (k != m - 2) {
                                b = A.getByIndex(i, k + 2);
                                p += f * b;
                                q += g * b;
                                z += z * b;
                                A.setByIndex(i, k + 2, r);
                            }
                            A.setByIndex(i, k + 1, q);
                            A.setByIndex(i, k, p);
                        }
                    }
                    k++;
                }
            }
        }
        return true;
    }

    template <class T>
    T BasicMatrix<T>::getRowMin(int row) {
        T min = this->getByIndex(row, 0);
        for (int i = 0; i < this->col; i++) {
            if (min > this->getByIndex(row, i)) {
                min = this->getByIndex(row, i);
            }
        }
        return min;
    }

    template <class T>
    T BasicMatrix<T>::getColMin(int col) {
        T min = this->getByIndex(0, col);
        for (int i = 0; i < this->row; i++) {
            if (min > this->getByIndex(i, col)) {
                min = this->getByIndex(i, col);
            }
        }
        return min;
    }

    template <class T>
    T BasicMatrix<T>::getRowMax(int row) {
        T max = this->getByIndex(this->row, 0);
        for (int i = 0; i < this->col; i++) {
            if (max < this->getByIndex(row, i)) {
                max = this->getByIndex(row, i);
            }
        }
        return max;
    }

    template <class T>
    T BasicMatrix<T>::getColMax(int col) {
        T max = this->getByIndex(0, col);
        for (int i = 0; i < this->row; i++) {
            if (max > this->getByIndex(i, col)) {
                max = this->getByIndex(i, col);
            }
        }
        return max;
    }

    template <class T>
    T BasicMatrix<T>::getRowSum(int row) {
        T sum = 0;
        for (int i = 0; i < this->col; i++) {
            sum += this->getByIndex(row, i);
        }
        return sum;
    }

    template <class T>
    T BasicMatrix<T>::getColSum(int col) {
        T sum = 0;
        for (int i = 0; i < this->row; i++) {
            sum += this->getByIndex(i, col);
        }
        return sum;
    }

    template <class T>
    T BasicMatrix<T>::getRowAvg(int row) {
        return this->getRowSum(row) / this->col;
    }

    template <class T>
    T BasicMatrix<T>::getColAvg(int col) {
        return this->getColSum(col) / this->row;
    }

    template <class T>
    BasicMatrix<T>& BasicMatrix<T>::getEigenvector(BasicMatrix<T> &eigenvector) {
        T iterator[this->col];
        for (int i = 0; i < this->col; i++)
        {
            iterator[i]=1;
        }


        T y=1.0;
        T z;

        do
        {
            z=y;
            this->loop(iterator);
            T max=iterator[0];
            for (int i = 0; i < this->col; i++)//求最大值
            {
                if(y<iterator[i])
                    max=iterator[i];
            }
            y=max;
            for (int i = 0; i < this->col; i++)
            {
                eigenvector.setByIndex(i,0,iterator[i]/y);
                iterator[i]=eigenvector.getByIndex(i,0);
            }




        } while (fabs(z-y)>=0.000001);
        cout<<"y= "<<y<<endl;
        return eigenvector;
    }


    template <class T>
    bool BasicMatrix<T>::getEigen(BasicMatrix<T> &eigenvector, BasicMatrix<T> &eigenvalue, double error, double iterator) {
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

        BasicMatrix<T> temp(this->col, this->col);
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
    }

    template <class T>
    T BasicMatrix<T>::getTrace() {
        if (this->col != this->row) {
            throw ex::NotSquareException(*this, "getting the trace");
        }
        T trace = 0;
        for (int i = 0; i < this->col; i++) {
            trace += getByIndex(i, i);
        }
        return trace;
    }

    template <class T>
    T BasicMatrix<T>::getDeterminant() {
        if (this->col != this->row) {
            throw ex::NotSquareException(*this, "getting the determinant");
        }
        T det = 0;
        if (this->col == 1) {
            det = getByIndex(0, 0);
        } else if (this->col == 2) {
            det = this->getByIndex(0, 0) * getByIndex(1, 1) - getByIndex(0, 1) * getByIndex(1, 0); //一阶二阶直接计算
        } else {
            for (int k = 0; k < this->col; k++) {
                BasicMatrix<T> M(this->row - 1, this->col - 1); //为代数余子式申请内存
                for (int i = 0; i < this->col - 1; i++) {
                    for (int j = 0; j < this->col - 1; j++) {
                        M.setByIndex(i, j, this->getByIndex(i + 1, j < k ? j : j + 1)); //为代数余子式赋值
                    }
                }
                if (this->getByIndex(0, k) != 0) //如果是零可以直接不继续算
                {
                    det += this->getByIndex(0, k) * M.getDeterminant() * (((2 + k) % 2) == 1 ? -1 : 1); //从一第行展开，采用递归算法计算行列式
                }
            }
        }
        return det;
    }

    template <class T>
    void BasicMatrix<T>::reshape(int row, int col) {
        long _size = row * col;
        if (this->getSize() == _size) {
            this->setRow(col);
            this->setCol(row);
        } else
            throw ex::MismatchedSizeException(this->getRow(), this->getCol(), row, col, "matrix reshaping");
    }

    template <class T>
    void BasicMatrix<T>::sliceRow(int row1, int row2) {
        T* _new_data_;
        if (row1 == row2) {
            this->setRow(1);
            this->setSize(this->getCol());
            _new_data_ = new T[this->getCol()];
            int start_index = row1 * this->getCol();
            for (size_t i = 0; i < this->getCol(); i++)
            {
                _new_data_[i] = this->getData()[start_index++];
            }
        }else if (row1 < row2){
            long temp = this->getCol() * (row2-row1+1);
            this->setRow(row2-row1+1);
            this->setSize(temp);
            _new_data_ = new T[temp];
            int start_index = row1 * this->getCol();
            for (size_t i = 0; i < temp; i++)
            {
                _new_data_[i] = this->getData()[start_index++];
            }
        } else{
            throw ex::InvalidSizeException("slicing the row", 2, row1, row2);
        }
        delete[] this->m_data;
        this->m_data = _new_data_;
    }

    template <class T>
    void BasicMatrix<T>::sliceCol(int col1, int col2) {
        T* _new_data_;
        if (col1 == col2){
            this->setCol(1);
            this->setSize(this->getRow());
            _new_data_ = new T[this->getRow()];
            for (size_t i = 0; i < this->getRow(); i++)
            {
                _new_data_[i] = col1 + i * this->getCol();
            }
        }else if (col1 < col2){
            long temp = this->getRow() * (col2-col1+1);
            _new_data_ = new T[temp];
            long k = 0;
            for (size_t i = 0; i < this->getRow(); i++)
            {
                for (size_t j = col1; j <= col2; j++)
                {
                    _new_data_[k++] = this->m_data[j + i *this-> getCol()];
                }
            }
            this->setCol(col2-col1+1);
            this->setSize(temp);
        } else{
            throw ex::InvalidSizeException("slicing the column", 3, col1, col2);
        }
        delete[] this->m_data;
        this->m_data = _new_data_;
    }

    template <class T>
    void BasicMatrix<T>::slice(int row1, int row2, int col1, int col2) {
        if (row1 > row2) throw ex::InvalidSizeException("slicing the row", 2, row1, row2);
        if (col1 > col2) throw ex::InvalidSizeException("slicing the column", 3, col1, col2);
        long temp = (row2-row1+1) * (col2-col1+1);
        T* _new_data_ = new T[temp];
        long k = 0;
        for (size_t i = row1;i <= row2; i++)
        {
            for (size_t j = col1; j <= col2; j++)
            {
                _new_data_[k++] = this->m_data[j + i * this->getCol()];
            }
        }
        this->setRow(row2-row1+1);
        this->setCol(col2-col1+1);
        this->setSize(temp);
    }

    template<class T>
    BasicMatrix<T>* BasicMatrix<T>::convolve(BasicMatrix<T> &right, int stride, int padding) {
        if (right.row != right.col) {
            throw ex::NotSquareException(right, "doing matrix convolution");
        }
        int r = (this->row - right.row + 2 * padding) / stride + 1;
        int c = (this->col - right.col + 2 * padding) / stride + 1;
        BasicMatrix<T> mat(r, c);
        BasicMatrix<T> rev(right); rev.reverse();

        BasicMatrix<T> ext(this->row + 2 * padding, this->col + 2 * padding);
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                ext.setByIndex(i+padding, j+padding, this->getByIndex(i, j));
            }
        }

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < rev.row; k++) {
                    for (int t = 0; t < rev.col; t++) {
                        mat.setByIndex(i, j, mat.getByIndex(i, j) + rev.getByIndex(k, t) * ext.getByIndex(k+i, t+j));
                    }
                }
            }
        }
        return new BasicMatrix<T>(mat);
    }

    template <class T>
    void BasicMatrix<T>::exponent(int exp) {
        if(this->getRow() != this->getCol())
            throw ex::NotSquareException(*this, "doing matrix exponential");
        BasicMatrix<T> temp(*this);
        for (int i = 0; i < exp; i++)
            this->crossProduct(temp);
    }

    template <class T>
    void BasicMatrix<T>::show() {
        using namespace std;
        cout << "Basic Matrix:" << endl;
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++)
                cout << getByIndex(i, j) << " ";
            cout << endl;
        }
    }

    template <class T>
    BasicMatrix<T> BasicMatrix<T>::operator+(const BasicMatrix<T> & right) {
        if (this->getRow() != right.getRow() || this->getCol() != right.getCol())
            throw ex::MismatchedSizeException(*this, right,
                                              "matrix addition");
        T new_data[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++) {
            new_data[i] = this->m_data[i] + right.m_data[i];
        }
        return BasicMatrix<T>(this->getRow(), this->getCol(), new_data);
    }

    template <class T>
    BasicMatrix<T> BasicMatrix<T>::operator-(const BasicMatrix<T> & right) {
        if (this->getRow() != right.getRow() || this->getCol() != right.getCol())
            throw ex::MismatchedSizeException(*this, right,
                                              "matrix addition");
        T new_data[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++) {
            new_data[i] = this->m_data[i] - right.m_data[i];
        }
        return BasicMatrix<T>(this->getRow(), this->getCol(), new_data);
    }

    template <class T>
    template <typename P>
    BasicMatrix<T> BasicMatrix<T>::operator*(P val) {
        T new_data[this->getSize()];
        for (size_t i = 0; i < this->getSize(); i++)
        {
            new_data[i] = this->m_data[i] * val;
        }
        return BasicMatrix<T>(this->getRow(), this->getCol(), new_data);
    }

    template <class T>
    void BasicMatrix<T>::loop(T iterator[]){
        T S,U[this->col];//loopstart
        for (int i = 0; i <this->col; i++)
        {
            U[i]=iterator[i];
        }
        for (int i = 0; i <this->col; i++)
        {
            S=0.0;
            for (int j = 0; j < this->row ;j++)
            {
                S+=this->getByIndex(i,j)*U[j];
            }
            iterator[i]=S;
        }//loop end;
    };
} // namespace mat

#endif