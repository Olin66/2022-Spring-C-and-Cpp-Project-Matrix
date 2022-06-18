#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "basic-matrix.h"
#include "matrix-ex.h"
#include "matrix.h"
#include <cstring>
#include <iostream>
#include <math.h>
#include <set>
#include <string>
#include <vector>

namespace mat {
namespace ex {
class MatrixException;
class MismatchedSizeException;
class DuplicatedTripleException;
class NotSquareException;
class NoInverseException;
class InvalidSizeException;
} // namespace ex
template <typename T>
struct Triple {
    long _row;
    long _col;
    T val;

    bool operator==(const Triple<T> &right) {
        if (this->_row == right._row && this->_col == right._col) {
            return true;
        } else {
            return false;
        }
    }
};

template <typename>
class Matrix;

template <typename>
class BasicMatrix;

template <class T>
class SparseMatrix : public Matrix<T> {
  private:
    long size;
    std::set<Triple<T>> triples;

    inline void addEle(Triple<T> tri) {
        triples.insert(tri);
    }

  public:
    SparseMatrix(int, int);

    // SparseMatrix(const cv::Mat &mat);

    SparseMatrix(std::vector<std::vector<T>>);

    SparseMatrix(int, int, T *);

    SparseMatrix(int, int, std::vector<Triple<T>>);

    SparseMatrix(int, int, std::set<Triple<T>>);

    SparseMatrix(const SparseMatrix<T> &);

    SparseMatrix<T> &operator=(const SparseMatrix<T> &);

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

    T getEigenvalue(int LoopNumber,  BasicMatrix<T> &result);

    void Gaussian_Eliminate(BasicMatrix<T>&ans,BasicMatrix<T>&eigenmatirx);

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

    void QR(BasicMatrix<T>&Q,BasicMatrix<T>&R);

    std::set<Triple<T>> getTriples() const {
        return this->triples;
    }
};

template <class T>
SparseMatrix<T>::SparseMatrix(int row, int col) : Matrix<T>(row, col) {}

// template<class T>
// SparseMatrix<T>::SparseMatrix(const cv::Mat &mat) {

// }

template <class T>
SparseMatrix<T>::SparseMatrix(std::vector<std::vector<T>> mat) : Matrix<T>(mat.size(), mat[0].size()) {
    for (size_t i = 0; i < mat.size(); i++) {
        for (size_t j = 0; j < mat[i].size(); j++) {
            if (mat[i][j] != 0) {
                Triple<T> triple(i, j, mat[i][j]);
                triples.insert(triple);
            }
        }
    }
}

template <class T>
SparseMatrix<T>::SparseMatrix(int row, int col, T *_data) : Matrix<T>(row, col) {
    for (size_t i = 0; i < this->getSize(); i++) {
        if (_data[i] != 0) {
            Triple<T> triple(i / this->getCol(), i % this->getCol(), _data[i]);
            triples.insert(triple);
        }
    }
}

template <class T>
SparseMatrix<T>::SparseMatrix(int row, int col, std::vector<Triple<T>> mat) : Matrix<T>(row, col) {
    std::set<Triple<T>> temp(mat.begin(), mat.end());
    if (mat.size() != temp.size()) {
        throw ex::DuplicatedTripleException();
    }
    mat.assign(temp.begin(), temp.end());
    this->triples(mat);
}

template <class T>
SparseMatrix<T>::SparseMatrix(int row, int col, std::set<Triple<T>> mat) : Matrix<T>(row, col) {
    this->triples(mat);
}

template <class T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix<T> &right) : Matrix<T>(right.getRow(), right.getCol()) {
    this->triples(right.triples);
}

template <class T>
SparseMatrix<T> &SparseMatrix<T>::operator=(const SparseMatrix<T> &right) {
    this->setSize(right.getSize());
    this->setRow(right.getRow());
    this->setCol(right.getCol());
    this->triples(right.triples);
}

template <class T>
T SparseMatrix<T>::getByIndex(int _row, int _col) const {
    T target = 0;
    for (int i = 0; i < this->triples.size(); i++) {
        if (this->triples[i]._col == _col && this->triples[i]._row == _row) {
            return target = this->triples[i].val;
        }
    }
    return target;
}

template <class T>
void SparseMatrix<T>::setByIndex(int _row, int _col, T val) {
    for (int i = 0; i < this->triples.size(); i++) {
        if (this->triples[i]._row == _row && this->triples[i]._col == _col) {
            this->triples[i].val = val;
            return;
        }
    }
    Triple<T> triple(_row, _col, val);
    this->addEle(triple);
    return;
}

template <class T>
void SparseMatrix<T>::add(const BasicMatrix<T> &) {
}

template <class T>
void SparseMatrix<T>::add(const SparseMatrix<T> &right) {
    SparseMatrix<T> mat(this->row, this->col);
    for (int i = 0; i < this->triples.size(); i++) {
        for (int j = 0; j < right.triples.size(); j++) {
            Triple<T> LP = triples[i];
            Triple<T> RP = triples[j];
            if (LP._row == RP._row && LP._row == RP._row) {
                mat.setByIndex(LP._row, LP._col, LP.val + RP.val);
            }
        }
    }
    *this = mat;
}

template <class T>
void SparseMatrix<T>::subtract(const BasicMatrix<T> &right) {
       
}

template <class T>
void SparseMatrix<T>::subtract(const SparseMatrix<T> &right) {
    SparseMatrix<T> mat(this->row, this->col);
    for (int i = 0; i < this->triples.size(); i++) {
        for (int j = 0; j < right.triples.size(); j++) {
            Triple<T> LP = triples[i];
            Triple<T> RP = triples[j];
            if (LP._row == RP._row && LP._row == RP._row) {
                mat.setByIndex(LP._row, LP._col, LP.val - RP.val);
            }
        }
    }
    *this = mat;
}

template <class T>
void SparseMatrix<T>::scalarMultiply(T x) {
    for (int i = 0; i < this->triples.size(); i++) {
        this->triples[i] *= x;
    }
}

template <class T>
void SparseMatrix<T>::scalarDivide(T x) {
    for (int i = 0; i < this->triples.size(); i++) {
        this->triples[i] /= x;
    }
}

template <class T>
void SparseMatrix<T>::dotProduct(const BasicMatrix<T> &) {
}

template <class T>
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
        for (int i = 0; i < this->triples.size(); i++) {
            for (int j = 0; j < right.triples.size(); j++) {
                mat::Triple<T> &l_point = this->triples[i];
                mat::Triple<T> &r_point = right.triples[j];
                if (l_point._row == r_point._row) {
                    mat.setByIndex(r_point._row, r_point._col, l_point.val * r_point.val);
                }
            }
        }
        *this = mat;
    } else {
        SparseMatrix<T> mat(r, c);
        for (int i = 0; i < this->triples.size(); i++) {
            for (int j = 0; j < right.triples.size(); j++) {
                mat::Triple<T> &l_point = this->triples[i];
                mat::Triple<T> &r_point = right.triples[i];
                if (l_point._row == r_point._row &&
                    l_point._col == r_point._col) {
                    mat.setByIndex(r_point._row, r_point._col, l_point.val * r_point.val);
                }
            }
        }
        *this = mat;
    }
}

template <class T>
void SparseMatrix<T>::crossProduct(const BasicMatrix<T> &) {
}

template <class T>
void SparseMatrix<T>::crossProduct(const SparseMatrix<T> &right) {
    if (this->col != right.col) {
        throw ex::MismatchedSizeException(*this, right, "matrix cross product");
    }
    int r = this->row;
    int c = right.col;
    SparseMatrix<T> mat(r, c);
    for (int i = 0; i < this->triples.size(); i++) {
        Triple<T> &l_point = this->triples[i];
        for (int j = 0; j < right.triples.size(); j++) {
            Triple<T> &r_point = right.triples[j];
            if (l_point._col == r_point._row) {
                mat.setByIndex(l_point._row, r_point._col, mat.getByIndex(l_point._row, r_point._col) + l_point.val * r_point.val);
            }
        }
    }
    *this = mat;
}

template <class T>
void SparseMatrix<T>::transpose() {
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

template <class T>
void SparseMatrix<T>::reverse() {
}

template <class T>
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

template <class T>
T SparseMatrix<T>::getMax() {
    T max = 0;
    for (int i = 0; i < this->triples.size(); i++) {
        if (max < this->triples[i].val) {
            max = this->triples[i].val;
        }
    }
    return max;
}

template <class T>
T SparseMatrix<T>::getMin() {
    T min = 0;
    for (int i = 0; i < this->triples.size(); i++) {
        if (min > this->triples[i].val) {
            min = this->triples[i].val;
        }
    }
    return min;
}

template <class T>
T SparseMatrix<T>::getSum() {
    T sum = 0;
    for (int i = 0; i < this->triples.size(); i++) {
        sum += this->triples[i].val;
    }
    return sum;
}

template <class T>
T SparseMatrix<T>::getAvg() {
    return this->getSum / (this->col * this->row);
}

template <class T>
T SparseMatrix<T>::getRowMax(int row) {
    T max = 0;
    for (int i = 0; i < this->triples.size(); i++) {
        if (max < this->triples[i].val && this->triples[i]._row == row) {
            max = this->triples[i].val;
        }
    }
    return max;
}

template <class T>
T SparseMatrix<T>::getColMax(int col) {
    T max = 0;
    for (int i = 0; i < this->triples.size(); i++) {
        if (max < this->triples[i].val && this->triples[i]._col == col) {
            max = this->triples[i].val;
        }
    }
    return max;
}

template <class T>
T SparseMatrix<T>::getRowMin(int row) {
    T min = 0;
    for (int i = 0; i < this->triples.size(); i++) {
        if (min > this->triples[i].val && this->triples[i]._row == row) {
            min = this->triples[i].val;
        }
    }
    return min;
}

template <class T>
T SparseMatrix<T>::getColMin(int col) {
    T min = 0;
    for (int i = 0; i < this->triples.size(); i++) {
        if (min > this->triples[i].val && this->triples[i]._col == col) {
            min = this->triples[i].val;
        }
    }
    return min;
}

template <class T>
T SparseMatrix<T>::getRowSum(int row) {
    T sum = 0;
    for (int i = 0; i < this->triples.size(); i++) {
        if (this->triples[i]._row == row) {
            sum += this->triples[i].val;
        }
    }
    return sum;
}

template <class T>
T SparseMatrix<T>::getColSum(int col) {
    T sum = 0;
    for (int i = 0; i < this->triples.size(); i++) {
        if (this->triples[i]._col == col) {
            sum += this->triples[i].val;
        }
    }
    return sum;
}

template <class T>
T SparseMatrix<T>::getRowAvg(int row) {
    return this->getRowSum(row) / this->col;
}

template <class T>
T SparseMatrix<T>::getColAvg(int col) {
    return this->getColSum(col) / this->row;
}

template <class T>
T SparseMatrix<T>::getEigenvalue(int LoopNumber,  BasicMatrix<T> &result) {
    if (this->col != this->row) {
        throw ex::NotSquareException(this->row, this->col, "eigen value");
    }
    BasicMatrix<T> tempA (*this);//这是一个临时的矩阵，用来保存每一次被QR分解迭代的对象
	BasicMatrix<T> tempR(this->row,this->col);
    BasicMatrix<T>  tempQ(this->row,this->col);
    
	for (int i = 0; i < LoopNumber; i++) {
		tempA.QR(tempQ, tempR);
        BasicMatrix<T>tempRR(tempR);
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
BasicMatrix<T> &BasicMatrix<T>::getEigenvector(BasicMatrix<T> &eigenvector,const T lamda) {
    BasicMatrix<T>eigenM(this->col,this->row);
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

template <class T>
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

template <class T>
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
                det += this->getByIndex(0, k) * M.getDeterminant() * (((2 + k) % 2) == 1 ? -1 : 1); //从一第行展开，采用递归算法计算行列式
            }
        }
    }
    return det;
}

template <class T>
void SparseMatrix<T>::reshape(int row, int col) {
    long _size = row * col;
    if (this->getSize() == _size) {
        this->setRow(col);
        this->setCol(row);


    } else
        throw ex::MismatchedSizeException(this->getRow(), this->getCol(), row, col, "matrix reshaping");
}

template <class T>
void SparseMatrix<T>::sliceRow(int row1, int row2) {
}

template <class T>
void SparseMatrix<T>::sliceCol(int col1, int col2) {
}

template <class T>
void SparseMatrix<T>::slice(int row1, int row2, int col1, int col2) {
}

template <class T>
Matrix<T> &SparseMatrix<T>::convolve(BasicMatrix<T> &, int stride, int padding) {
}

template <class T>
Matrix<T> &SparseMatrix<T>::convolve(SparseMatrix<T> &, int stride, int padding) {
}

template <class T>
void SparseMatrix<T>::exponent(int exp) {
}

template <class T>
void BasicMatrix<T>::QR(BasicMatrix<T>&Q,BasicMatrix<T>&R){
     int i,j,k,r,m;
  T temp,sum,dr,cr,hr;
  BasicMatrix<T>ur(this->row*this->row,1);
  BasicMatrix<T>pr(this->row*this->row,1);
  BasicMatrix<T>wr(this->row*this->row,1);

  BasicMatrix<T>q1(this->row,this->row);
  BasicMatrix<T>emp(this->row,this->row);
 
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
void BasicMatrix<T>::Gaussian_Eliminate(BasicMatrix<T>&ans,BasicMatrix<T>&eigenmatirx){
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


} // namespace mat

#endif