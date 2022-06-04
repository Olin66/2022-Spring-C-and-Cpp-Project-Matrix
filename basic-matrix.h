#ifndef BASIC_MATRIX_H
#define BASIC_MATRIX_H

#include <string>
#include <cstring>
#include <vector>
#include "matrix.h"
#include "sparse-matrix.h"

template<typename>
class Matrix;

template<typename>
class SparseMatrix;

template<class T>
class BasicMatrix : public Matrix<T> {
private:
    T *data;
public:
    BasicMatrix(int row, int col);

    BasicMatrix(const cv::Mat &mat);

    BasicMatrix(std::vector<std::vector<T>>);

    Matrix<T> &add(const BasicMatrix<T> &);

    Matrix<T> &add(const SparseMatrix<T> &);

    Matrix<T> &subtract(const BasicMatrix<T> &);

    Matrix<T> &subtract(const SparseMatrix<T> &);

    Matrix<T> &scalarMultiply(short);

    Matrix<T> &scalarMultiply(int);

    Matrix<T> &scalarMultiply(long);

    Matrix<T> &scalarMultiply(long long);

    Matrix<T> &scalarMultiply(float);

    Matrix<T> &scalarMultiply(double);

    Matrix<T> &scalarMultiply(long double);

    Matrix<T> &scalarMultiply(std::complex<short>);

    Matrix<T> &scalarMultiply(std::complex<int>);

    Matrix<T> &scalarMultiply(std::complex<long>);

    Matrix<T> &scalarMultiply(std::complex<long long>);

    Matrix<T> &scalarMultiply(std::complex<float>);

    Matrix<T> &scalarMultiply(std::complex<double>);

    Matrix<T> &scalarMultiply(std::complex<long double>);

    Matrix<T> &scalarDivide(short);

    Matrix<T> &scalarDivide(int);

    Matrix<T> &scalarDivide(long);

    Matrix<T> &scalarDivide(long long);

    Matrix<T> &scalarDivide(float);

    Matrix<T> &scalarDivide(double);

    Matrix<T> &scalarDivide(long double);

    Matrix<T> &scalarDivide(std::complex<short>);

    Matrix<T> &scalarDivide(std::complex<int>);

    Matrix<T> &scalarDivide(std::complex<long>);

    Matrix<T> &scalarDivide(std::complex<long long>);

    Matrix<T> &scalarDivide(std::complex<float>);

    Matrix<T> &scalarDivide(std::complex<double>);

    Matrix<T> &scalarDivide(std::complex<long double>);

};

template<class T>
BasicMatrix<T>::BasicMatrix(int row, int col):Matrix<T>(row, col) {
    this->data = new T[row * col];
    std::memset(data, 0, sizeof(T)*row*col);
}

template<class T>
BasicMatrix<T>::BasicMatrix(const cv::Mat &mat):Matrix<T>(mat) {

}

template<class T>
BasicMatrix<T>::BasicMatrix(std::vector<std::vector<T>> mat): Matrix<T>(mat.size(), mat[0].size()) {

}

template<class T>
Matrix<T> &BasicMatrix<T>::add(const BasicMatrix<T> &right) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::add(const SparseMatrix<T> &) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::subtract(const BasicMatrix<T> &) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::subtract(const SparseMatrix<T> &) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(short) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(int) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(long) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(long long int) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(float) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(double) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(long double) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(std::complex<short>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(std::complex<int>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(std::complex<long>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(std::complex<long long int>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(std::complex<float>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(std::complex<double>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarMultiply(std::complex<long double>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(short) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(int) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(long) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(long long int) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(float) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(double) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(long double) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(std::complex<short>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(std::complex<int>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(std::complex<long>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(std::complex<long long int>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(std::complex<float>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(std::complex<double>) {
    return (*this);
}

template<class T>
Matrix<T> &BasicMatrix<T>::scalarDivide(std::complex<long double>) {
    return (*this);
}

#endif