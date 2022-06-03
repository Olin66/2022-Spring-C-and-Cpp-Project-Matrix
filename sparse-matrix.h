#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <cstring>
#include <vector>
#include "matrix.h"
#include "basic-matrix.h"

template<typename T>
struct Triple {
    long row;
    long col;
    T val;
};

template <typename >
class Matrix;

template <typename >
class BasicMatrix;

template<class T>
class SparseMatrix : public Matrix<T> {
private:
    std::vector<Triple<T>> triples;
public:
    SparseMatrix(int row, int col);

    SparseMatrix(const cv::Mat &mat);

    SparseMatrix(std::vector<std::vector<T>>);

    Matrix <T> &add(const BasicMatrix<T> &);

    Matrix <T> &add(const SparseMatrix<T> &);

    Matrix <T> &subtract(const BasicMatrix<T> &);

    Matrix <T> &subtract(const SparseMatrix<T> &);

    Matrix <T> &scalarMultiply(short);

    Matrix <T> &scalarMultiply(int);

    Matrix <T> &scalarMultiply(long);

    Matrix <T> &scalarMultiply(long long);

    Matrix <T> &scalarMultiply(float);

    Matrix <T> &scalarMultiply(double);

    Matrix <T> &scalarMultiply(long double);

    Matrix <T> &scalarMultiply(std::complex<short>);

    Matrix <T> &scalarMultiply(std::complex<int>);

    Matrix <T> &scalarMultiply(std::complex<long>);

    Matrix <T> &scalarMultiply(std::complex<long long>);

    Matrix <T> &scalarMultiply(std::complex<float>);

    Matrix <T> &scalarMultiply(std::complex<double>);

    Matrix <T> &scalarMultiply(std::complex<long double>);

    Matrix <T> &scalarDivide(short);

    Matrix <T> &scalarDivide(int);

    Matrix <T> &scalarDivide(long);

    Matrix <T> &scalarDivide(long long);

    Matrix <T> &scalarDivide(float);

    Matrix <T> &scalarDivide(double);

    Matrix <T> &scalarDivide(long double);

    Matrix <T> &scalarDivide(std::complex<short>);

    Matrix <T> &scalarDivide(std::complex<int>);

    Matrix <T> &scalarDivide(std::complex<long>);

    Matrix <T> &scalarDivide(std::complex<long long>);

    Matrix <T> &scalarDivide(std::complex<float>);

    Matrix <T> &scalarDivide(std::complex<double>);

    Matrix <T> &scalarDivide(std::complex<long double>);

};

template<class T>
SparseMatrix<T>::SparseMatrix(int row, int col) {

}

template<class T>
SparseMatrix<T>::SparseMatrix(const cv::Mat &mat) {

}

template<class T>
SparseMatrix<T>::SparseMatrix(std::vector<std::vector<T>>) {

}

template<class T>
Matrix <T> &SparseMatrix<T>::add(const BasicMatrix<T> &) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::add(const SparseMatrix<T> &) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::subtract(const BasicMatrix<T> &) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::subtract(const SparseMatrix<T> &) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(short) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(int) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(long) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(long long int) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(float) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(double) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(long double) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(std::complex<short>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(std::complex<int>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(std::complex<long>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(std::complex<long long int>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(std::complex<float>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(std::complex<double>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarMultiply(std::complex<long double>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(short) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(int) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(long) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(long long int) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(float) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(double) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(long double) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(std::complex<short>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(std::complex<int>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(std::complex<long>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(std::complex<long long int>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(std::complex<float>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(std::complex<double>) {
    return (*this);
}

template<class T>
Matrix <T> &SparseMatrix<T>::scalarDivide(std::complex<long double>) {
    return (*this);
}

#endif