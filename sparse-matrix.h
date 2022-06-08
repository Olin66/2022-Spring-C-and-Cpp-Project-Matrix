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

template<typename>
class Matrix;

template<typename>
class BasicMatrix;

template<class T>
class SparseMatrix : public Matrix<T> {
private:
    std::vector<Triple<T>> triples;
public:
    SparseMatrix(int row, int col);

    SparseMatrix(const cv::Mat &mat);

    SparseMatrix(std::vector<std::vector<T>>);

    void add(const BasicMatrix<T> &);

    void add(const SparseMatrix<T> &);

    void subtract(const BasicMatrix<T> &);

    void subtract(const SparseMatrix<T> &);

    void scalarMultiply(short);

    void scalarMultiply(int);

    void scalarMultiply(long);

    void scalarMultiply(long long);

    void scalarMultiply(float);

    void scalarMultiply(double);

    void scalarMultiply(long double);

    void scalarMultiply(std::complex<short>);

    void scalarMultiply(std::complex<int>);

    void scalarMultiply(std::complex<long>);

    void scalarMultiply(std::complex<long long>);

    void scalarMultiply(std::complex<float>);

    void scalarMultiply(std::complex<double>);

    void scalarMultiply(std::complex<long double>);

    void scalarDivide(short);

    void scalarDivide(int);

    void scalarDivide(long);

    void scalarDivide(long long);

    void scalarDivide(float);

    void scalarDivide(double);

    void scalarDivide(long double);

    void scalarDivide(std::complex<short>);

    void scalarDivide(std::complex<int>);

    void scalarDivide(std::complex<long>);

    void scalarDivide(std::complex<long long>);

    void scalarDivide(std::complex<float>);

    void scalarDivide(std::complex<double>);

    void scalarDivide(std::complex<long double>);

    void dotProduct(const BasicMatrix<T> &);

    void dotProduct(const SparseMatrix<T> &);

    void crossProduct(const BasicMatrix<T> &);

    void crossProduct(const SparseMatrix<T> &);
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
void SparseMatrix<T>::scalarMultiply(short) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(int) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(long) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(long long int) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(float) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(double) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(long double) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(std::complex<short>) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(std::complex<int>) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(std::complex<long>) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(std::complex<long long int>) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(std::complex<float>) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(std::complex<double>) {
}

template<class T>
void SparseMatrix<T>::scalarMultiply(std::complex<long double>) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(short) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(int) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(long) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(long long int) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(float) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(double) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(long double) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(std::complex<short>) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(std::complex<int>) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(std::complex<long>) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(std::complex<long long int>) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(std::complex<float>) {
}

template<class T>
void SparseMatrix<T>::scalarDivide(std::complex<double>) {
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

#endif