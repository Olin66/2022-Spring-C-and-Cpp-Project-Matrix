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
void BasicMatrix<T>::add(const BasicMatrix<T> &right) {
    for (size_t i = 0;i < this->getSize();i++){
        this->data[i] = this->data[i] + right.data[i];
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
void BasicMatrix<T>::scalarMultiply(short) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(int) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(long) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(long long int) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(float) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(double) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(long double) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(std::complex<short>) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(std::complex<int>) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(std::complex<long>) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(std::complex<long long int>) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(std::complex<float>) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(std::complex<double>) {
}

template<class T>
void BasicMatrix<T>::scalarMultiply(std::complex<long double>) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(short) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(int) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(long) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(long long int) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(float) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(double) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(long double) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(std::complex<short>) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(std::complex<int>) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(std::complex<long>) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(std::complex<long long int>) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(std::complex<float>) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(std::complex<double>) {
}

template<class T>
void BasicMatrix<T>::scalarDivide(std::complex<long double>) {
}

#endif