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
        T *data;
    public:
        BasicMatrix(int row, int col);

        // BasicMatrix(const cv::Mat &mat);

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

        T getTrace();

        T getDeterminant();

        void reshape(int row, int col);

        void sliceRow(int row1, int row2);

        void sliceCol(int col1, int col2);

        void slice(int row1, int row2, int col1, int col2);

        Matrix<T> &convolve(BasicMatrix<T> &);

        Matrix<T> &convolve(SparseMatrix<T> &);

    };

    template<class T>
    BasicMatrix<T>::BasicMatrix(int row, int col):Matrix<T>(row, col) {
        this->data = new T[row * col];
        std::memset(data, 0, sizeof(T) * row * col);
    }

    // template<class T>
    // BasicMatrix<T>::BasicMatrix(const cv::Mat &mat):Matrix<T>(mat) {

    // }

    template<class T>
    BasicMatrix<T>::BasicMatrix(std::vector<std::vector<T>> mat): Matrix<T>(mat.size(), mat[0].size()) {
        for (size_t i = 0; i < this->getSize(); i++)
        {
            // this->data[i] = mat[i/this->col][i%this->col];
        }
        
        
    }

    template<class T>
    void BasicMatrix<T>::add(const BasicMatrix<T> &right) {
        if (this->getRow() != right.getRow() || this->getCol() != right.getCol())
            throw ex::MismatchedSizeException(this->getRow(), this->getCol(), right.getRow(), right.getCol(), "matrix addition");
        for (size_t i = 0; i < this->getSize(); i++) {
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

    template<class T>
    void BasicMatrix<T>::dotProduct(const BasicMatrix<T> &) {

    }

    template<class T>
    void BasicMatrix<T>::dotProduct(const SparseMatrix<T> &) {

    }

    template<class T>
    void BasicMatrix<T>::crossProduct(const BasicMatrix<T> &) {

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
    }

    template<class T>
    T BasicMatrix<T>::getMin() {
    }

    template<class T>
    T BasicMatrix<T>::getSum() {
    }

    template<class T>
    T BasicMatrix<T>::getAvg() {
    }

    template<class T>
    T BasicMatrix<T>::getEigenvalue() {
    }

    template<class T>
    T BasicMatrix<T>::getRowMin(int) {
    }

    template<class T>
    T BasicMatrix<T>::getColMin(int) {
    }

    template<class T>
    T BasicMatrix<T>::getRowMax(int) {
    }

    template<class T>
    T BasicMatrix<T>::getColMax(int) {
    }

    template<class T>
    T BasicMatrix<T>::getRowSum(int) {
    }

    template<class T>
    T BasicMatrix<T>::getColSum(int) {
    }

    template<class T>
    T BasicMatrix<T>::getRowAvg(int) {
    }

    template<class T>
    T BasicMatrix<T>::getColAvg(int) {
    }

    template<class T>
    Matrix<T> &BasicMatrix<T>::getEigenvector() {
    }

    template<class T>
    T BasicMatrix<T>::getTrace() {
    }

    template<class T>
    T BasicMatrix<T>::getDeterminant() {
    }

    template<class T>
    void BasicMatrix<T>::reshape(int row, int col) {

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

}

#endif