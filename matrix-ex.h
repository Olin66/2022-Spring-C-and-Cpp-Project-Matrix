#ifndef MATRIX_EXCEPTION
#define MATRIX_EXCEPTION

#include <exception>
#include <string>
#include "basic-matrix.h"
#include "matrix.h"
#include "sparse-matrix.h"

namespace mat {
    namespace ex {
        class MatrixException : public std::exception {
        public:
            MatrixException() = default;
            virtual const char* what(){
                return "Matrix Exception!";
            }
        };

        class MismatchedSizeException : public MatrixException {
        private:
            int l_row, l_col, r_row, r_col;
            std::string operation;
        public:
            MismatchedSizeException(int l_row, int l_col, int r_row, int r_col, std::string message) : l_row(l_row),
                                                                                                     l_col(l_col),
                                                                                                     r_row(r_row),
                                                                                                     r_col(r_col),
                                                                                                     operation(
                                                                                                             message) {};
            template<class T>                                                                                                 
            MismatchedSizeException(const Matrix<T> &l_mat, const Matrix<T> &r_mat, std::string message): l_row(l_mat.getRow()), l_col(l_mat.getCol()), r_row(r_mat.getRow()), r_col(r_mat.getCol()){};

            template<class T>
            MismatchedSizeException(const BasicMatrix<T> &l_mat, const BasicMatrix<T> &r_mat, std::string message) {
                const Matrix<T>& l = l_mat;
                const Matrix<T>& r = r_mat;
                MismatchedSizeException(l, r, message);
            }

            template<class T>
            MismatchedSizeException(const SparseMatrix<T> &l_mat, const SparseMatrix<T> &r_mat, std::string message) {
                const Matrix<T>& l = l_mat;
                const Matrix<T>& r = r_mat;
                MismatchedSizeException(l, r, message);
            }

            const char *what() override {
                return (new std::string("The matrix exception occurs when doing " + operation + "\nThe sizes " +
                                        std::to_string(l_row) + "*" + std::to_string(l_col) + " and " +
                                        std::to_string(r_row) + "*" + std::to_string(r_col) +
                                        " don't match!"))->c_str();
            }
        };
        class NotSquareException : public MatrixException{
            private:
                int col,row;
                std::string operation;
            public:
            NotSquareException(int row,int col,std::string message) :row(row),col(col),operation(message){};
            template<class T>
            NotSquareException(const BasicMatrix<T> &mat,std::string message):row(mat.getRow()),col(mat.getCol()),operation(message){};
             const char *what() override {
                return (new std::string("The matrix exception occurs when doing " + operation + "\nThe matrix is not square " 
                             
                                       ))->c_str();
            }
        };
        class NoInverseException :public MatrixException{
            private:
                int col,row;
                std::string operation;
            public:
                NoInverseException(int row, int col, std::string message):row(row),col(col),operation(message){};
            template<class T>
            NoInverseException(const BasicMatrix<T> &mat,std::string message):row(mat.getRow()),col(mat.getCol()),operation(message){};
            const char *what() override {
                return (new std::string("The matrix exception occurs when doing " + operation + "\nThe matrix dose not have inverse. " 
                                       ))->c_str();
            }
        };
    }
}

#endif