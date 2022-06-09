#ifndef MATRIX_EXCEPTION
#define MATRIX_EXCEPTION

#include <exception>
#include <string>

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

            const char *what() override {
                return (new std::string("The matrix exception occurs when doing " + operation + "\nThe sizes " +
                                        std::to_string(l_row) + "*" + std::to_string(l_col) + " and " +
                                        std::to_string(r_row) + "*" + std::to_string(r_col) +
                                        " don't match!"))->c_str();
            }
        };
    }
}

#endif