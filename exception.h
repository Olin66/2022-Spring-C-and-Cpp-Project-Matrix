#ifndef MATRIX_EXCEPTION
#define MATRIX_EXCEPTION
#include <exception>

namespace mat{
    namespace exception{
        class MatrixException: public std::exception
        {
        public:
            MatrixException(){}
        };
    }
}

#endif