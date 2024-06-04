#include "../../include/safetics.hpp"
#include "../external/Eigen/Dense"
#include "../external/Eigen/Eigenvalues"
#include <cmath>
#include <cstdio>
using namespace Eigen;
namespace Safetics {
    // Define the Matrix2f type
    typedef Matrix<float, 2, 2> Matrix2f;
    // Define the MatrixXf type
    typedef Matrix<float, Dynamic, Dynamic> MatrixXf;
    int main(int argc, char const* argv[]) {
    /* code */
    return 0;
}
EM_PORT_API(int)
add(int a, int b) {
    return a + b;
}
EM_PORT_API(Matrix2f)
addMatrix2f(Matrix2f a, Matrix2f b) {
    return a + b;
}
// fractional_matrix_power 함수 정의
EM_PORT_API(MatrixXf)
fractionalMatrixPower(const MatrixXf& matrix, float power) {
    // ComplexSchur<MatrixXf> schur(matrix);
    // if (schur.info() != Success) {
    //     // Handle error
    //     std::printf("Schur decomposition failed.\n");
    //     return MatrixXf::Identity(matrix.rows(), matrix.cols());
    // }

    // MatrixXcf T = schur.matrixT();
    // MatrixXcf U = schur.matrixU();
    // for (int i = 0; i < T.rows(); ++i) {
    //     T(i, i) = std::pow(T(i, i), power);
    // }

    // MatrixXcf result = U * T * U.adjoint();
    // if ((result.imag().array() != 0).any()) {
    //     // If any element is complex, return identity matrix
    //     return MatrixXf::Identity(matrix.rows(), matrix.cols());
    // }
    // return result.real();
    return MatrixXf::Identity(matrix.rows(), matrix.cols());

}
}
