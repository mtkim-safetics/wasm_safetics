#include "../../include/safetics.hpp"
#include "../external/Eigen/Dense"
#include <cstdio>
using namespace Eigen;
int main(int argc, char const* argv[])
{
    /* code */
    return 0;
}


EM_PORT_API(int) add(int a, int b)
{
    return a + b;
}
EM_PORT_API(Matrix2f) addMatrix2f(Matrix2f a, Matrix2f b)
{
    return a + b;
}
