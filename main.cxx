#include <iostream>
#include "Matrix.h"

int main(int argc, char **argv) {
    auto size = (argc == 1) ? 4 : std::atoi(argv[1]);
    Matrix A(size), B(size), C(size);
    A.set(-1.0, 1.0, 5555);
    B.set(-2.0, 2.0, 7777);
    C = 0.0;
    // to have the same order as fortran
    A = A.transpose();
    B = B.transpose();

    // test for pivoting
    for (auto i=0; i<size; i++)
	B.set(i, i, 0.0);

    std::cout << "A:\n";
    A.show();
    std::cout << '\n';
    std::cout << "B:\n";
    B.show();
    std::cout << '\n';
    C += A;
    C -= A;
    std::cout << "C: \n";
    C.show();
    std::cout << '\n';
    C = A*B;
    std::cout << "C = A*B\n";
    C.show();
    std::cout << '\n';

    C /= B;
    std::cout << "C /= B(= A)\n";
    C.show();
    std::cout << '\n';
    if (C == A)
        std::cout << "C == A\n";
    if (A != B)
        std::cout << "A != B\n";
    std::cout << '\n';
    
    auto Ainv = A.inverse();
    std::cout << "Ainv = A^-1\n";
    Ainv.show();
    std::cout << '\n';
    
    auto I = A*Ainv;
    std::cout << "A*Ainv (should equal to unit matrix)\n";
    I.show();
    std::cout << '\n';
    
    auto trA = A.trace();
    std::cout << "tr(A)                                  = " << trA << '\n';
    auto detA = A.determinant();
    std::cout << "det(A)                                 = " << detA << '\n';
    auto detAinv = Ainv.determinant();
    std::cout << "det(Ainv)                              = " << detAinv << '\n';
    std::cout << "det(A)*det(Ainv) (should equal to one) = " << detA*detAinv << '\n';
    std::cout << '\n';
    
    C = C - C;
    std::cout << "trivial test: C = C - C (=0)\n";
    C.show();
    std::cout << '\n';
    
    if (C == 0.0)
        std::cout << "C == 0.0\n";

    Matrix D = 1.0;
    if (D != 2.0)
        std::cout << "D != 2.0\n";
    
    return 0;
}
