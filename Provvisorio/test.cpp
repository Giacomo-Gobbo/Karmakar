#include "LinearConstrainSystem.hpp"

int main()
{
    using namespace boost::numeric::ublas;

    matrix<double> A{identity_matrix<double>(2)};
    A(0, 1) = -1;
    A(1, 0) = 1;
    std::cout << "A:" << A << std::endl;
    vector<double> b(2, 2);
    b(1) = 6;
    vector<double> c(2, 1);
    c(1) = 0.5;
    vector<double> x0(2, 0.5);
    x0(1) = 0;

    LinearConstrainSystem<double>::ConstrainType cstType{LinearConstrainSystem<double>::ConstrainType::LE};

    LinearConstrainSystem<double> lcs(A, b, cstType);

    if (lcs.karmakar(x0, c, 0.5, LinearConstrainSystem<double>::OptimizationType::MAX) == LinearConstrainSystem<double>::SolutionType::BOUNDED)
    {
        std::cout << "BOUNDED" << std::endl;
    } else {
        std::cout << "UNBOUNDED" << std::endl;   
    }

    return 0;
}