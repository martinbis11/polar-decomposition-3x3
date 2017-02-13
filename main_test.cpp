// This is a example test of the polar decomposition algorithm implementation.
#include "include/polar_decomposition_3x3.h"

#include <iostream>


namespace polar
{


    namespace detail
    {


        template < typename TReal >
        void multiply(
            matrix<TReal, 3, 3>& result,
            const matrix<TReal, 3, 3>& a,
            const matrix<TReal, 3, 3>& b
            )
        {
            result(0,0) = a(0,0) * b(0,0) + a(1,0) * b(0,1) + a(2,0) * b(0,2);
            result(0,1) = a(0,1) * b(0,0) + a(1,1) * b(0,1) + a(2,1) * b(0,2);
            result(0,2) = a(0,2) * b(0,0) + a(1,2) * b(0,1) + a(2,2) * b(0,2);
            result(1,0) = a(0,0) * b(1,0) + a(1,0) * b(1,1) + a(2,0) * b(1,2);
            result(1,1) = a(0,1) * b(1,0) + a(1,1) * b(1,1) + a(2,1) * b(1,2);
            result(1,2) = a(0,2) * b(1,0) + a(1,2) * b(1,1) + a(2,2) * b(1,2);
            result(2,0) = a(0,0) * b(2,0) + a(1,0) * b(2,1) + a(2,0) * b(2,2);
            result(2,1) = a(0,1) * b(2,0) + a(1,1) * b(2,1) + a(2,1) * b(2,2);
            result(2,2) = a(0,2) * b(2,0) + a(1,2) * b(2,1) + a(2,2) * b(2,2);
        }


        template < typename TReal >
        void subtract(
            matrix<TReal, 3, 3>& result,
            const matrix<TReal, 3, 3>& a,
            const matrix<TReal, 3, 3>& b
            )
        {
            result(0) = a(0) - b(0);
            result(1) = a(1) - b(1);
            result(2) = a(2) - b(2);
            result(3) = a(3) - b(3);
            result(4) = a(4) - b(4);
            result(5) = a(5) - b(5);
            result(6) = a(6) - b(6);
            result(7) = a(7) - b(7);
            result(8) = a(8) - b(8);
        }


        template < typename TReal >
        inline TReal norm(const matrix<TReal, 3, 3>& m)
        {
            TReal length = m(0) * m(0);
            length += m(1) * m(1);
            length += m(2) * m(2);
            length += m(3) * m(3);
            length += m(4) * m(4);
            length += m(5) * m(5);
            length += m(6) * m(6);
            length += m(7) * m(7);
            length += m(8) * m(8);

            return math_utils<TReal>::sqrt(length);
        }


        template < typename TReal >
        inline matrix<TReal, 3, 3> identity()
        {
            matrix<TReal, 3, 3> m;
            m(0) = 1;
            m(1) = 0;
            m(2) = 0;
            m(3) = 0;
            m(4) = 1;
            m(5) = 0;
            m(6) = 0;
            m(7) = 0;
            m(8) = 1;
            return m;
        }


    }; // End of namespace detail.


}; // End of namespace polar.




namespace
{


    // These methods validate the output of the algorithm.
    template <typename TReal>
    void Check(
        const TReal valuesA[],
        const TReal valuesQ[],
        const TReal valuesH[]
        )
    {
        const polar::detail::matrix<TReal, 3, 3>& A = *reinterpret_cast<const polar::detail::matrix<TReal, 3, 3>*>(valuesA);
        const polar::detail::matrix<TReal, 3, 3>& Q = *reinterpret_cast<const polar::detail::matrix<TReal, 3, 3>*>(valuesQ);
        const polar::detail::matrix<TReal, 3, 3>& H = *reinterpret_cast<const polar::detail::matrix<TReal, 3, 3>*>(valuesH);

        polar::detail::matrix<TReal, 3, 3> temp;
        polar::detail::multiply(temp, Q, H);
        polar::detail::subtract(temp, A, temp);
        const TReal residual = polar::detail::norm(temp) / (polar::detail::norm(A) + std::numeric_limits<TReal>::min());

        polar::detail::transpose_multiply(temp, Q, Q);
        polar::detail::subtract(temp, temp, polar::detail::identity<TReal>());
        const TReal orthogonality = norm(temp);

        std::cout << "Relative residual = " << residual << " ,  orthogonality = " << orthogonality << std::endl;
    }


    template <typename TReal>
    void CheckExactQ(
        const TReal y,
        const TReal valuesQ[]
        )
    {
        const polar::detail::matrix<TReal, 3, 3>& Q = *reinterpret_cast<const polar::detail::matrix<TReal, 3, 3>*>(valuesQ);

        static const TReal valuesQ1[] = {
            139 / static_cast<TReal>(255), 466 / static_cast<TReal>(1275), 962 / static_cast<TReal>(1275),
            -14 / static_cast<TReal>(51), -197 / static_cast<TReal>(255), 146 / static_cast<TReal>(255),
            202 / static_cast<TReal>(255), -662 / static_cast<TReal>(1275), -409 / static_cast<TReal>(1275),
        };
        const polar::detail::matrix<TReal, 3, 3>& Q1 = *reinterpret_cast<const polar::detail::matrix<TReal, 3, 3>*>(valuesQ1);
        // Condition number of U.
        const TReal kappa = polar::detail::math_utils<TReal>::sqrt((1 + 2 * y * y) / (3 * y * y));

        polar::detail::matrix<TReal, 3, 3> temp;
        polar::detail::subtract(temp, Q, Q1);
        const TReal error = polar::detail::norm(temp) / (polar::detail::norm(Q1) * kappa);
        std::cout << "Condition number of U = " << kappa << ",  scaled relative error in Q = " << error << std::endl;
    }


    template <typename TReal>
    void RunTest()
    {
        // This implementation is column-major.
        std::cout << "Test (5.1) from paper:" << std::endl;
        {
            const TReal A[] = {
                static_cast<TReal>(0.1), static_cast<TReal>(0.1), static_cast<TReal>(0.3),
                static_cast<TReal>(0.2), static_cast<TReal>(0.1), static_cast<TReal>(0.2),
                static_cast<TReal>(0.3), static_cast<TReal>(0.0), static_cast<TReal>(0.1),
            };
            TReal Q[9];
            TReal H[9];
            polar::polar_decomposition(Q, H, A);
            Check(A, Q, H);
        }

        std::cout << "Test (5.2) from paper:" << std::endl;
        static const TReal yValues[] = {
            static_cast<TReal>(1),
            static_cast<TReal>(1.0e-4),
            static_cast<TReal>(1.0e-8),
            static_cast<TReal>(1.0e-12),
            static_cast<TReal>(1.0e-16),
        };
        static const int yValuesCount = sizeof(yValues) / sizeof(yValues[0]);
        for (int i = 0; i < yValuesCount; ++i)
        {
            const TReal y = yValues[i];
            std::cout << i << " : y = " << y << std::endl;

            const TReal A[] = {
                (720*y - 25) / 1275, (396*y + 70) / 1275, (972*y - 10) / 1275,
                (-650*y + 300) / 1275, (-145*y - 840) / 1275, (610*y + 120) / 1275,
                (710*y + 300) / 1275, (178*y - 840) / 1275, (-529*y + 120) / 1275,
            };
            TReal Q[9];
            TReal H[9];
            polar::polar_decomposition(Q, H, A);
            Check(A, Q, H);
            CheckExactQ(y, Q);
        }
    }


}; // End of anonymous namespace.




int main(void)
{
    RunTest<float>();
    RunTest<double>();

    return 0;
}
