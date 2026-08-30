#include "ZonoOpt.hpp"
#include "unit_test_utilities.hpp"

using namespace ZonoOpt;

namespace
{
    template <typename SetT>
    void check_form_conversion(SetT Z)
    {
        const SetT original = Z;

        // to_0_1_form: sets flag, idempotent
        Z.to_0_1_form();
        EXPECT_TRUE(Z.is_0_1_form());
        const SetT after_first_call = Z;
        Z.to_0_1_form();
        EXPECT_TRUE(hz_eq(Z, after_first_call)) << "to_0_1_form should be idempotent";

        // to_canonical_form: clears flag, idempotent, recovers original fields
        Z.to_canonical_form();
        EXPECT_FALSE(Z.is_0_1_form());
        const SetT after_second_call = Z;
        Z.to_canonical_form();
        EXPECT_TRUE(hz_eq(Z, after_second_call)) << "to_canonical_form should be idempotent";
        EXPECT_TRUE(hz_eq(Z, original)) << "round trip should recover the original fields";
    }
}

TEST(FormConversion, HybZonoRoundTrip)
{
    std::mt19937 rand_gen(42);
    check_form_conversion(random_hybzono(3, 4, 2, 2, 0.7, -2., 2., rand_gen));
}

TEST(FormConversion, ConZonoRoundTrip)
{
    std::mt19937 rand_gen(42);
    check_form_conversion(random_conzono(3, 4, 2, 0.7, -2., 2., rand_gen));
}

TEST(FormConversion, ZonoRoundTrip)
{
    std::mt19937 rand_gen(42);
    check_form_conversion(random_zono(3, 4, 0.7, -2., 2., rand_gen));
}

TEST(FormConversion, PointIsNoOp)
{
    Eigen::Vector<zono_float, -1> c(2);
    c << 1., 2.;
    Point P(c);

    P.to_0_1_form();
    EXPECT_FALSE(P.is_0_1_form());
    P.to_canonical_form();
    EXPECT_FALSE(P.is_0_1_form());
    EXPECT_TRUE(vector_eq(P.get_c(), c));
}
