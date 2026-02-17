#include <gtest/gtest.h>
#include <cmath>
#include "../src/Polinoms.h"

int32_t make_ijk(int i, int j, int k) {
    return i | (j << 10) | (k << 20);
}

TEST(PolinomTest, CreationAndOutput) {
    int32_t deg = make_ijk(2, 3, 1);
    polijop p(vector<pint>{{deg, 2.5}});

    std::stringstream ss;
    ss << p;
    EXPECT_TRUE(ss.str().find("2.5x^2y^3z") != std::string::npos);
}

TEST(PolinomTest, Addition) {
    polijop p1(vector<pint>{{make_ijk(1, 0, 0), 1.0}});
    polijop p2(vector<pint>{{make_ijk(1, 0, 0), 2.0}, {0, 5.0}});

    polijop res = p1 + p2;
    polijop expected(vector<pint>{{make_ijk(1, 0, 0), 3.0}, {0, 5.0}});

    EXPECT_EQ(res, expected);

    polijop p3 = p1 + 10.0;
    EXPECT_EQ(p3, polijop(vector<pint>{{make_ijk(1, 0, 0), 1.0}, {0, 10.0}}));
}

TEST(PolinomTest, Multiplication) {
    polijop p1(vector<pint>{{make_ijk(1, 1, 0), 1.0}});
    polijop p2(vector<pint>{{make_ijk(1, 0, 1), 2.0}});

    polijop res = p1 * p2;
    polijop expected(vector<pint>{{make_ijk(2, 1, 1), 2.0}});

    EXPECT_EQ(res, expected);

    polijop zero = p1 * 0.0;
    EXPECT_EQ(zero, polijop(vector<pint>{}));
}

TEST(PolinomTest, Overflow) {
    polijop p1(vector<pint>{{make_ijk(1000, 0, 0), 1.0}});
    polijop p2(vector<pint>{{make_ijk(30, 0, 0), 1.0}});

    EXPECT_THROW(p1 * p2, std::overflow_error);
}

TEST(PolinomTest, SubtractionAndDivision) {
    polijop p1(vector<pint>{{0, 10.0}});
    polijop res = p1 - 4.0;
    EXPECT_EQ(res, polijop(vector<pint>{{0, 6.0}}));

    polijop p2(vector<pint>{{make_ijk(1, 1, 1), 10.0}});
    polijop div = p2 / 2.0;
    EXPECT_EQ(div, polijop(vector<pint>{{make_ijk(1, 1, 1), 5.0}}));

    EXPECT_THROW(p2 / 0.0, std::runtime_error);
}

TEST(PolinomTest, Compression) {
    polijop p1(vector<pint>{{make_ijk(1, 0, 0), 1.0}});
    polijop p2(vector<pint>{{make_ijk(1, 0, 0), -1.0}});
    polijop res = p1 + p2;

    polijop zero(vector<pint>{});
    EXPECT_EQ(res, zero);
}


TEST(PolinomAdvanced, Distributivity) {
    polijop a(vector<pint>{{make_ijk(1, 0, 0), 2.0}});
    polijop b(vector<pint>{{make_ijk(0, 1, 0), 1.0}});
    polijop c(vector<pint>{{0, 3.0}});

    polijop left = a * (b + c);
    polijop right = (a * b) + (a * c);

    EXPECT_EQ(left, right);
}

TEST(PolinomAdvanced, DifferenceOfSquares) {
    polijop x(vector<pint>{{make_ijk(1, 0, 0), 1.0}});
    polijop y(vector<pint>{{make_ijk(0, 1, 0), 1.0}});

    polijop res = (x - y) * (x + y);
    polijop expected(vector<pint>{{make_ijk(2, 0, 0), 1.0}, {make_ijk(0, 2, 0), -1.0}});

    EXPECT_EQ(res, expected);
}

TEST(PolinomAdvanced, ChainedOpsWithConstants) {
    polijop p(vector<pint>{{make_ijk(0, 0, 1), 1.0}});
    polijop res = ((p + 2.0) * 3.0 - 6.0) / 3.0;

    EXPECT_EQ(res, p);
}

TEST(PolinomAdvanced, TripleVariables) {
    polijop p(vector<pint>{
        {make_ijk(1, 0, 0), 1.0},
        {make_ijk(0, 1, 0), 1.0},
        {make_ijk(0, 0, 1), 1.0}
    });

    polijop res = p * p;
    polijop expected(vector<pint>{
        {make_ijk(2, 0, 0), 1.0},
        {make_ijk(0, 2, 0), 1.0},
        {make_ijk(0, 0, 2), 1.0},
        {make_ijk(1, 1, 0), 2.0},
        {make_ijk(1, 0, 1), 2.0},
        {make_ijk(0, 1, 1), 2.0}
    });

    EXPECT_EQ(res, expected);
}

TEST(PolinomAdvanced, ComplexZero) {
    polijop p(vector<pint>{{make_ijk(5, 5, 5), 10.0}});
    polijop res = (p * 0.5) - (p * 0.5);

    EXPECT_EQ(res, polijop(vector<pint>{}));
}

TEST(PolinomAdvanced, OrderCheck) {
    polijop p1(vector<pint>{{make_ijk(1, 0, 0), 1.0}});
    polijop p2(vector<pint>{{make_ijk(2, 0, 0), 1.0}});
    polijop p3(vector<pint>{{0, 1.0}});

    polijop sum = p1 + p2 + p3;
    polijop expected(vector<pint>{
        {make_ijk(2, 0, 0), 1.0},
        {make_ijk(1, 0, 0), 1.0},
        {0, 1.0}
    });

    EXPECT_EQ(sum, expected);
}

TEST(PolinomAdvanced, Commutativity) {
    polijop p1(vector<pint>{{make_ijk(1, 2, 3), 4.5}});
    polijop p2(vector<pint>{{make_ijk(3, 2, 1), 1.2}});

    EXPECT_EQ(p1 + p2, p2 + p1);
    EXPECT_EQ(p1 * p2, p2 * p1);
}

TEST(PolinomAdvanced, ScalarMultiplication) {
    polijop p(vector<pint>{{make_ijk(1, 1, 1), 1.0}});
    EXPECT_EQ(p * 2.5, 2.5 * p);
    EXPECT_EQ(p * 0.0, polijop(vector<pint>{}));
}