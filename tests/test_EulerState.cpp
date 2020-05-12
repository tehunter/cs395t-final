#define _USE_MATH_DEFINES

#include "gtest/gtest.h"
#include "../TH_CS395T_final/lib/EulerState.h"
#include <iostream>

TEST(basic_tests, TestCompiles) {
  int a = 1;
  EXPECT_EQ(a, 1);
}

namespace FluidSimulation
{
  class EulerStateStencilTest : public ::testing::Test {
  protected:
    EulerStateStencilTest() : state(Eigen::Vector3i(1, 1, 4), Eigen::Vector3d(0.1, 0.1, 0.1)) { };

    void SetUp() override {

    }

    EulerState state;
  };

  TEST_F(EulerStateStencilTest, StencilZ)
  {
    Eigen::SparseMatrix<double, Eigen::ColMajor> stencil;

    std::cout << stencil.size() << std::endl;

    state.getCentralDifferenceStencil(Z, true, stencil);
    std::cout << stencil.size() << std::endl;

    EXPECT_DOUBLE_EQ(stencil.coeff(0, 0), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 1), 1);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 2), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 3), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 4), 0);

    EXPECT_DOUBLE_EQ(stencil.coeff(1, 0), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 1), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 2), 1);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 3), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 4), 0);

    EXPECT_DOUBLE_EQ(stencil.coeff(2, 0), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 1), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 2), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 3), 1);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 4), 0);

    EXPECT_DOUBLE_EQ(stencil.coeff(3, 0), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 1), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 2), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 3), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 4), 1);
  }

  TEST_F(EulerStateStencilTest, StencilY)
  {
    Eigen::SparseMatrix<double, Eigen::ColMajor> stencil;

    std::cout << stencil.size() << std::endl;

    state.getCentralDifferenceStencil(Y, true, stencil);
    std::cout << stencil.size() << std::endl;

    EXPECT_DOUBLE_EQ(stencil.coeff(0, 0), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 1), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 2), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 3), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 4), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 5), 1);

    EXPECT_DOUBLE_EQ(stencil.coeff(1, 0), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 1), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 2), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 3), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 4), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 5), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 6), 1);

    EXPECT_DOUBLE_EQ(stencil.coeff(2, 0), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 1), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 2), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 3), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 4), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 7), 1);

    EXPECT_DOUBLE_EQ(stencil.coeff(3, 0), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 1), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 2), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 3), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 8), 1);
  }

  TEST_F(EulerStateStencilTest, StencilX)
  {
    Eigen::SparseMatrix<double, Eigen::ColMajor> stencil;

    std::cout << stencil.size() << std::endl;

    state.getCentralDifferenceStencil(X, true, stencil);
    std::cout << stencil.size() << std::endl;

    EXPECT_DOUBLE_EQ(stencil.coeff(0, 0), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 1), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 2), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 3), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 4), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 5), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(0, 10), 1);

    EXPECT_DOUBLE_EQ(stencil.coeff(1, 0), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 1), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 2), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 3), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 4), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 5), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(1, 11), 1);

    EXPECT_DOUBLE_EQ(stencil.coeff(2, 0), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 1), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 2), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 3), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 4), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(2, 12), 1);

    EXPECT_DOUBLE_EQ(stencil.coeff(3, 0), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 1), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 2), 0);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 3), -1);
    EXPECT_DOUBLE_EQ(stencil.coeff(3, 13), 1);
  }

  class EulerStateVelocityGradientTest : public ::testing::Test {
  protected:
    EulerStateVelocityGradientTest() : state(Eigen::Vector3i(1, 1, 4), Eigen::Vector3d(0.1, 0.1, 0.1)) { };

    void SetUp() override {
      //Eigen::SparseMatrix<double> v = state.getVelocityGrid();
      state.m_velocity.resize(state.getGridMatrixSize(true), 3);

      state.m_velocity.insert(1, 2) = 1;
      state.m_velocity.insert(2, 2) = 2;

      state.m_velocity.insert(0, 1) = 1;
      state.m_velocity.insert(5, 1) = 3;

      state.m_velocity.insert(0, 0) = 1;
      state.m_velocity.insert(10, 0) = 4;

      std::cout << "v updated\n";
      std::cout << state.m_velocity << "\n";
    }

    EulerState state;
  };

  TEST_F(EulerStateVelocityGradientTest, VelocityGradient)
  {
    Eigen::SparseMatrix<double, Eigen::ColMajor> dv;
    state.getQuantityGradient(dv, true, state.m_velocity);

    EXPECT_EQ(dv.cols(), 3);
    EXPECT_EQ(dv.rows(), 4);

    EXPECT_DOUBLE_EQ(dv.coeff(0, 2), 1 /0.1);
    EXPECT_DOUBLE_EQ(dv.coeff(1, 2), 1 / 0.1);
    EXPECT_DOUBLE_EQ(dv.coeff(2, 2), -2 / 0.1);
    EXPECT_DOUBLE_EQ(dv.coeff(3, 2), 0);

    EXPECT_DOUBLE_EQ(dv.coeff(0, 1), 2 / 0.1);
    EXPECT_DOUBLE_EQ(dv.coeff(1, 1), 0);
    EXPECT_DOUBLE_EQ(dv.coeff(2, 1), 0);
    EXPECT_DOUBLE_EQ(dv.coeff(3, 1), 0);

    EXPECT_DOUBLE_EQ(dv.coeff(0, 0), 3 / 0.1);
    EXPECT_DOUBLE_EQ(dv.coeff(1, 0), 0);
    EXPECT_DOUBLE_EQ(dv.coeff(2, 0), 0);
    EXPECT_DOUBLE_EQ(dv.coeff(3, 0), 0);
  }
}
