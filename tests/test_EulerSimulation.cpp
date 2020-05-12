#define _USE_MATH_DEFINES

#include <iostream>
#include "gtest/gtest.h"
#include "../TH_CS395T_final/lib/EulerSimulation.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

namespace FluidSimulation
{
  class EulerSimulationTest : public ::testing::Test {
  protected:
    EulerSimulationTest() : simulation(Eigen::Vector3i(4,4,4), Eigen::Vector3d(0.1, 0.1, 0.1))
    {

    }

    void SetUp() override {

    }

    EulerSimulation simulation;
  };

  class EulerSimulationLevelSetTest : public ::testing::Test {
  protected:
    EulerSimulationLevelSetTest() : simulation(Eigen::Vector3i(4, 4, 4), Eigen::Vector3d(0.1, 0.1, 0.1))
    {
      Eigen::Vector3d sphereCenter(0.25, 0.25, 0.25);
      double radius = 0.05;
      for (int x = 0; x < simulation.m_currentState.m_dims(X)+1; x++)
      {
        for (int y = 0; y < simulation.m_currentState.m_dims(Y)+1; y++)
        {
          for (int z = 0; z < simulation.m_currentState.m_dims(Z)+1; z++)
          {
            int i = z + y * (simulation.m_currentState.m_dims(Z) + 1) + x * (simulation.m_currentState.m_dims(Z) + 1) * (simulation.m_currentState.m_dims(Y) + 1);
            Eigen::Vector3d position(x + 0.5, y + 0.5, z + 0.5);
            position = position.cwiseProduct(simulation.m_currentState.m_gridSizeHorizontal);
            Eigen::Vector3d direction = (position - sphereCenter);
            double distance = (direction).norm() - radius;

            // Set component 0 to the distance
            simulation.m_currentState.m_signedDistance(i) = distance;
            // Set components 1-3 to the direction
            if (direction.norm() >= 0)
            {
              simulation.m_currentState.m_dSignedDistance.row(i) = direction.transpose().normalized();
            }
            else
            {
              // If it's equidistant than pick a random direction
              simulation.m_currentState.m_dSignedDistance.row(i) = Eigen::Vector3d::Random().transpose();
            }
          }
        }
      }
      //simulation.m_currentState.m_velocity.insert(3 + 5 * 3 + 5 * 3 * 3, 2) = -1;
      //simulation.m_currentState.m_velocity.insert(2 + 5 * 3 + 5 * 3 * 3, 2) = -1;

      spdlog::info("Signed Distance = \n{}", simulation.m_currentState.m_signedDistance);
    }

    void SetUp() override {

    }

    EulerSimulation simulation;
  };

  TEST_F(EulerSimulationLevelSetTest, SphereAdvection)
  {
    simulation.m_enablePressure = false;
    simulation.m_enableGravity = false;
    spdlog::set_level(spdlog::level::debug);
    spdlog::info("Running SphereAdvection");

    simulation.m_currentState.m_velocity.insert(3 + 5 * 3 + 5 * 3 * 3, 2) = -1;
    simulation.m_currentState.m_velocity.insert(2 + 5 * 3 + 5 * 3 * 3, 2) = -1;

    // 0.35
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_signedDistance(3 + 5 * 2 + 5 * 5 * 2), 0.05);
    // 0.25
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_signedDistance(2 + 5 * 2 + 5 * 5 * 2), -0.05);
    // 0.15
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_signedDistance(1 + 5 * 2 + 5 * 5 * 2), 0.05);
    // 0.05
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_signedDistance(0 + 5 * 2 + 5 * 5 * 2), 0.15);

    simulation.step(0.05);
    spdlog::info("Stepped");
    spdlog::info("New Signed Distance = \n{}", simulation.m_currentState.m_signedDistance);
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_signedDistance(3 + 5 * 2 + 5 * 5 * 2), 0.1);
    EXPECT_NEAR(simulation.m_currentState.m_signedDistance(2 + 5 * 2 + 5 * 5 * 2), 0, 1e-9);
    EXPECT_NEAR(simulation.m_currentState.m_signedDistance(1 + 5 * 2 + 5 * 5 * 2), 0, 1e-9);
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_signedDistance(0 + 5 * 2 + 5 * 5 * 2), 0.1);
  }

  TEST_F(EulerSimulationLevelSetTest, ConstantVelocityAdvection)
  {
    simulation.m_enablePressure = false;
    simulation.m_enableGravity = false;
    spdlog::set_level(spdlog::level::debug);
    spdlog::info("Running ConstantVelocityAdvection");

    simulation.m_currentState.m_velocity.insert(3 + 5 * 3 + 5 * 3 * 3, 2) = -1;
    simulation.m_currentState.m_velocity.insert(2 + 5 * 3 + 5 * 3 * 3, 2) = -1;

    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(3 + 5 * 2 + 5 * 5 * 2, 2), -1);
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(2 + 5 * 2 + 5 * 5 * 2, 2), -1);
    simulation.step(0.05);
    spdlog::info("Stepped");
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(2 + 5 * 2 + 5 * 5 * 2, 2), -1);
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(1 + 5 * 2 + 5 * 5 * 2, 2), -1);
    //EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(2, 0), -0.5);
    //EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(2, 0), -1.0);
  }

  TEST_F(EulerSimulationLevelSetTest, GravityEnabled)
  {
    simulation.m_enablePressure = false;
    simulation.m_enableGravity = true;
    spdlog::set_level(spdlog::level::debug);
    spdlog::info("Running ConstantVelocityAdvection");

    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(3 + 5 * 2 + 5 * 5 * 2, 2), 0);
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(2 + 5 * 2 + 5 * 5 * 2, 2), 0);
    simulation.step(0.05);
    spdlog::info("Stepped");
    ASSERT_LT(simulation.m_currentState.m_velocity.coeff(2 + 5 * 2 + 5 * 5 * 2, 2), 0);
    ASSERT_LT(simulation.m_currentState.m_velocity.coeff(1 + 5 * 2 + 5 * 5 * 2, 2), 0);

    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(2 + 5 * 2 + 5 * 5 * 2, 2), -0.05 * 9.8);
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(1 + 5 * 2 + 5 * 5 * 2, 2), -0.05 * 9.8);

    // 0.35
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_signedDistance(3 + 5 * 2 + 5 * 5 * 2), 0.05);
    // 0.25
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_signedDistance(2 + 5 * 2 + 5 * 5 * 2), -0.05);
    // 0.15
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_signedDistance(1 + 5 * 2 + 5 * 5 * 2), 0.05);
    // 0.05
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_signedDistance(0 + 5 * 2 + 5 * 5 * 2), 0.15);

    simulation.step(0.05);
    ASSERT_LT(simulation.m_currentState.m_velocity.coeff(2 + 5 * 2 + 5 * 5 * 2, 2), 0);
    ASSERT_LT(simulation.m_currentState.m_velocity.coeff(1 + 5 * 2 + 5 * 5 * 2, 2), 0);
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(2 + 5 * 2 + 5 * 5 * 2, 2), -0.1 * 9.8);
    EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(1 + 5 * 2 + 5 * 5 * 2, 2), -0.1 * 9.8);



    //EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(2, 0), -0.5);
    //EXPECT_DOUBLE_EQ(simulation.m_currentState.m_velocity.coeff(2, 0), -1.0);
  }
}
