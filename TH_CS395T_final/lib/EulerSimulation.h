#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "EulerState.h"
#include "spdlog/fmt/ostr.h" // must be included

namespace FluidSimulation
{
	static Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	static Eigen::IOFormat LongFmt(3, 0, ", ", "; ", "[", "]");

	class EulerSimulation
	{
	public:
		EulerSimulation(Eigen::Vector3i dimSize, Eigen::Vector3d gridWidth) : 
			m_density(1.0), m_gravity(-9.8), m_currentState(dimSize, gridWidth) { };
		void step(double h);
	private:
		void getAdvectedVelocityField(double h, Eigen::SparseMatrix<double>& velocityField,
			Eigen::SparseMatrix<double>& interpolateX, Eigen::SparseMatrix<double>& interpolateY, Eigen::SparseMatrix<double>& interpolateZ,
			Eigen::SparseMatrix<double, Eigen::ColMajor>& oldVelocityMid);
		void advectSignedDistance(double h, const Eigen::SparseMatrix<double, Eigen::ColMajor> oldVelocityMid);

		// TODO: Calculate timestep from formula
		double getTimeStep() { return 0.001; };

		// State of current fluid simulation
		EulerState m_currentState;

		// Position of fluid box (0,0,0) in world coordinates
		Eigen::Vector3d m_position;
		// Axis-angle rotation of fluid box
		Eigen::Vector3d m_theta;

		// Assume incompressible fluid and constant temperature/concentration - density is constant
		double m_density;
		double m_gravity;
		bool m_enablePressure;

		friend class EulerSimulationLevelSetTest;
		friend class EulerSimulationTest_TestConstantAdvection_Test;
		friend class EulerSimulationLevelSetTest_SphereAdvection_Test;
		friend class EulerSimulationLevelSetTest_ConstantVelocityAdvection_Test;
	};

}