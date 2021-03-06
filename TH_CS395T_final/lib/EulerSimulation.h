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

	enum Solver { CG, Simplicial, CholMod };

	class EulerSimulation
	{
	public:
		EulerSimulation(Eigen::Vector3i dimSize, Eigen::Vector3d gridWidth) :
			m_density(1.0), m_gravity(-9.8), m_currentState(dimSize, gridWidth),
			m_position(0, 0, 0), m_theta(0, 0, 0)
		{
			m_enableGravity = true;
			m_enablePressure = true;
			m_recalculateLevelSet = false;
			m_hasTrackerParticles = false;
			m_cutoffDistance = 4.0;
			m_solver = Simplicial;
		};
		double step(double h);
		void reset();
		void reset(EulerState state);
		void addTrackerParticles();
		EulerState const* getState() { return &m_currentState; };
		EulerState* updateState() { return &m_currentState; };

		bool getRecalculateLevelSet() { return m_recalculateLevelSet; };
		bool getGravityEnabled() { return m_enableGravity; };
		bool getPressureEnabled() { return m_enablePressure; };
		void toggleGravity() { m_enableGravity = !m_enableGravity; };
		void togglePressure() { m_enablePressure = !m_enablePressure; };
		void toggleRecalculateLevelSet() { m_recalculateLevelSet = !m_recalculateLevelSet; };

		double m_cutoffDistance;
		Solver m_solver;

		Eigen::MatrixX3d m_previousVelocities;
		Eigen::MatrixX3d m_trackerParticles;
		bool m_hasTrackerParticles;
	private:
		void getExtrapolatedVelocityField(double h, Eigen::SparseMatrix<double>& velocityField, Eigen::SparseMatrix<double, Eigen::ColMajor>& velocityMid);
		void getAdvectedVelocityField(double h, Eigen::SparseMatrix<double>& velocityField);
		void advectSignedDistance(double h, const Eigen::SparseMatrix<double, Eigen::ColMajor> oldVelocityMid);
		void advectTrackerParticles(double h, const Eigen::SparseMatrix<double>& velocityField);
		void updateGravity(double h, Eigen::SparseMatrix<double>& velocity, const Eigen::VectorXd& signedDistance);

		void updatePressure(double h, Eigen::SparseMatrix<double>& velocity, Eigen::SparseVector<double>& pressure);
		void updateVelocityFromPressureGradient(double h, Eigen::SparseMatrix<double> velocity);

		void recalculateLevelSet();


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
		bool m_enableGravity;
		bool m_recalculateLevelSet;

		friend class EulerSimulationLevelSetTest;
		friend class EulerSimulationTest_TestConstantAdvection_Test;
		friend class EulerSimulationLevelSetTest_SphereAdvection_Test;
		friend class EulerSimulationLevelSetTest_ConstantVelocityAdvection_Test;
		friend class EulerSimulationLevelSetTest_GravityEnabled_Test;
		friend class EulerSimulationPressureTest;
		friend class EulerSimulationPressureTest_HydrostaticPressure_Test;
		friend class EulerSimulationPressureTest_SteadyStateHydrostaticPressure_Test;

	};

}