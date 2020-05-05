#pragma once

#include <Eigen/Sparse>
#include "EulerState.h"

namespace FluidSimulation
{
	class EulerSimulation
	{
	public:
		int step(double h);
	private:
		// TODO: Calculate timestep from formula
		double getTimeStep() { return 0.001; };

		// State of current fluid simulation
		EulerState m_currentState;

		// Position of fluid box in world coordinates
		Eigen::Vector3d m_position;
		// Axis-angle rotation of fluid box
		Eigen::Vector3d m_theta;

		// Assume incompressible fluid and constant temperature/concentration - density is constant
		double m_density = 1.0;
		double m_gravity = -9.8;
	};

}