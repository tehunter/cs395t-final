#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "EulerSimulation.h"
#include "EulerState.h"

namespace FluidSimulation
{
	void EulerSimulation::step(double h)
	{
		Eigen::SparseMatrix<double> v, p, f, dv;
		std::vector<Eigen::Triplet<double>> velocityTriplets, pressureTriplets;

		m_currentState.getForcesGrid(f);
		m_currentState.getVelocityGrid(v);

		/* SOLVE VELOCITY */
		velocityTriplets.reserve(v.nonZeros());

		dv = v;

		v + h * p;

		//forces = p / m_density - m_gravity * ;
		
	}
}