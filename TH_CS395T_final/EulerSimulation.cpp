#include <vector>

#include <Eigen/Base>
#include <Eigen/Sparse>
#include "EulerSimulation.h"
#include "EulerState.h"

namespace FluidSimulation
{
	void EulerSimulation::step(double h)
	{
		Eigen::SparseMatrix<double> v, p, f, dv;
		std::vector<Eigen::Triplet<double>> velocityTriplets, pressureTriplets;

		f = m_currentState.getForcesGrid();
		v = m_currentState.getVelocityGrid();

		/* SOLVE VELOCITY */
		velocityTriplets.reserve(v.nonZeros());

		dv = v;

		v + h * p;

		forces = p / m_density - m_gravity * ;
		
	}
}