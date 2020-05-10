#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "EulerSimulation.h"
#include "EulerState.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

namespace FluidSimulation
{
	void EulerSimulation::step(double h)
	{
		Eigen::SparseMatrix<double, Eigen::ColMajor> interpolateX, interpolateY, interpolateZ;
		Eigen::SparseMatrix<double, Eigen::ColMajor> newVelocity, oldVelocityMid;
		std::vector<Eigen::Triplet<double>> triplets[3];






		// Interpolate a 
		// m_velocity = 0  -> w(i) = 1;
		// m_velocity = 1  -> w(i) = 0; w(i-1) = 1;
		// m_velocity = -1 -> w(i) = 0; w(i+1) = 1;
		/*interpolateX.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		interpolateY.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		interpolateZ.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		newVelocity.resize(m_currentState.getGridMatrixSize(true), 3);

		spdlog::debug("m_velocity size = {} x {}", m_currentState.m_velocity.rows(), m_currentState.m_velocity.cols());
		//interpolateX.setIdentity();
		//interpolateX.diagonal() -= -m_currentState.m_velocity.col(X) * h / m_currentState.m_gridSizeHorizontal(X);
		// TODO: This is Forward Euler
		for (int k = 0; k < m_currentState.m_velocity.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(m_currentState.m_velocity, k); it; ++it)
			{
				int direction;
				double value = it.value() * h / m_currentState.m_gridSizeHorizontal(it.col());
				spdlog::debug("Inserting Interpolation Value {} at ({}, {}) into dimension {}", 1 - std::abs(value), it.row(), it.row(), it.col());
				triplets[it.col()].emplace_back(it.row(), it.row(), 1 - std::abs(value));

				if (value != 0.0)
				{
					direction = (value < 0) ? 1 : -1;
					spdlog::debug("Inserting Interpolation Value {} at ({}, {}) into dimension {}", std::abs(value), it.row(), it.row() + direction, it.col());
					triplets[it.col()].emplace_back(it.row(), it.row() + direction, std::abs(value));
				}
			}
		}
		interpolateX.setFromTriplets(triplets[0].begin(), triplets[0].end());
		interpolateY.setFromTriplets(triplets[1].begin(), triplets[1].end());
		interpolateZ.setFromTriplets(triplets[2].begin(), triplets[2].end());

		// Backtrace velocity
		spdlog::debug("Interpolate size = {} x {}", interpolateX.rows(), interpolateX.cols());
		Eigen::SparseMatrix<double, Eigen::ColMajor> temp;
		temp = (interpolateX * m_currentState.m_velocity);
		spdlog::debug("Temp size = {} x {}", temp.rows(), temp.cols());
		newVelocity.col(0) = temp.col(0);
		temp = (interpolateY * m_currentState.m_velocity);
		newVelocity.col(1) = temp.col(1);
		temp = interpolateZ * m_currentState.m_velocity;
		newVelocity.col(2) = temp.col(2);*/

		getAdvectedVelocityField(h, newVelocity, interpolateX, interpolateY, interpolateZ, oldVelocityMid);
		advectSignedDistance(h, oldVelocityMid);


		spdlog::debug("Setting m_currentState.m_velocity");
		m_currentState.m_velocity = newVelocity;

		/*m_currentState.m_signedDistance = interpolateX * interpolateY * interpolateZ * m_currentState.m_signedDistance.col(0);
		m_currentState.m_signedDistance.col(1) = interpolateX * m_currentState.m_signedDistance.col(1);
		m_currentState.m_signedDistance.col(2) = interpolateY * m_currentState.m_signedDistance.col(2);
		m_currentState.m_signedDistance.col(3) = interpolateZ * m_currentState.m_signedDistance.col(3);*/

		// Calculate the position of cells to place the new values at
		//prevPosMid = m_currentState.m_positionsMid - m_currentState.m_velocity;

		// "Interpolate" prevPos into a matrix of forcing values aligned to known values grid index locations
		// Gives a matrix of gridPosition.offset
		//prevPosMid.cwiseQuotient(m_currentState.m_gridSizeHorizontal);
		// (n x n) * known values (n x 1) -> new values (n x 1)

		// Calculate the forces of the current cells



		//std::vector<Eigen::Triplet<double>> velocityTriplets, pressureTriplets;

		//m_currentState.getForcesGrid(f);
		//m_currentState.getVelocityGrid(v);

		/* SOLVE VELOCITY */
		//velocityTriplets.reserve(v.nonZeros());

		//dv = v;

		//v + h * p;

		//forces = p / m_density - m_gravity * ;
		
	}

	void EulerSimulation::getAdvectedVelocityField(double h, Eigen::SparseMatrix<double>& velocityField, 
		Eigen::SparseMatrix<double>& interpolateX, Eigen::SparseMatrix<double>& interpolateY, Eigen::SparseMatrix<double>& interpolateZ, 
		Eigen::SparseMatrix<double, Eigen::ColMajor>& oldVelocityMid)
	{
		//Eigen::SparseMatrix<double, Eigen::ColMajor> interpolateX, interpolateY, interpolateZ;
		std::vector<Eigen::Triplet<double>> triplets[3], tripletsMid;

		interpolateX.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		interpolateY.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		interpolateZ.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		oldVelocityMid.resize(m_currentState.getGridMatrixSize(true), 3);
		velocityField.resize(m_currentState.getGridMatrixSize(true), 3);

		//interpolateX.setIdentity(); interpolateY.setIdentity(); interpolateZ.setIdentity();
		
		// Extend the velocity field into the air region		
		for (int i = 0; i < m_currentState.getGridMatrixSize(true); i++)
		{
			// Check for grid positions that are within 3 grids of the edge
			//bool withinRange = m_currentState.m_signedDistance(i,0) < 5.0 * m_currentState.m_gridSizeHorizontal.maxCoeff();
			bool withinRange = true;
			if (withinRange && m_currentState.m_signedDistance(i) > 0)
			{
				Eigen::Vector3d closestPosition = m_currentState.m_positionsMid.row(i) - m_currentState.m_signedDistance(i) * m_currentState.m_dSignedDistance.row(i);
				
				double closestIndexX, closestIndexY, closestIndexZ;
				double weights[3];
				weights[0] = std::modf(closestPosition(0) / m_currentState.m_gridSizeHorizontal(0), &closestIndexX);
				weights[1] = std::modf(closestPosition(1) / m_currentState.m_gridSizeHorizontal(1), &closestIndexY);
				weights[2] = std::modf(closestPosition(2) / m_currentState.m_gridSizeHorizontal(2), &closestIndexZ);
				size_t closestIndex = (size_t)closestIndexZ + (size_t)closestIndexY * (m_currentState.m_dims(Z) + 1) + (size_t)closestIndexX * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);

				for (int dim = 0; dim < 3; dim++)
				{
					//double closestIndex;
					//double weight = std::modf(closestPosition(dim) / m_currentState.m_gridSizeHorizontal(dim), &closestIndex);
					double weight = weights[dim];
					int offset = (dim == Z) * 1 + (dim == Y) * (m_currentState.m_dims(Z) + 1) + (dim == X) * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);

					if (closestIndex < 0 || closestIndex >= m_currentState.getGridMatrixSize(true))
					{
						spdlog::critical("Closest within-volume index is outside of the grid domain:");
						spdlog::debug("Position of unknown cell = {}", m_currentState.m_positionsMid.row(i));
						spdlog::debug("Distance to water volume = {} in direction {}", m_currentState.m_signedDistance(i), m_currentState.m_dSignedDistance.row(i).format(LongFmt));
						continue;
					}

					// If both points are outside water volume
					if (m_currentState.m_signedDistance(closestIndex) > 0 && closestIndex + offset < m_currentState.getGridMatrixSize(true) && m_currentState.m_signedDistance(closestIndex + offset) > 0)
					{
						spdlog::warn("Want to extrapolate the velocity field but both points are outside of the water volume");
					}
					// If just one point is outside water volume, take the value of the point that is inside
					if (m_currentState.m_signedDistance(closestIndex) > 0)
					{
						triplets[dim].emplace_back(i, closestIndex + offset, 1);
					}
					else if (closestIndex + offset < m_currentState.getGridMatrixSize(true) && m_currentState.m_signedDistance(closestIndex + offset) > 0)
					{
						triplets[dim].emplace_back(i, closestIndex, 1);
					}
					// Otherwise interpolate the two points
					else
					{
						//spdlog::debug("Inserting (air) Interpolation Value {} at ({}, {}) into dimension {}", 1 - weight, i, closestIndex, dim);
						triplets[dim].emplace_back(i, closestIndex, 1 - weight);
						if (weight > 0)
						{
							if (closestIndex + offset < 0 || closestIndex + offset >= m_currentState.getGridMatrixSize(true))
							{
								spdlog::debug("Boundary Edge");
								continue;
							}
							//spdlog::debug("Inserting (air) Interpolation Value {} at ({}, {}) into dimension {}", weight, i, closestIndex + offset, dim);
							triplets[dim].emplace_back(i, closestIndex + offset, weight);
						}
					}
				}
			}
			else if (m_currentState.m_signedDistance(i) <= 0)
			{
				for (int dim = 0; dim < 3; dim++)
				{
					triplets[dim].emplace_back(i, i, 1);
				}
			}
		}

		interpolateX.setFromTriplets(triplets[0].begin(), triplets[0].end());
		interpolateY.setFromTriplets(triplets[1].begin(), triplets[1].end());
		interpolateZ.setFromTriplets(triplets[2].begin(), triplets[2].end());

		spdlog::debug("Setting Velocity from Extrapolation");
		Eigen::SparseMatrix<double, Eigen::ColMajor> temp;
		temp = interpolateX * m_currentState.m_velocity;
		velocityField.col(0) = temp.col(0);
		temp = interpolateY * m_currentState.m_velocity;
		velocityField.col(1) = temp.col(1);
		temp = interpolateZ * m_currentState.m_velocity;
		velocityField.col(2) = temp.col(2);

		//spdlog::debug("InterpolateX = \n{}", interpolateX);
		//spdlog::debug("InterpolateY = \n{}", interpolateY);
		//spdlog::debug("InterpolateZ = \n{}", interpolateZ);
		spdlog::debug("Extrapolated Velocity Field = \n{}", velocityField);
		
		triplets[0].clear();
		triplets[1].clear();
		triplets[2].clear();

		for (int k = 0; k < velocityField.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(velocityField, k); it; ++it)
			{
				double value = it.value() * h / m_currentState.m_gridSizeHorizontal(it.col());
				int direction = (value < 0) ? 1 : -1;
				int offset = (it.col() == Z) * 1 + (it.col() == Y) * (m_currentState.m_dims(Z) + 1) + (it.col() == X) * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);
				bool withinBoundaryOutgoing = (it.row() - offset >= 0) && (it.row() - offset < m_currentState.getGridMatrixSize(true));
				bool withinBoundaryIncoming = (it.row() + offset >= 0) && (it.row() + offset < m_currentState.getGridMatrixSize(true));

				spdlog::debug("Inserting Interpolation Value {} at ({}, {}) into dimension {}", 1 - std::abs(value), it.row(), it.row(), it.col());
				triplets[it.col()].emplace_back(it.row(), it.row(), 1 - std::abs(value));
				tripletsMid.emplace_back(it.row(), it.col(), (1 - 0.5 * withinBoundaryIncoming) * it.value());

				if (value != 0.0 && it.row() + direction * offset >= 0 && it.row() + direction * offset < m_currentState.getGridMatrixSize(true))
				{
					spdlog::debug("Inserting Interpolation Value {} at ({}, {}) into dimension {}", std::abs(value), it.row(), it.row() + direction * offset, it.col());
					triplets[it.col()].emplace_back(it.row(), it.row() + direction * offset, std::abs(value));
				}

				/* CALCULATE MID INTERPOLATION */
				if (withinBoundaryOutgoing)
				{
					tripletsMid.emplace_back(it.row() - offset, it.col(), 0.5 * it.value());
				}
			}
		}
		interpolateX.setZero(); interpolateY.setZero(); interpolateZ.setZero();
		
		interpolateX.setFromTriplets(triplets[0].begin(), triplets[0].end());
		interpolateY.setFromTriplets(triplets[1].begin(), triplets[1].end());
		interpolateZ.setFromTriplets(triplets[2].begin(), triplets[2].end());
		oldVelocityMid.setFromTriplets(tripletsMid.begin(), tripletsMid.end(), [](const double& a, const double& b) { return a + b; });

		//spdlog::debug("Interpolate size = {} x {}", interpolateX.rows(), interpolateX.cols());
		//Eigen::SparseMatrix<double, Eigen::ColMajor> temp;
		//temp = (interpolateX * m_currentState.m_velocity);
		//spdlog::debug("Temp size = {} x {}", temp.rows(), temp.cols());
		velocityField.col(0) = (interpolateX * velocityField).col(0);
		//temp = (interpolateY * m_currentState.m_velocity);
		velocityField.col(1) = (interpolateY * velocityField).col(1);
		//temp = interpolateZ * m_currentState.m_velocity;
		velocityField.col(2) = (interpolateZ * velocityField).col(2);
	}

	void EulerSimulation::advectSignedDistance(double h, const Eigen::SparseMatrix<double, Eigen::ColMajor> oldVelocityMid)
	{
		Eigen::SparseMatrix<double> newPositions(m_currentState.getGridMatrixSize(true), 3);
		Eigen::MatrixXd interpolateMatrix(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		interpolateMatrix.setIdentity();

		std::vector<Eigen::Triplet<double>> triplets;
		
		spdlog::debug("Velocity Midpoints = \n{}", oldVelocityMid);

		// Columns
		for (int k = 0; k < oldVelocityMid.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(oldVelocityMid, k); it; ++it)
			{
				double value = it.value() * h / m_currentState.m_gridSizeHorizontal(it.col());
				int direction = (value < 0) ? 1 : -1;
				int offset = (it.col() == Z) * 1 + (it.col() == Y) * (m_currentState.m_dims(Z) + 1) + (it.col() == X) * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);
				int blockSize = (it.col() == Z) * 1 + (it.col() == Y) * 2 + (it.col() == X) * (m_currentState.m_dims(Z) + 1);

				if (value != 0.0 && it.row() + direction * offset >= 0 && it.row() + direction * offset < m_currentState.getGridMatrixSize(true))
				{
					interpolateMatrix.block(it.row(), it.row() + direction * offset, 1, blockSize) = interpolateMatrix.block(it.row(), it.row(), 1, blockSize) * std::abs(value);
				}
				
				interpolateMatrix.block(it.row(), it.row(), 1, blockSize) = interpolateMatrix.block(it.row(), it.row(), 1, blockSize) * (1.0 - std::abs(value));
				//(it.row(), it.row()) -= std::abs(value);
			}
		}

		//spdlog::debug("Signed Distance Advection: Interpolate Matrix = \n{}", interpolateMatrix);
		/*spdlog::debug("Old Velocity Mid ({} -> {}, {} -> {}, {} -> {}, {} -> {})", 
			3, oldVelocityMid.coeff(3 + 5 * 3 + 5 * 3 * 3, 2),
			2, oldVelocityMid.coeff(2 + 5 * 3 + 5 * 3 * 3, 2),
			1, oldVelocityMid.coeff(1 + 5 * 3 + 5 * 3 * 3, 2),
			0, oldVelocityMid.coeff(0 + 5 * 3 + 5 * 3 * 3, 2));*/
		/*spdlog::debug("Interpolation Matrix Block = \n{}", interpolateMatrix.block<5, 5>(5 * 3 + 5 * 3 * 3, 5 * 3 + 5 * 3 * 3));
		spdlog::debug("Signed Distance Block = \n{}", m_currentState.m_signedDistance.segment<5>(5 * 3 + 5 * 3 * 3));*/
		m_currentState.m_signedDistance = (interpolateMatrix * m_currentState.m_signedDistance).eval();
		//oldVelocityMid.cwiseQuotient(m_currentState.m_gridSizeHorizontal) * h;
	}


}