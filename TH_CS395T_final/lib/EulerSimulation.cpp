#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Sparse>
//#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include "EulerSimulation.h"
#include "EulerState.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

namespace FluidSimulation
{
	void EulerSimulation::step(double h)
	{
		Eigen::SparseMatrix<double, Eigen::ColMajor> newVelocity, oldVelocityMid;
		Eigen::SparseMatrix<double, Eigen::ColMajor> pressure;
		std::vector<Eigen::Triplet<double>> triplets[3];

		getExtrapolatedVelocityField(h, newVelocity, oldVelocityMid);
		advectSignedDistance(h, oldVelocityMid);
		if (m_recalculateLevelSet)
		{
			recalculateLevelSet();
		}

		if (m_enableGravity)
		{
			// What if the velocity is zero...
			updateGravity(h, newVelocity, m_currentState.m_signedDistance);
		}

		getAdvectedVelocityField(h, newVelocity);

		if (m_enablePressure)
		{
			Eigen::SparseVector<double> pressure;
			updatePressure(h, newVelocity, pressure);
			updateVelocityFromPressureGradient(h, newVelocity);
		}

		//spdlog::info("Setting m_currentState.m_velocity = \n{}", newVelocity);
		m_currentState.m_velocity = newVelocity;
		
	}

	void EulerSimulation::getExtrapolatedVelocityField(double h, Eigen::SparseMatrix<double>& velocityField,
		Eigen::SparseMatrix<double, Eigen::ColMajor>& velocityMid)
	{
		//Eigen::SparseMatrix<double, Eigen::ColMajor> interpolateX, interpolateY, interpolateZ;
		std::vector<Eigen::Triplet<double>> triplets[3], tripletsMid[3];
		Eigen::SparseMatrix<double> interpolateX, interpolateY, interpolateZ;

		interpolateX.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		interpolateY.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		interpolateZ.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		velocityMid.resize(m_currentState.getGridMatrixSize(true), 3);
		velocityField.resize(m_currentState.getGridMatrixSize(true), 3);

		//interpolateX.setIdentity(); interpolateY.setIdentity(); interpolateZ.setIdentity();

		// Extend the velocity field into the air region		
		for (int i = 0; i < m_currentState.getGridMatrixSize(true); i++)
		{
			// Check for grid positions that are within 3 grids of the edge
			//bool withinRange = m_currentState.m_signedDistance(i,0) < 5.0 * m_currentState.m_gridSizeHorizontal.maxCoeff();
			Eigen::Vector3d pos = m_currentState.m_positionsMid.row(i);
			bool withinRange = true;
			if (withinRange && m_currentState.m_signedDistance(i) > 0)
			{
				Eigen::Vector3d closestPosition = m_currentState.m_positionsMid.row(i) - m_currentState.m_signedDistance(i) * m_currentState.m_dSignedDistance.row(i);				
				
				//Eigen::Vector3d closestPosition = closestPoint;
				double closestIndexX, closestIndexY, closestIndexZ;
				double weights[3];
				weights[0] = std::modf(closestPosition(0) / m_currentState.m_gridSizeHorizontal(0) - 0.5, &closestIndexX);
				weights[1] = std::modf(closestPosition(1) / m_currentState.m_gridSizeHorizontal(1) - 0.5, &closestIndexY);
				weights[2] = std::modf(closestPosition(2) / m_currentState.m_gridSizeHorizontal(2) - 0.5, &closestIndexZ);
				size_t closestIndex = (size_t)closestIndexZ + (size_t)closestIndexY * (m_currentState.m_dims(Z) + 1) + (size_t)closestIndexX * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);
				size_t originalClosestIndex = closestIndex;
				Eigen::Vector3d gradient = -m_currentState.m_dSignedDistance.row(i);
				bool definedInitially = true;
				// If we're not inside yet, then step a little further in
				while (closestIndex >= 0 && closestIndex < m_currentState.getGridMatrixSize(true) && 
							 m_currentState.m_signedDistance(closestIndex) > 0 &&
							 !gradient.isZero())
				{
					definedInitially = false;
					int dirToMove, amount;
					amount = abs(gradient.array()).maxCoeff(&dirToMove);

					if (gradient(dirToMove) > 0)
					{
						closestIndex += ((dirToMove == 2) * 1 + (dirToMove == 1) * (m_currentState.m_dims(2) + 1) + (dirToMove == 0) * (m_currentState.m_dims(2) + 1) * (m_currentState.m_dims(1) + 1));
					}
					else
					{
						closestIndex -= ((dirToMove == 2) * 1 + (dirToMove == 1) * (m_currentState.m_dims(2) + 1) + (dirToMove == 0) * (m_currentState.m_dims(2) + 1) * (m_currentState.m_dims(1) + 1));
					}
					weights[dirToMove] = 1 - weights[dirToMove];

					gradient(dirToMove) = 0;

					/*closestPosition -= 0.1 * m_currentState.m_signedDistance(i) * m_currentState.m_dSignedDistance.row(i);
					weights[0] = std::modf(closestPosition(0) / m_currentState.m_gridSizeHorizontal(0), &closestIndexX);
					weights[1] = std::modf(closestPosition(1) / m_currentState.m_gridSizeHorizontal(1), &closestIndexY);
					weights[2] = std::modf(closestPosition(2) / m_currentState.m_gridSizeHorizontal(2), &closestIndexZ);
					closestIndex = (size_t)closestIndexZ + (size_t)closestIndexY * (m_currentState.m_dims(Z) + 1) + (size_t)closestIndexX * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);*/
				}

				for (int dim = 0; dim < 3; dim++)
				{
					//double closestIndex;
					//double weight = std::modf(closestPosition(dim) / m_currentState.m_gridSizeHorizontal(dim), &closestIndex);
					double weight = weights[dim];
					int offset = (dim == Z) * 1 + (dim == Y) * (m_currentState.m_dims(Z) + 1) + (dim == X) * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);
					
					// TODO: Update this. I'm not properly checking if we're within bounds
					int dimIndex = (dim == 2) * (i % (m_currentState.m_dims(2) + 1)) +
						(dim == 1) * ((i / (m_currentState.m_dims(2) + 1)) % (m_currentState.m_dims(1) + 1)) +
						(dim == 0) * (i / ((m_currentState.m_dims(2) + 1) * (m_currentState.m_dims(1) + 1)));

					bool withinBoundaryIncoming = (dimIndex - 1 >= 0);
					bool withinBoundaryOutgoing = (dimIndex + 1 < m_currentState.m_dims(dim));
					
					//bool withinBoundaryOutgoing = (i + offset >= 0) && (i + offset < m_currentState.getGridMatrixSize(true)) && m_currentState.m_velocity.coeff(i + offset, dim);

					//tripletsMid[dim].emplace_back(i, i, (1 - 0.5 * withinBoundaryOutgoing));
					if (withinBoundaryOutgoing)
					{
						//tripletsMid[dim].emplace_back(i, i + offset, 0.5);
					}

					if (closestIndex < 0 || closestIndex >= m_currentState.getGridMatrixSize(true))
					{
						spdlog::critical("Closest within-volume index is outside of the grid domain:");
						spdlog::debug("Position of unknown cell = {}", m_currentState.m_positionsMid.row(i));
						spdlog::debug("Distance to water volume = {} in direction {}", m_currentState.m_signedDistance(i), m_currentState.m_dSignedDistance.row(i).format(LongFmt));
						
						try
						{
							if (originalClosestIndex >= m_currentState.getGridMatrixSize(true))
							{
								spdlog::warn("originalClosestIndex too large ({}), resetting to zero", originalClosestIndex);
								originalClosestIndex = 0;
							}

							if (closestIndex >= m_currentState.getGridMatrixSize(true))
							{
								spdlog::warn("closestIndex too large ({}), resetting to zero", closestIndex);
								closestIndex = 0;
							}

							size_t gridSize = m_currentState.getGridMatrixSize(true) - 1;
							spdlog::debug("Index of unknown cell = {}; Index of original 'closest' cell = {}; Index of 'closest' cell = {}", i, originalClosestIndex, closestIndex);
							spdlog::debug("Distance of unknown cell from surface = {}; distance of original 'closest' cell to surface = {}; distance of 'closest' cell to surface = {}", m_currentState.m_signedDistance(i), m_currentState.m_signedDistance(originalClosestIndex), m_currentState.m_signedDistance(closestIndex));
							spdlog::debug("Gradient of signed distance = {}", m_currentState.m_dSignedDistance.row(i).format(LongFmt));
							/*spdlog::info("Signed distance near original 'closest cell' = \n \t\t\t {} \n {} \t {} \t {} \t | X = \t {} \t {} \t {} \n \t\t\t {}",
													 m_currentState.m_signedDistance(std::min(gridSize, originalClosestIndex + 1)),
													 m_currentState.m_signedDistance(std::max((size_t)0, originalClosestIndex - m_currentState.m_dims(Z) - 1)),
													 m_currentState.m_signedDistance(originalClosestIndex),
													 m_currentState.m_signedDistance(std::min(gridSize, originalClosestIndex + m_currentState.m_dims(Z) + 1)),
													 m_currentState.m_signedDistance(std::max((size_t)0, originalClosestIndex - (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1))),
													 m_currentState.m_signedDistance(originalClosestIndex),
													 m_currentState.m_signedDistance(std::min(gridSize, originalClosestIndex + (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1))),
													 m_currentState.m_signedDistance(std::max((size_t)0, originalClosestIndex - 1)));
							spdlog::info("Signed distance near 'closest cell' = \n \t\t\t {} \n {} \t {} \t {} \t | X = \t {} \t {} \t {} \n \t\t\t {}",
													 m_currentState.m_signedDistance(std::min(gridSize, closestIndex + 1)),
													 m_currentState.m_signedDistance(std::max((size_t)0, closestIndex - m_currentState.m_dims(Z) - 1)),
													 m_currentState.m_signedDistance(closestIndex),
													 m_currentState.m_signedDistance(std::min(gridSize, closestIndex + m_currentState.m_dims(Z) + 1)),
													 m_currentState.m_signedDistance(std::max((size_t)0, closestIndex - (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1))),
													 m_currentState.m_signedDistance(closestIndex),
													 m_currentState.m_signedDistance(std::min(gridSize, closestIndex + (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1))),
													 m_currentState.m_signedDistance(std::max((size_t)0, closestIndex - 1)));*/
							spdlog::debug("Velocity at original 'closest cell' = {}, Velocity at 'closest cell' = {}", m_currentState.m_velocity.row(originalClosestIndex), m_currentState.m_velocity.row(closestIndex));
						}
						catch (const std::exception& e)
						{
							
							spdlog::debug("Caught exception:\n{}", e.what());
						}
						
						continue;
					}

					// If both points are outside water volume
					if (m_currentState.m_signedDistance(closestIndex) > 0 && closestIndex + offset < m_currentState.getGridMatrixSize(true) && m_currentState.m_signedDistance(closestIndex + offset) > 0)
					{
						triplets[dim].emplace_back(i, closestIndex, 1);

						/*spdlog::warn("Want to extrapolate the velocity field but both points are outside of the water volume");
						spdlog::info("Index of unknown cell = {}; Index of original 'closest' cell = {}; Index of 'closest' cell = {}", i, originalClosestIndex, closestIndex);
						spdlog::info("Distance of unknown cell from surface = {}; distance of original 'closest' cell to surface = {}; distance of 'closest' cell to surface = {}", m_currentState.m_signedDistance(i), m_currentState.m_signedDistance(originalClosestIndex), m_currentState.m_signedDistance(closestIndex));
						spdlog::info("Gradient of signed distance = {}", m_currentState.m_dSignedDistance.row(i).format(LongFmt));
						spdlog::info("Signed distance near original 'closest cell' = \n \t\t\t {} \n {} \t {} \t {} \t | X = \t {} \t {} \t {} \n \t\t\t {}",
												 m_currentState.m_signedDistance(originalClosestIndex + 1),
												 m_currentState.m_signedDistance(originalClosestIndex - m_currentState.m_dims(Z) - 1),
												 m_currentState.m_signedDistance(originalClosestIndex),
												 m_currentState.m_signedDistance(originalClosestIndex + m_currentState.m_dims(Z) + 1),
												 m_currentState.m_signedDistance(originalClosestIndex - (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1)),
												 m_currentState.m_signedDistance(originalClosestIndex),
												 m_currentState.m_signedDistance(originalClosestIndex + (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1)),
												 m_currentState.m_signedDistance(originalClosestIndex - 1));
						spdlog::info("Signed distance near 'closest cell' = \n \t\t\t {} \n {} \t {} \t {} \t | X = \t {} \t {} \t {} \n \t\t\t {}",
												 m_currentState.m_signedDistance(closestIndex + 1),
												 m_currentState.m_signedDistance(closestIndex - m_currentState.m_dims(Z) - 1),
												 m_currentState.m_signedDistance(closestIndex),
												 m_currentState.m_signedDistance(closestIndex + m_currentState.m_dims(Z) + 1),
												 m_currentState.m_signedDistance(closestIndex - (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1)),
												 m_currentState.m_signedDistance(closestIndex),
												 m_currentState.m_signedDistance(closestIndex + (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1)),
												 m_currentState.m_signedDistance(closestIndex - 1));
						spdlog::info("Velocity at original 'closest cell' = {}, Velocity at 'closest cell' = {}", m_currentState.m_velocity.row(originalClosestIndex), m_currentState.m_velocity.row(closestIndex));*/
					}
					// If just one point is outside water volume, take the value of the point that is inside
					if (closestIndex + offset < m_currentState.getGridMatrixSize(true) && m_currentState.m_signedDistance(closestIndex) > 0)
					{
						triplets[dim].emplace_back(i, closestIndex + offset, 1);
					}
					else if (closestIndex < m_currentState.getGridMatrixSize(true) && m_currentState.m_signedDistance(closestIndex) > 0 || !definedInitially)
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
					/*int offset = (dim == Z) * 1 + (dim == Y) * (m_currentState.m_dims(Z) + 1) + (dim == X) * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);
					bool withinBoundaryOutgoing = (i + offset >= 0) && (i + offset < m_currentState.getGridMatrixSize(true));

					tripletsMid[dim].emplace_back(i, i, (1 - 0.5 * withinBoundaryOutgoing));
					if (withinBoundaryOutgoing)
					{
						tripletsMid[dim].emplace_back(i, i + offset, 0.5);
					}*/
				}
			}

			for (int dim = 0; dim < 3; dim++)
			{
				int offset = (dim == Z) * 1 + (dim == Y) * (m_currentState.m_dims(Z) + 1) + (dim == X) * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);
				int dimIndex = (dim == 2) * (i % (m_currentState.m_dims(2) + 1)) +
					(dim == 1) * ((i / (m_currentState.m_dims(2) + 1)) % (m_currentState.m_dims(1) + 1)) +
					(dim == 0) * (i / ((m_currentState.m_dims(2) + 1) * (m_currentState.m_dims(1) + 1)));

				bool withinBoundaryIncoming = (dimIndex - 1 >= 0);
				bool withinBoundaryOutgoing = (dimIndex + 1 < m_currentState.m_dims(dim));

				
				//bool withinBoundaryOutgoing = (i + offset >= 0) && (i + offset < m_currentState.getGridMatrixSize(true));

				tripletsMid[dim].emplace_back(i, i, (1 - 0.5 * withinBoundaryOutgoing));
				if (withinBoundaryOutgoing)
				{
					tripletsMid[dim].emplace_back(i, i + offset, 0.5);
				}
			}
		}



		interpolateX.setZero(); interpolateY.setZero(); interpolateZ.setZero();
		interpolateX.setFromTriplets(triplets[0].begin(), triplets[0].end());
		interpolateY.setFromTriplets(triplets[1].begin(), triplets[1].end());
		interpolateZ.setFromTriplets(triplets[2].begin(), triplets[2].end());

		spdlog::debug("Setting Velocity from Extrapolation");
		Eigen::SparseMatrix<double, Eigen::ColMajor> temp;
		temp = interpolateX * m_currentState.m_velocity.col(0);
		velocityField.col(0) = temp;
		temp = interpolateY * m_currentState.m_velocity.col(1);
		velocityField.col(1) = temp;
		temp = interpolateZ * m_currentState.m_velocity.col(2);
		velocityField.col(2) = temp;
		//velocityMid.setFromTriplets(tripletsMid.begin(), tripletsMid.end(), [](const double& a, const double& b) { return a + b; });

		interpolateX.setZero(); interpolateY.setZero(); interpolateZ.setZero();
		interpolateX.setFromTriplets(tripletsMid[0].begin(), tripletsMid[0].end());
		interpolateY.setFromTriplets(tripletsMid[1].begin(), tripletsMid[1].end());
		interpolateZ.setFromTriplets(tripletsMid[2].begin(), tripletsMid[2].end());

		interpolateX.makeCompressed();
		interpolateY.makeCompressed();
		interpolateZ.makeCompressed();
		velocityMid.col(0) = interpolateX * velocityField.col(0);
		velocityMid.col(1) = interpolateY * velocityField.col(1);
		velocityMid.col(2) = interpolateZ * velocityField.col(2);


		/*spdlog::debug("InterpolateX = \n{}", interpolateX);
		spdlog::debug("InterpolateY = \n{}", interpolateY);
		spdlog::debug("InterpolateZ = \n{}", interpolateZ);*/
		spdlog::debug("Extrapolated Velocity Field = \n{}", velocityField);
	}

	void EulerSimulation::getAdvectedVelocityField(double h, Eigen::SparseMatrix<double>& velocityField)
	{
		Eigen::SparseMatrix<double> interpolateX, interpolateY, interpolateZ;
		std::vector<Eigen::Triplet<double>> triplets[3];

		interpolateX.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		interpolateY.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));
		interpolateZ.resize(m_currentState.getGridMatrixSize(true), m_currentState.getGridMatrixSize(true));

		triplets[0].clear();
		triplets[1].clear();
		triplets[2].clear();

		// OK to skip if velocity is 0
		for (int k = 0; k < velocityField.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(velocityField, k); it; ++it)
			{
				int row = it.row();
				double value = it.value() * h / m_currentState.m_gridSizeHorizontal(it.col());
				int direction = (value < 0) ? 1 : -1;
				int offset = (it.col() == Z) * 1 + (it.col() == Y) * (m_currentState.m_dims(Z) + 1) + (it.col() == X) * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);
				int dimIndex = (it.col() == 2) * (it.row() % (m_currentState.m_dims(2) + 1)) +
					(it.col() == 1) * ((it.row() / (m_currentState.m_dims(2) + 1)) % (m_currentState.m_dims(1) + 1)) +
					(it.col() == 0) * (it.row() / ((m_currentState.m_dims(2) + 1) * (m_currentState.m_dims(1) + 1)));

				bool withinBoundaryIncoming = (dimIndex - 1 >= 0);
				bool withinBoundaryOutgoing = (dimIndex + 1 < m_currentState.m_dims(it.col()));
				bool within2BoundaryOutgoing = (dimIndex < m_currentState.m_dims(it.col()));

				/*bool withinBoundaryIncoming = (it.row() - offset >= 0) && (it.row() - offset < m_currentState.getGridMatrixSize(true));
				bool withinBoundaryOutgoing = (it.row() + offset >= 0) && (it.row() + offset < m_currentState.getGridMatrixSize(true));
				bool within2BoundaryOutgoing = (it.row() + 2 * offset >= 0) && (it.row() + 2 * offset < m_currentState.getGridMatrixSize(true));*/

				// If the cell is at an edge, keep velocity at zero per Neumann condition
				if (withinBoundaryIncoming && (withinBoundaryOutgoing || it.col() == Z))
				{
					spdlog::debug("Inserting Interpolation Value {} at ({}, {}) into dimension {}", 1 - std::abs(value), it.row(), it.row(), it.col());
					triplets[it.col()].emplace_back(it.row(), it.row(), 1 - std::abs(value));

					if (value != 0.0 && it.row() + direction * offset >= 0 && it.row() + direction * offset < m_currentState.getGridMatrixSize(true))
					{
						spdlog::debug("Inserting Interpolation Value {} at ({}, {}) into dimension {}", std::abs(value), it.row(), it.row() + direction * offset, it.col());
						triplets[it.col()].emplace_back(it.row(), it.row() + direction * offset, std::abs(value));
					}
					else
					{
						spdlog::debug("Skipping advection interpolation value at ({}, {})", it.row(), it.row() + direction * offset);
					}
				}
			}
		}
		interpolateX.setZero(); interpolateY.setZero(); interpolateZ.setZero();

		interpolateX.setFromTriplets(triplets[0].begin(), triplets[0].end());
		interpolateY.setFromTriplets(triplets[1].begin(), triplets[1].end());
		interpolateZ.setFromTriplets(triplets[2].begin(), triplets[2].end());

		//spdlog::debug("Interpolate size = {} x {}", interpolateX.rows(), interpolateX.cols());
		Eigen::SparseMatrix<double, Eigen::ColMajor> temp;
		velocityField.makeCompressed();
		if (interpolateX.nonZeros() > 0)
		{
			interpolateX.makeCompressed();

			temp = (interpolateX * velocityField);
			temp = temp.pruned();
			//spdlog::debug("Temp size = {} x {}", temp.rows(), temp.cols());
			velocityField.col(0) = temp.pruned().col(0);
		}
		if (interpolateY.nonZeros() > 0)
		{
			interpolateY.makeCompressed();

			temp = (interpolateY * velocityField);
			velocityField.col(1) = temp.pruned().col(1);
		}
		if (interpolateZ.nonZeros() > 0)
		{
			interpolateZ.makeCompressed();
			temp = interpolateZ * velocityField;
			velocityField.col(2) = temp.pruned().col(2);
		}
		velocityField.makeCompressed();
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

				//spdlog::info("Row = {}; Col = {}; Offset = {}; Block Size = {}", it.row(), it.col(), offset, blockSize);

				if (it.row() + blockSize < m_currentState.getGridMatrixSize(true))
				{
					if (value != 0.0 &&
						it.row() >= 0 &&
						it.row() + direction * offset >= 0 &&
						it.row() + direction * offset + blockSize >= 0 &&
						it.row() + blockSize >= 0 &&
						it.row() + direction * offset < m_currentState.getGridMatrixSize(true) &&
						it.row() + direction * offset + blockSize < m_currentState.getGridMatrixSize(true) &&
						it.row() + blockSize < m_currentState.getGridMatrixSize(true))
					{
						interpolateMatrix.block(it.row(), it.row() + direction * offset, 1, blockSize) = (interpolateMatrix.block(it.row(), it.row(), 1, blockSize) * std::abs(value)).eval();
					}


					interpolateMatrix.block(it.row(), it.row(), 1, blockSize) = (interpolateMatrix.block(it.row(), it.row(), 1, blockSize) * (1.0 - std::abs(value))).eval();
					//(it.row(), it.row()) -= std::abs(value);
				}
			}
		}

		//spdlog::info("Signed Distance Advection: Interpolate Matrix = \n{}", interpolateMatrix);
		/*spdlog::debug("Old Velocity Mid ({} -> {}, {} -> {}, {} -> {}, {} -> {})", 
			3, oldVelocityMid.coeff(3 + 5 * 3 + 5 * 3 * 3, 2),
			2, oldVelocityMid.coeff(2 + 5 * 3 + 5 * 3 * 3, 2),
			1, oldVelocityMid.coeff(1 + 5 * 3 + 5 * 3 * 3, 2),
			0, oldVelocityMid.coeff(0 + 5 * 3 + 5 * 3 * 3, 2));*/
		//spdlog::info("Interpolation Matrix Block = \n{}", interpolateMatrix.block<5, 5>(5 * 3 + 5 * 3 * 3, 5 * 3 + 5 * 3 * 3));
		//spdlog::info("Signed Distance Block = \n{}", m_currentState.m_signedDistance.segment<5>(5 * 3 + 5 * 3 * 3));
		spdlog::debug("Prior Signed Distance = \n{}", m_currentState.m_signedDistance);
		m_currentState.m_signedDistance = (interpolateMatrix * m_currentState.m_signedDistance).eval();
		spdlog::debug("Updated Signed Distance = \n{}", m_currentState.m_signedDistance);
		m_currentState.m_dSignedDistance += h * oldVelocityMid;
		m_currentState.m_dSignedDistance.rowwise().normalize();
		//oldVelocityMid.cwiseQuotient(m_currentState.m_gridSizeHorizontal) * h;
	}

	void EulerSimulation::recalculateLevelSet()
	{

	}

	void EulerSimulation::updateGravity(double h, Eigen::SparseMatrix<double>& velocity, const Eigen::VectorXd& signedDistance)
	{
		Eigen::SparseMatrix<double>::InnerIterator it(velocity, Z);
		std::vector<int> valuesToInsert;

		for (int x = 0; x < m_currentState.m_dims(X) + 1; x++)
		{
			for (int y = 0; y < m_currentState.m_dims(Y) + 1; y++)
			{
				for (int z = 0; z < m_currentState.m_dims(Z) + 1; z++)
				{
					int i = z + y * (m_currentState.m_dims(Z) + 1) + x * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);

					if (it && it.row() == i)
					{
						it.valueRef() += h * m_gravity;
						++it;
					} 
					// TODO: OVERRIDE AND WRITE THE VELOCITY EVERYWHERE
					else if (true || signedDistance(i) <= 3 * m_currentState.m_gridSizeHorizontal.maxCoeff())
					{
							valuesToInsert.push_back(i);
					}
				}
			}
		}
		for (std::vector<int>::const_iterator it = valuesToInsert.begin(); it != valuesToInsert.end(); it++)
		{
			velocity.insert(*it, Z) = h * m_gravity;
		}
		velocity.makeCompressed();
	}

	void EulerSimulation::updatePressure(double h, Eigen::SparseMatrix<double>& velocity, Eigen::SparseVector<double>& pressure)
	{
		// Divergence of velocity
		Eigen::VectorXd p = m_currentState.m_pressure;
		Eigen::SparseMatrix<double> dv, d2p, dp;

		m_currentState.getLaplacianOperator(d2p);

		for (int x = 0; x < m_currentState.m_dims(X); x++)
		{
			for (int y = 0; y < m_currentState.m_dims(Y); y++)
			{
				for (int z = 0; z < m_currentState.m_dims(Z); z++)
				{
					int iElement = z + y * (m_currentState.m_dims(Z) + 1) + x * (m_currentState.m_dims(Z) + 1) * (m_currentState.m_dims(Y) + 1);
					int i = z + y * (m_currentState.m_dims(Z)) + x * (m_currentState.m_dims(Z)) * (m_currentState.m_dims(Y));

					if (m_currentState.m_signedDistance(iElement) > 0)
					{
						d2p.coeffRef(i, i) = 1.0 / std::pow(m_currentState.m_gridSizeHorizontal(X), 2) + 
							1.0 / std::pow(m_currentState.m_gridSizeHorizontal(Y), 2) +
							1.0 / std::pow(m_currentState.m_gridSizeHorizontal(Z), 2);
						// Make sure this air cell doesn't factor into any other cells
						// Also make sure this cell is stand-alone so it receives the value of the divergence at that point
						if (i - 1 >= 0)
						{
							d2p.coeffRef(i - 1, i) = 0;
							d2p.coeffRef(i, i - 1) = 0;
						}
						if (i - m_currentState.m_dims(Z) >= 0)
						{
							d2p.coeffRef(i - m_currentState.m_dims(Z), i) = 0;
							d2p.coeffRef(i, i - m_currentState.m_dims(Z)) = 0;
						}
						if (i - m_currentState.m_dims(Z) * m_currentState.m_dims(Y) >= 0)
						{
							d2p.coeffRef(i - m_currentState.m_dims(Z) * m_currentState.m_dims(Y), i) = 0;
							d2p.coeffRef(i, i - m_currentState.m_dims(Z) * m_currentState.m_dims(Y)) = 0;
						}
						if (i + 1 < m_currentState.getGridMatrixSize(false))
						{
							d2p.coeffRef(i + 1, i) = 0;
							d2p.coeffRef(i, i + 1) = 0;
						}
						if (i + m_currentState.m_dims(Z) < m_currentState.getGridMatrixSize(false))
						{
							d2p.coeffRef(i + m_currentState.m_dims(Z), i) = 0;
							d2p.coeffRef(i, i + m_currentState.m_dims(Z)) = 0;
						}
						if (i + m_currentState.m_dims(Z) * m_currentState.m_dims(Y) < m_currentState.getGridMatrixSize(false))
						{
							d2p.coeffRef(i + m_currentState.m_dims(Z) * m_currentState.m_dims(Y), i) = 0;
							d2p.coeffRef(i, i + m_currentState.m_dims(Z) * m_currentState.m_dims(Y)) = 0;
						}


					}
				}
			}
		}


		d2p *= h / m_density;

		//spdlog::debug("Solving for pressure with the following inputs");
		//spdlog::debug("velocity = \n{}", velocity);
		spdlog::debug("d2p = \n{}", (Eigen::MatrixXd)d2p);
		spdlog::debug("d2p (block) = \n{}", ((Eigen::MatrixXd)d2p).block(0 + 0 + (m_currentState.m_dims(X) / 2) * m_currentState.m_dims(Z) * m_currentState.m_dims(Y), 0, 
																																		 4, m_currentState.getGridMatrixSize(false)));

		/* SOLVE */
		// Solve with modified Cholesky preconditioner
		//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cgSolver;

		//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double, Eigen::Lower | Eigen::Upper>> cgSolver;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(d2p);
		
		spdlog::debug("velocity = \n{}", (Eigen::MatrixXd)velocity);
		m_currentState.getQuantityDivergence(dv, velocity);

		for (int k = 0; k < dv.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(dv, k); it; ++it)
			{
				int z = (it.row() % (m_currentState.m_dims(2)));
				int y = ((it.row() / (m_currentState.m_dims(2))) % (m_currentState.m_dims(1)));
				int x =	(it.row() / ((m_currentState.m_dims(2)) * (m_currentState.m_dims(1))));

				if (m_currentState.m_signedDistance(z + y * (m_currentState.m_dims(2) + 1) + x * (m_currentState.m_dims(2) + 1) * (m_currentState.m_dims(1) + 1)) > 0)
				{
					it.valueRef() = 0;
				}
			}
		}

		spdlog::debug("div v = \n{}", (Eigen::MatrixXd)dv);
		//solver.compute(d2p);
		pressure = solver.solve(-dv);
		//pressure = p.sparseView();
		m_currentState.m_pressure = pressure;
		spdlog::debug("pressure = \n{}", (Eigen::MatrixXd)pressure);

		m_currentState.getPressureGradient(dp);
		spdlog::debug("dp = \n{}", (Eigen::MatrixXd)dp);
		spdlog::debug("change in velocity = \n{}", h / m_density * m_currentState.m_midToElement.transpose() * dp);
		// Note: The way I have defined velocity, 
		velocity -= h / m_density * m_currentState.m_midToElement.transpose() * dp;
		//p = cgSolver.solveWithGuess(-(Eigen::VectorXd)dv, p);

		spdlog::debug("Solver info: {} in {} iterations", solver.info() == Eigen::Success ? "Successfully converged" : "Did not successfully converge", "n/a");
		//spdlog::info("Solver error: {}", solver)

		spdlog::debug("Results");
		//spdlog::info("p = \n{}", pressure);


		//pressure = p.sparseView();
		//m_currentState.m_pressure = pressure;
	}

	void EulerSimulation::updateVelocityFromPressureGradient(double h, Eigen::SparseMatrix<double> velocity)
	{
		Eigen::SparseMatrix<double> dp;

		//m_currentState.getPressureGradient(dp);

		//velocity -= h / m_density * m_currentState.m_midToElement.transpose() * dp;
	}


}