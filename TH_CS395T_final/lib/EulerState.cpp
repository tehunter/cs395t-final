#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "EulerState.h"


namespace FluidSimulation
{
	EulerState::EulerState(Eigen::Vector3i dims, Eigen::Vector3d gridSize) : m_dims(dims), m_gridSizeHorizontal(gridSize) 
	{ 
		calculateCentralDifferenceStencil(X, false, m_stencilX);
		calculateCentralDifferenceStencil(Y, false, m_stencilY);
		calculateCentralDifferenceStencil(Z, false, m_stencilZ);

		calculateCentralDifferenceStencil(X, true, m_stencilXMid);
		calculateCentralDifferenceStencil(Y, true, m_stencilYMid);
		calculateCentralDifferenceStencil(Z, true, m_stencilZMid);

		calculateLaplacian(m_laplacian);

		m_positionsMid.resize(getGridMatrixSize(true), 3);
		for (int x = 0; x < m_dims(X) + 1; x++)
		{
			for (int y = 0; y < m_dims(Y) + 1; y++)
			{
				for (int z = 0; z < m_dims(Z) + 1; z++)
				{
					// faux index in the element-aligned matrix
					int i = z + y * m_dims(Z) + x * m_dims(Z) * m_dims(Y);
					// actual index in the position grid
					int iMid = z + y * (m_dims(Z) + 1) + x * (m_dims(Z) + 1) * (m_dims(Y) + 1);
					// Note that this function doesn't actually care if it's within bounds
					m_positionsMid.row(iMid) = getLocalCoordinatesOfElement(iMid, true);
				}
			}
		}
		
		
		// (n x m) (m x 3) -> (n x 3)
		m_midToElement.resize(getGridMatrixSize(false), getGridMatrixSize(true));
		Eigen::MatrixXd helper(getGridMatrixSize(false), getGridMatrixSize(true));
		spdlog::info("Helper Size = {} x {}", getGridMatrixSize(false), getGridMatrixSize(true));
		for (int x = 0; x < m_dims(X); x++)
		{
			for (int y = 0; y < m_dims(Y); y++)
			{
				// index in the element-aligned matrix
				int i = y * m_dims(Z) + x * m_dims(Z) * m_dims(Y);
				// index in the midgrid-aligned matrix
				int iMid = y * (m_dims(Z) + 1) + x * (m_dims(Z) + 1) * (m_dims(Y) + 1);

				spdlog::info("Setting m_midToElement ({}, {}) diagonal block", i, iMid);
				helper.block(i, iMid, m_dims(Z), m_dims(Z)) = Eigen::MatrixXd::Identity(m_dims(Z), m_dims(Z));
			}
		}
		m_midToElement = helper.sparseView();
		m_midToElement.makeCompressed();

		/* SET SIZES */
		m_velocity.resize(getGridMatrixSize(true), 3);
		m_pressure.resize(getGridMatrixSize(false), 1);
		m_dSignedDistance.resize(getGridMatrixSize(true), 3);
		m_signedDistance.resize(getGridMatrixSize(true), 1);
		spdlog::info("Finished Constructing EulerState");
	};

	const size_t EulerState::getGridMatrixSize(bool midGrid)
	{ 
		int size = 0;
		if (!midGrid)
		{
			size = m_dims.prod();
		}
		else
		{
			size = (m_dims(X) + 1) * (m_dims(Y) + 1) * (m_dims(Z) + 1);
		}
		return size;
	}

	void EulerState::calculateCentralDifferenceStencil(Dimension dim, bool midGrid, Eigen::SparseMatrix<double, Eigen::ColMajor>& stencil)
	{
		std::vector<Eigen::Triplet<double>> triplets;

		// If we want the central difference for a midpoint defined value, e.g. v_(1/2), we need a different grid size for the inputs
		if (midGrid)
		{
			stencil.resize(getGridMatrixSize(false), getGridMatrixSize(true));
		}
		else
		{
			stencil.resize(getGridMatrixSize(false), getGridMatrixSize(false));
		}

		int spacing = 1;

		// Indices stored as X * Y * Z (vertical)
		// spacing(Z) -> 1
		// spacing(Y) -> m_dim(Z)
		// spacing(X) -> m_dim(Z) * m_dim(Y)
		for (int i = Z; i > dim; i--)
		{
			spacing = spacing * (m_dims(i) + 1 * midGrid);
		}

		triplets.reserve(getGridMatrixSize(false) * (size_t)2);

		for (int x = 0; x < m_dims( X ); x++)
		{
			for (int y = 0; y < m_dims( Y ); y++)
			{
				for (int z = 0; z < m_dims( Z ); z++)
				{
					// Index
					int i = z + y * m_dims( Z ) + x * m_dims( Z ) * m_dims( Y );

					// If it's a mid-point grid location, the value on both sides are known
					if (midGrid)
					{
						triplets.emplace_back(i, i, -1.0);
						triplets.emplace_back(i, i + spacing, 1.0);
					}
					// Within bounds
					else if (i - spacing >= 0 && i < getGridMatrixSize(false))
					{

						triplets.emplace_back(i, i - spacing, 1.0);
						triplets.emplace_back(i, i, -1.0);
					}
					// Along an edge - no gradient (e.g. assume pressure is equalized on both sides of the grid cell)
					else if (i - spacing < 0)
					{
						//triplets.emplace_back(i, i, 
					} 
					else if (i + spacing >= getGridMatrixSize(false))
					{
						//triplets.emplace_back(i, i - spacing, -2.0);
					}
				}
			}
		}

		stencil.setFromTriplets(triplets.begin(), triplets.end());
	}

	void EulerState::calculateLaplacian(Eigen::SparseMatrix<double>& stencil)
	{
		std::vector<Eigen::Triplet<double>> triplets;
		int spacing[3];
		spacing[0] = m_dims(Z) * m_dims(Y);
		spacing[1] = m_dims(Z);
		spacing[2] = 1;

		stencil.resize(getGridMatrixSize(false), getGridMatrixSize(false));

		// Loop through all of the grid cells
		for (int x = 0; x < m_dims(X); x++)
		{
			for (int y = 0; y < m_dims(Y); y++)
			{
				for (int z = 0; z < m_dims(Z); z++)
				{
					int i = z + y * m_dims(Z) + x * m_dims(Z) * m_dims(Y);

					// For each grid cell, insert the central difference on each side for each dimension
					for (int dim = 0; dim < 3; dim++)
					{
						for (int dir = -1; dir <= 1; dir += 2)
						{
							// Only add this central difference if it is defined
							if (i + dir * spacing[dim] >= 0 && i + dir * spacing[dim] < getGridMatrixSize(false))
							{
								triplets.emplace_back(i, i, 1 / m_gridSizeHorizontal(dim));
								triplets.emplace_back(i, i + dir * spacing[dim], 1 / m_gridSizeHorizontal(dim));
							}
						}
					}
				}
			}
		}
		stencil.setFromTriplets(triplets.begin(), triplets.end());
	}

	void EulerState::getCentralDifferenceStencil(Dimension dim, bool midGrid, Eigen::SparseMatrix<double, Eigen::ColMajor>& stencil)
	{
		if (!midGrid)
		{
			switch (dim)
			{
			case X:
				stencil = m_stencilX;
				break;
			case Y:
				stencil = m_stencilY;
				break;
			case Z:
				stencil = m_stencilZ;
				break;
			}
		}
		else
		{
			switch (dim)
			{
			case X:
				stencil = m_stencilXMid;
				break;
			case Y:
				stencil = m_stencilYMid;
				break;
			case Z:
				stencil = m_stencilZMid;
				break;
			}
		}
	}

	void EulerState::getQuantityGradient(Eigen::SparseMatrix<double, Eigen::ColMajor>& dv, bool midGrid, const Eigen::SparseMatrix<double>& quantity)
	{
		//Eigen::SparseMatrix<double, Eigen::ColMajor> xStencil, yStencil, zStencil;
		//Eigen::SparseMatrix<double, Eigen::ColMajor> dv(3, getGridMatrixSize());

		dv.resize(getGridMatrixSize(false), 3);

		/*getCentralDifferenceStencil(X, midGrid, xStencil);
		getCentralDifferenceStencil(Y, midGrid, yStencil);
		getCentralDifferenceStencil(Z, midGrid, zStencil);*/

		//spdlog::info("Stencil = {} x {}", xStencil.rows(), xStencil.cols());
		//spdlog::info("quantity = {} x {}", quantity.rows(), quantity.cols());
		//std::cout << "Stencil = " << xStencil.rows() << " x " << xStencil.cols() << std::endl;
		//std::cout << "m_velocity = " << m_velocity.rows() << " x " << m_velocity.cols() << std::endl;
		//std::cout << "dv = " << dv.rows() << " x " << dv.cols() << std::endl;

		if (!midGrid)
		{
			dv.col(0) = m_stencilX * quantity / m_gridSizeHorizontal(0);
			dv.col(1) = m_stencilY * quantity / m_gridSizeHorizontal(1);
			dv.col(2) = m_stencilZ * quantity / m_gridSizeHorizontal(2);
		}
		else
		{
			dv.col(0) = m_stencilXMid * quantity.col(0) / m_gridSizeHorizontal(0);
			dv.col(1) = m_stencilYMid * quantity.col(1) / m_gridSizeHorizontal(1);
			dv.col(2) = m_stencilZMid * quantity.col(2) / m_gridSizeHorizontal(2);
		}
	}

	void EulerState::getPressureGradient(Eigen::SparseMatrix<double>& dp)
	{
		getQuantityGradient(dp, false, m_pressure);
	}

	void EulerState::getQuantityDivergence(Eigen::SparseMatrix<double>& dv, const Eigen::SparseMatrix<double> quantity)
	{
		Eigen::SparseMatrix<double> dvDir;

		dv.resize(getGridMatrixSize(false), getGridMatrixSize(false));

		getQuantityGradient(dvDir, true, quantity);

		Eigen::SparseVector<double, Eigen::ColMajor> ones(3, 1);
		ones.insert(0, 0) = 1;
		ones.insert(1, 0) = 1;
		ones.insert(2, 0) = 1;

		dv = dvDir * ones;
	}

	void EulerState::getLaplacianOperator(Eigen::SparseMatrix<double>& d2v)
	{
		d2v = m_laplacian;
	}

	Eigen::Vector3d EulerState::getLocalCoordinatesOfElement(size_t gridIndex, bool midGrid)
	{
		size_t x, y, z;
		z = gridIndex % (m_dims(Z) + midGrid);
		y = (gridIndex / (m_dims(Z) + midGrid)) % (m_dims(Y) + midGrid);
		x = gridIndex / (((size_t)m_dims(Y) + midGrid) * ((size_t)m_dims(Z) + midGrid));

		Eigen::Vector3d pos(x + 0.5, y + 0.5, z + 0.5);
		
		pos = pos.cwiseProduct(m_gridSizeHorizontal);

		return pos;
	}


}