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
		m_pressure.resize(getGridMatrixSize(false), 3);
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
					else if (i - spacing >= 0 && i + spacing < getGridMatrixSize(false))
					{

						triplets.emplace_back(i, i - spacing, -1.0);
						triplets.emplace_back(i, i + spacing, 1.0);
					}
					// Along an edge - no gradient (e.g. assume pressure is equalized on both sides of the grid cell)
					else if (i - spacing < 0)
					{
						//triplets.emplace_back(i, i + spacing, 2.0);
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
		Eigen::SparseMatrix<double, Eigen::ColMajor> xStencil, yStencil, zStencil;
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
			dv.col(0) = m_stencilX * quantity.col(0) / 2.0 / m_gridSizeHorizontal(0);
			dv.col(1) = m_stencilY * quantity.col(1) / 2.0 / m_gridSizeHorizontal(1);
			dv.col(2) = m_stencilZ * quantity.col(2) / 2.0 / m_gridSizeHorizontal(2);
		}
		else
		{
			dv.col(0) = m_stencilXMid * quantity.col(0) / 2.0 / m_gridSizeHorizontal(0);
			dv.col(1) = m_stencilYMid * quantity.col(1) / 2.0 / m_gridSizeHorizontal(1);
			dv.col(2) = m_stencilZMid * quantity.col(2) / 2.0 / m_gridSizeHorizontal(2);
		}
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