#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "EulerState.h"

namespace FluidSimulation
{
	void EulerState::getCentralDifferenceStencil(Dimension dim, Eigen::SparseMatrix<int>& stencil)
	{
		std::vector<Eigen::Triplet<int>> triplets;

		// The stencil is only as long as the desired direction, in the other dimensions it just repeats
		stencil.resize(getGridMatrixSize(), getGridMatrixSize());

		int spacing = 1;

		// Indices stored as X * Y * Z (vertical)
		// spacing(Z) -> 1
		// spacing(Y) -> m_dim(Z)
		// spacing(X) -> m_dim(Z) * m_dim(Y)
		for (int i = Z; i > dim; i--)
		{
			spacing = spacing * m_dims(i);
		}

		triplets.reserve(getGridMatrixSize());
		for (int x = 0; x < m_dims(X); x++)
		{
			for (int y = 0; y < m_dims(Y); y++)
			{
				for (int z = 0; z < m_dims(Z); z++)
				{
					// Index
					int i = z + y * m_dims(Z) + x * m_dims(Z) * m_dims(Y);

					// Within bounds
					if (i - spacing >= 0 && i + spacing < getGridMatrixSize())
					{
						triplets.emplace_back(i, i - spacing, -1);
						triplets.emplace_back(i, i + spacing, 1);
					}
					// Along an edge
					else if (i - spacing < 0)
					{
						triplets.emplace_back(i, i + spacing, 2);
					} 
					else if (i + spacing >= getGridMatrixSize())
					{
						triplets.emplace_back(i, i - spacing, -2);
					}
				}
			}
		}

		stencil.setFromTriplets(triplets.begin(), triplets.end());
	}

	void EulerState::getVelocityGradient(Eigen::SparseMatrix<double, Eigen::ColMajor>& dv)
	{
		Eigen::SparseMatrix<int> xStencil, yStencil, zStencil;
		//Eigen::SparseMatrix<double, Eigen::ColMajor> dv(3, getGridMatrixSize());

		dv.resize(3, getGridMatrixSize());

		getCentralDifferenceStencil(X, xStencil);
		getCentralDifferenceStencil(Y, yStencil);
		getCentralDifferenceStencil(Z, zStencil);

		dv.col(0) = xStencil.cast<double>() * m_velocity / 2.0 / m_gridSizeHorizontal(0);
		dv.col(1) = yStencil.cast<double>() * m_velocity / 2.0 / m_gridSizeHorizontal(1);
		dv.col(2) = zStencil.cast<double>() * m_velocity / 2.0 / m_gridSizeHorizontal(2);
	}
}