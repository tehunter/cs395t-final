#pragma once

#include <Eigen/Sparse>

namespace FluidSimulation
{
	enum Dimension {X, Y, Z};

	class EulerState
	{
	public:
		EulerState(Eigen::Vector3d dims) : m_dims(dims) { };

		void getWaterSurface(Eigen::MatrixXd& surface);
		void getGridStates(Eigen::MatrixXi& grid);

		Eigen::SparseMatrix<double>& getVelocityGrid();
		Eigen::SparseMatrix<double>& getPressureGrid();
		Eigen::SparseMatrix<double>& getForcesGrid();

		Eigen::SparseMatrix<double>& getPressureGradient();
		Eigen::SparseMatrix<double>& getVelocityGradient();

		Eigen::Vector3i getDimensions() { return m_dims; };
		int getGridMatrixSize() { return m_dims.prod(); };

		void getCentralDifferenceStencil(Dimension dim, Eigen::SparseMatrix<int>& stencil)
	private:
		Eigen::Vector3d getVelocityAtPoint(const Eigen::Vector3d point);
		Eigen::Vector3d getPressureAtPoint(const Eigen::Vector3d point);
		
		Eigen::SparseMatrix<double> m_velocity;
		Eigen::SparseMatrix<double> m_pressure;
		Eigen::SparseMatrix<double> m_forces;

		const Eigen::Vector3i m_dims;
		// deltaX, deltaY, deltaZ
		const Eigen::Vector3d m_gridSizeHorizontal;
	};
}