#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>

namespace FluidSimulation
{
	enum Dimension {X, Y, Z};

	class EulerState
	{
	public:
		EulerState(Eigen::Vector3d dims) : m_dims(dims) { };

		void getWaterSurface(Eigen::MatrixXd& surface);
		void getGridStates(Eigen::MatrixXi& grid);

		void getVelocityGrid(Eigen::SparseMatrix<double>& v);
		void getPressureGrid(Eigen::SparseMatrix<double>& p);
		void getForcesGrid(Eigen::SparseMatrix<double>& f);

		void getPressureGradient(Eigen::SparseMatrix<double>& dp);
		void EulerState::getVelocityGradient(Eigen::SparseMatrix<double, Eigen::ColMajor>& dv);

		Eigen::Vector3i getDimensions() { return m_dims; };
		int getGridMatrixSize() { return m_dims.prod(); };

		void getCentralDifferenceStencil(Dimension dim, Eigen::SparseMatrix<int>& stencil);

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