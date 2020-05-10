#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>

namespace FluidSimulation
{
	enum Dimension {X, Y, Z};

	class EulerState
	{
	public:
		EulerState(Eigen::Vector3i dims, Eigen::Vector3d gridSize);

		/* DATA MEMBERS */
		Eigen::VectorXd m_signedDistance;
		Eigen::MatrixX3d m_dSignedDistance;
		Eigen::SparseMatrix<double> m_velocity;
		Eigen::MatrixX3d m_positionsMid;
		Eigen::SparseMatrix<double> m_midToElement;
		const Eigen::Vector3i m_dims;
		const Eigen::Vector3d m_gridSizeHorizontal;


		void getWaterSurface(Eigen::MatrixXd& surface);
		void getGridStates(Eigen::MatrixXi& grid);

		void getPressureGrid(Eigen::SparseMatrix<double>& p);
		void getForcesGrid(Eigen::SparseMatrix<double>& f);

		void getPressureGradient(Eigen::SparseMatrix<double>& dp);
		void getQuantityGradient(Eigen::SparseMatrix<double, Eigen::ColMajor>& dv, bool midGrid, const Eigen::SparseMatrix<double>& quantity);

		const Eigen::Vector3i getDimensions() { return m_dims; };
		const size_t getGridMatrixSize(bool midGrid);

		void getCentralDifferenceStencil(Dimension dim, bool midGrid, Eigen::SparseMatrix<double, Eigen::ColMajor>& stencil);
		Eigen::Vector3d getLocalCoordinatesOfElement(size_t gridIndex, bool midGrid);

	private:
		void calculateCentralDifferenceStencil(Dimension dim, bool midGrid, Eigen::SparseMatrix<double, Eigen::ColMajor>& stencil);

		Eigen::Vector3d getVelocityAtPoint(const Eigen::Vector3d point);
		Eigen::Vector3d getPressureAtPoint(const Eigen::Vector3d point);
		
		//Eigen::SparseMatrix<double> m_velocity;
		Eigen::SparseMatrix<double> m_pressure;
		Eigen::SparseMatrix<double> m_forces;

		Eigen::SparseMatrix<double, Eigen::ColMajor> m_stencilX;
		Eigen::SparseMatrix<double, Eigen::ColMajor> m_stencilY;
		Eigen::SparseMatrix<double, Eigen::ColMajor> m_stencilZ;

		Eigen::SparseMatrix<double, Eigen::ColMajor> m_stencilXMid;
		Eigen::SparseMatrix<double, Eigen::ColMajor> m_stencilYMid;
		Eigen::SparseMatrix<double, Eigen::ColMajor> m_stencilZMid;
	};
}