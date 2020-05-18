#include "FluidScenarios.h"
#include "../lib/EulerSimulation.h"
#include <Eigen/Core>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"


namespace FluidVisualizer
{
	void setupScenario(FluidSimulation::EulerSimulation* m_simulation, Scenario m_scenario)
	{
		Eigen::Vector3d sphereCenter;
		sphereCenter = m_simulation->getState()->m_gridSizeHorizontal.cwiseProduct(m_simulation->getState()->m_dims.cast<double>()) / 2.0;
		double radius = sphereCenter(0) / 3.0;
		for (int x = 0; x < m_simulation->getState()->m_dims(X) + 1; x++)
		{
			for (int y = 0; y < m_simulation->getState()->m_dims(Y) + 1; y++)
			{
				for (int z = 0; z < m_simulation->getState()->m_dims(Z) + 1; z++)
				{
					int i = z + y * (m_simulation->getState()->m_dims(Z) + 1) + x * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1);
					Eigen::Vector3d position(x + 0.5, y + 0.5, z + 0.5);
					position = position.cwiseProduct(m_simulation->getState()->m_gridSizeHorizontal);
					Eigen::Vector3d direction;
					double distance = 999;

					switch (m_scenario)
					{
					case Hydrostatic:
						direction = Eigen::Vector3d(0, 0, 1.0);
						distance = position(2) - sphereCenter(2);
						break;
					case Sphere:
						direction = (position - sphereCenter);
						if (direction.norm() > 0)
							distance = (direction).norm() - radius;
						else
							distance = -radius;
						break;
					case SuspendedColumn:
						/*if (x == m_simulation->getState()->m_dims(X))
						{
							direction = Eigen::Vector3d(0.5 * m_simulation->getState()->m_gridSizeHorizontal(X), 0, std::max(sphereCenter(2), position(2)));
							distance = direction.norm();
							direction.normalize();
						}
						else if (y == m_simulation->getState()->m_dims(X))
						{
							direction = Eigen::Vector3d(0, 0.5 * m_simulation->getState()->m_gridSizeHorizontal(Y), std::max(sphereCenter(2), position(2)));
							distance = direction.norm();
							direction.normalize();
						}
						else */
						if (z == m_simulation->getState()->m_dims(Z))
						{
							direction = direction = Eigen::Vector3d(0, 0, 1.0);
							distance = position(2) - m_simulation->getState()->m_dims(Z) * m_simulation->getState()->m_gridSizeHorizontal(Z);
						}
						else
						{
							bool closerToTop = (position(2) - sphereCenter(2)) >= m_simulation->getState()->m_dims(Z) * m_simulation->getState()->m_gridSizeHorizontal(Z) / 4.0;
							direction = direction = Eigen::Vector3d(0, 0, closerToTop ? 1.0 : -1.0);
							if (closerToTop)
								distance = position(2) - m_simulation->getState()->m_dims(Z) * m_simulation->getState()->m_gridSizeHorizontal(Z);
							else
								distance = sphereCenter(2) - position(2);
						}
						break;

					case DamBreak:
						if (y == m_simulation->getState()->m_dims(Y))
						{
							direction = Eigen::Vector3d(std::max(position(0) - sphereCenter(0), 0.0), 0.5 * m_simulation->getState()->m_gridSizeHorizontal(Y), 0);
							distance = direction.norm();
							direction.normalize();
						}
						else if (z == m_simulation->getState()->m_dims(Z))
						{
							direction = Eigen::Vector3d(std::max(position(0) - sphereCenter(0), 0.0), 0, 0.5 * m_simulation->getState()->m_gridSizeHorizontal(Z));
							distance = direction.norm();
							direction.normalize();
						}
						else if (x < (m_simulation->getState()->m_dims(X) / 2))
						{
							bool closerToTop = (m_simulation->getState()->m_dims(Z) * m_simulation->getState()->m_gridSizeHorizontal(Z) - (position(2))) < (m_simulation->getState()->m_dims(X) * m_simulation->getState()->m_gridSizeHorizontal(X) / 2.0 - (position(0)));
							if (closerToTop)
							{
								direction = Eigen::Vector3d(0, 0, 1.0);
								distance = -(m_simulation->getState()->m_dims(Z) * m_simulation->getState()->m_gridSizeHorizontal(Z) - (position(2)));
							}
							else
							{
								direction = Eigen::Vector3d(1.0, 0, 0);
								distance = position(0) - sphereCenter(0);
							}
							direction.normalize();
						}
						else
						{
							direction = Eigen::Vector3d(1.0, 0, 0);
							distance = position(0) - sphereCenter(0);
							spdlog::debug("({}, {}, {}) Position = {}, Center = {}", x, y, z, position, sphereCenter);
						}
						break;
					}

					// Set component 0 to the distance
					m_simulation->updateState()->m_signedDistance(i) = distance;
					// Set components 1-3 to the direction
					if (direction.norm() >= 0)
					{
						m_simulation->updateState()->m_dSignedDistance.row(i) = direction.transpose().normalized();
					}
					else
					{
						// If it's equidistant than pick a random direction
						m_simulation->updateState()->m_dSignedDistance.row(i) = Eigen::Vector3d::Random().transpose();
					}
				}
			}
		}
	}
}