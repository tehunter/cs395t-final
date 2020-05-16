#include "ParticleVisualizer.h"
#include <igl/unproject.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <Eigen/Core>
#include <thread>

enum Dimension{X, Y, Z};

namespace FluidVisualizer
{
	ParticleVisualizer::ParticleVisualizer(FluidSimulation::EulerSimulation* simulation) : 
		m_simulation(simulation), sim_thread(NULL), please_pause(false), please_die(false), running(false), m_scenario(Hydrostatic)
	{
		/* FLOOR */
		if (true)
		{
			double floory = 0.0;
			renderQ.resize(5, 3);
			renderC.resize(5, 3);
			// floor
			renderQ.row(0) << 0, 0, floory; renderC.row(0) << 0.5, 0.5, 0.5;
			renderQ.row(1) << 1e6, 1e6, floory; renderC.row(1) << 0.5, 0.5, 0.5;
			renderQ.row(2) << -1e6, 1e6, floory; renderC.row(2) << 0.5, 0.5, 0.5;
			renderQ.row(3) << -1e6, -1e6, floory; renderC.row(3) << 0.5, 0.5, 0.5;
			renderQ.row(4) << 1e6, -1e6, floory; renderC.row(4) << 0.5, 0.5, 0.5;

			renderF.resize(4, 3);
			renderF.row(0) << 0, 2, 1;
			renderF.row(1) << 0, 3, 2;
			renderF.row(2) << 0, 4, 3;
			renderF.row(3) << 0, 1, 4;
		}
	}
	void ParticleVisualizer::toggleSimulation()
	{
		if (!m_simulation)
			return;

		if (isPaused())
			run();
		else
			pause();
	}

	/*
	 * Runs the simulation, if it has been paused (or never started).
	 */
	void ParticleVisualizer::run()
	{
		status_mutex.lock();
		please_pause = false;
		status_mutex.unlock();
	}
	void ParticleVisualizer::pause()
	{
		status_mutex.lock();
		please_pause = true;
		status_mutex.unlock();
	}
	bool ParticleVisualizer::isPaused()
	{
		bool ret = false;
		status_mutex.lock();
		if (running && please_pause)
			ret = true;
		status_mutex.unlock();
		return ret;
	}

	void ParticleVisualizer::initSimulation()
	{
		int nverts = 5, nfaces = 4;

		/* WALLS */
		if (true)
		{
			renderQ.conservativeResize(nverts +6, 3);
			renderC.conservativeResize(nverts +6, 3);
			renderQ.row(nverts) << 0, 0, 0; renderC.row(nverts) << 0.1, 0.1, 0.1;
			renderQ.row(nverts + 1) << m_simulation->getState()->m_dims(0) * m_simulation->getState()->m_gridSizeHorizontal(0), 0, 0; renderC.row(nverts + 1) << 0.1, 0.1, 0.1;
			renderQ.row(nverts + 2) << m_simulation->getState()->m_dims(0) * m_simulation->getState()->m_gridSizeHorizontal(0), 0, m_simulation->getState()->m_dims(2)* m_simulation->getState()->m_gridSizeHorizontal(2); renderC.row(nverts + 2) << 0.1, 0.1, 0.1;
			renderQ.row(nverts + 3) << 0, m_simulation->getState()->m_dims(1) * m_simulation->getState()->m_gridSizeHorizontal(1), 0; renderC.row(nverts + 3) << 0.1, 0.1, 0.1;
			renderQ.row(nverts + 4) << 0, m_simulation->getState()->m_dims(1)* m_simulation->getState()->m_gridSizeHorizontal(1), m_simulation->getState()->m_dims(2)* m_simulation->getState()->m_gridSizeHorizontal(2); renderC.row(nverts + 4) << 0.1, 0.1, 0.1;
			renderQ.row(nverts + 5) << 0, 0, m_simulation->getState()->m_dims(2) * m_simulation->getState()->m_gridSizeHorizontal(2); renderC.row(nverts + 5) << 0.1, 0.1, 0.1;
			

			renderF.conservativeResize(nfaces + 4, 3);
			// x-z
			renderF.row(nfaces) << nverts, nverts + 5, nverts + 1;
			// y-z
			renderF.row(nfaces+1) << nverts, nverts + 5, nverts + 3;
			// x-z
			renderF.row(nfaces+2) << nverts + 5, nverts + 2, nverts + 1;
			// y-z
			renderF.row(nfaces+3) << nverts + 5, nverts + 3, nverts + 4;
			nfaces += 4;
			nverts += 6;
		}

		Eigen::Vector3d sphereCenter;
		sphereCenter = m_simulation->getState()->m_gridSizeHorizontal.cwiseProduct(m_simulation->getState()->m_dims.cast<double>()) / 2.0;
		double radius = sphereCenter(0) / 3.0;
		for (int x = 0; x < m_simulation->getState()->m_dims(1) + 1; x++)
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

	void ParticleVisualizer::resetSimulation()
	{
		if (!m_simulation)
			return;

		killSimThread();
		please_die = running = false;
		please_pause = true;
		initSimulation();
		updateRenderGeometry();
		sim_thread = new std::thread(&ParticleVisualizer::runSimThread, this);
	}

	bool ParticleVisualizer::drawCallback(igl::opengl::glfw::Viewer& viewer)
	{
		if (!m_simulation)
			return false;

		//render(viewer);

		render_mutex.lock();

		Eigen::Matrix3d p;
		p << 1, 0, 0, 0, 0, 1, 0, 1, 0;
		
		viewer.data(0).clear();
		viewer.data(1).clear();
		viewer.data(0).set_mesh(renderQ * p, renderF);
		viewer.data(0).set_colors(renderC);
		if (renderP.rows() > 0)
			viewer.data(1).add_points(renderP * p, pointColors);

		render_mutex.unlock();


		return false;
	}

	void ParticleVisualizer::updateRenderGeometry()
	{
		double radius = m_simulation->getState()->m_gridSizeHorizontal.minCoeff() / 2.0 / 2.0;
		int numcirclewedges = 10;

		std::vector<Eigen::Vector3d> verts;
		std::vector<Eigen::Vector3d> points, colors;
		std::vector<Eigen::Vector3d> vertexColors;
		std::vector<Eigen::Vector3i> faces;

		int idx = 0;
		double eps = 1e-4;

		int gridCells = m_simulation->getState()->getGridMatrixSize(false);


		for (int i = 0; i < gridCells; i++)
		{
			Eigen::Vector3d color(0, 0, 1.0);

			for (int scatterBalls = 0; scatterBalls < 1; scatterBalls++)
			{


				Eigen::Vector3d pos = m_simulation->getState()->getLocalCoordinatesOfElement(i, true);
				if (m_simulation->getState()->m_signedDistance(i) <= 0)
				{
					//pos += Eigen::Vector3d::Random() * m_simulation->getState()->m_gridSizeHorizontal.minCoeff() / 10.0;
					points.push_back(Eigen::Vector3d(pos[0], pos[1], pos[2]));
					colors.push_back(Eigen::Vector3d(0.8, 0.0, 0.0) * std::min(m_simulation->getState()->m_velocity.row(i).norm() / 3.0, 1.0) + Eigen::Vector3d(0.2, 0.2, 0.2));
				}

				const double PI = 3.1415926535898;
				/*for (int j = 0; j <= numcirclewedges; j++)
				{
					for (int k = 0; k < numcirclewedges; k++)
					{
						verts.push_back(Eigen::Vector3d(pos[0] + radius * cos(2 * PI * j / numcirclewedges) * sin(PI * k / numcirclewedges - PI / 2.0),
							pos[1] + radius * sin(2 * PI * j / numcirclewedges) * sin(PI * k / numcirclewedges - PI / 2.0), 
							pos[2] + radius * cos(PI * k / numcirclewedges - PI / 2.0)));
					}
				}*/

				/*for (int j = 0; j <= numcirclewedges; j++)
				{
					for (int k = 0; k <= numcirclewedges; k++)
					{
						faces.push_back(Eigen::Vector3i(idx, idx + k + 1, idx + k * numcirclewedges * 1 + ((j + 1) % (numcirclewedges + 1))));
					}
				}*/

				//idx += numcirclewedges + 2;
			}
		}

		//renderQ.resize(verts.size(), 3);
		//renderC.resize(vertexColors.size(), 3);
		//for (int i = 0; i < verts.size(); i++)
		//{
			//renderQ.row(i) = verts[i];
			//renderC.row(i) = vertexColors[i];
		//}
		//renderF.resize(faces.size(), 3);
		//for (int i = 0; i < faces.size(); i++)
			//renderF.row(i) = faces[i];

		renderP.setZero();
		pointColors.setZero();
		renderP.resize(points.size(), 3);
		pointColors.resize(points.size(), 3);
		for (int i = 0; i < points.size(); i++)
		{
			renderP.row(i) = points[i].transpose();
			pointColors.row(i) = colors[i].transpose();
		}
	}

	bool ParticleVisualizer::keyCallback(igl::opengl::glfw::Viewer& viewer, unsigned int key, int modifiers)
	{
		if (key == ' ')
		{
			toggleSimulation();
			return true;
		}
		Eigen::Vector4f look4 = viewer.core().view.inverse() * Eigen::Vector4f(0, 0, 1.0, 0.0);
		Eigen::Vector4f left4 = viewer.core().view.inverse() * Eigen::Vector4f(1.0, 0.0, 0.0, 0.0);
		Eigen::Vector3f look(look4[0], look4[1], look4[2]);
		Eigen::Vector3f left(left4[0], left4[1], left4[2]);
		if (key == 'w')
		{
			viewer.core().camera_base_translation += look;
		}
		if (key == 's')
		{
			viewer.core().camera_base_translation -= look;
		}
		if (key == 'a')
		{
			viewer.core().camera_base_translation += left;
		}
		if (key == 'd')
		{
			viewer.core().camera_base_translation -= left;
		}
		return false;
	}

	bool ParticleVisualizer::mouseCallback(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
	{
		/*Eigen::Vector3f pos(viewer.down_mouse_x, viewer.core.viewport[3] - viewer.down_mouse_y, 1);
		Eigen::Matrix4f model = viewer.core.view;
		Eigen::Vector3f unproj = igl::unproject(pos, model, viewer.core.proj, viewer.core.viewport);
		Eigen::Vector4f eye = viewer.core.view.inverse() * Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
		Eigen::Vector3d dir;
		for (int i = 0; i < 3; i++)
			dir[i] = unproj[i] - eye[i];
		dir.normalize();*/
		return false;
		//return mouseClicked(viewer, dir, button);
	}

	void ParticleVisualizer::mouseClicked(double x, double y, int button)
	{
		return;
	}

	bool mouseScroll(igl::opengl::glfw::Viewer& viewer, float delta)
	{
		return true;
	}

	bool ParticleVisualizer::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu& menu)
	{
		if (ImGui::CollapsingHeader("Simulation Control", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::Button("Run/Pause Sim", ImVec2(-1, 0)))
			{
				toggleSimulation();
			}
			if (ImGui::Button("Reset Sim", ImVec2(-1, 0)))
			{
				resetSimulation();
			}
		}
		if (ImGui::CollapsingHeader("Simulation Scenario", ImGuiTreeNodeFlags_DefaultOpen))
		{
			std::string strings[2] = { "Hydrostatic Volume", "Falling Sphere" };
			for (int i = 0; i < 2; i++)
			{
				if (ImGui::RadioButton(strings[i].c_str(), m_scenario == i))
				{
					m_scenario = (Scenario)i;
				}
			}
		}
		//m_simulation->drawGUI(menu);
		return false;
	}

	void ParticleVisualizer::runSimThread()
	{
		status_mutex.lock();
		running = true;
		status_mutex.unlock();

		bool done = false;
		while (!done)
		{
			tick();

			status_mutex.lock();
			bool pausenow = please_pause;
			status_mutex.unlock();
			if (pausenow)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
			}
			else
			{
				m_simulation->step(0.01);
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
			}
			render_mutex.lock();
			updateRenderGeometry();
			render_mutex.unlock();
			status_mutex.lock();
			if (please_die)
				done = true;
			status_mutex.unlock();
		}

		status_mutex.lock();
		running = false;
		status_mutex.unlock();
	}

	void ParticleVisualizer::killSimThread()
	{
		if (sim_thread)
		{
			status_mutex.lock();
			please_die = true;
			status_mutex.unlock();
			sim_thread->join();
			delete sim_thread;
			sim_thread = NULL;
		}
	}

}