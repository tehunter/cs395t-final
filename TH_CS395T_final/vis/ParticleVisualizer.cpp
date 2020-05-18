#include "ParticleVisualizer.h"
#include <igl/unproject.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <Eigen/Core>
#include <thread>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "FluidScenarios.h"

enum Dimension{X, Y, Z};

namespace FluidVisualizer
{
	ParticleVisualizer::ParticleVisualizer(FluidSimulation::EulerSimulation* simulation) : 
		m_simulation(simulation), sim_thread(NULL), please_pause(false), please_die(false), please_stepOnce(false), running(false), m_scenario(Hydrostatic)
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
	void ParticleVisualizer::stepOnce()
	{
		if (!m_simulation)
			return;

		if (!isPaused())
			return;

		status_mutex.lock();
		please_stepOnce = true;
		status_mutex.unlock();
	}
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

		m_simulation->updateState()->m_velocity.setZero();

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


		
		setupScenario(m_simulation, m_scenario);
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

		int gridCells = m_simulation->getState()->getGridMatrixSize(true);

		Eigen::Vector3d offsets[8];
		offsets[0] << 0.25, 0.25, 0.25;
		offsets[1] << 0.25, 0.25, -0.25;
		offsets[2] << 0.25, -0.25, 0.25;
		offsets[3] << -0.25, 0.25, 0.25;
		offsets[4] << 0.25, -0.25, -0.25;
		offsets[5] << -0.25, -0.25, 0.25;
		offsets[6] << -0.25, 0.25, -0.25;
		offsets[7] << -0.25, -0.25, -0.25;

		bool refinedParticles = false;

		for (int i = 0; i < gridCells; i++)
		{
			Eigen::Vector3d color(0, 0, 1.0);
			
			if (i > (m_simulation->getState()->getGridMatrixSize(true) - 1 - m_simulation->getState()->m_dims(2)))
				continue;

			int z = (i % (m_simulation->getState()->m_dims(2) + 1));
			int y = ((i / (m_simulation->getState()->m_dims(2) + 1)) % (m_simulation->getState()->m_dims(1) + 1));
			int x = (i / ((m_simulation->getState()->m_dims(2) + 1) * (m_simulation->getState()->m_dims(1) + 1)));

			if (z == m_simulation->getState()->m_dims(Z) || y == m_simulation->getState()->m_dims(Y) || x == m_simulation->getState()->m_dims(X))
				continue;


			for (int scatterBalls = 0; scatterBalls < 1; scatterBalls++)
			{

				
				Eigen::Vector3d pos = m_simulation->getState()->getLocalCoordinatesOfElement(i, true);
				if (refinedParticles)
				{
					for (int j = 0; j < 8; j++)
					{
						Eigen::Vector3d point = pos + offsets[j].cwiseProduct(m_simulation->getState()->m_gridSizeHorizontal);

						double interpZ, interpZ_Y, interpZ_X, interpZ_XY;
						interpZ = interpZ_Y = interpZ_X = interpZ_XY = m_simulation->getState()->m_signedDistance(i);

						if ((i + (offsets[j](2) > 0 ? 1 : -1)) < gridCells && (i + (offsets[j](2) > 0 ? 1 : -1)) >= 0)
							interpZ = m_simulation->getState()->m_signedDistance(i) * 0.75 + m_simulation->getState()->m_signedDistance(i + (offsets[j](2) > 0 ? 1 : -1)) * 0.25;

						if ((i + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1)) < gridCells &&
								((i + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1))) >= 0 &&
								(i + (offsets[j](2) > 0 ? 1 : -1) + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1)) < gridCells &&
								(i + (offsets[j](2) > 0 ? 1 : -1) + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1)) >= 0)
							interpZ_Y = m_simulation->getState()->m_signedDistance(i + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1)) * 0.75 + m_simulation->getState()->m_signedDistance(i + (offsets[j](2) > 0 ? 1 : -1) + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1)) * 0.25;

						if ((i + (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) < gridCells &&
								(i + (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) >= 0 &&
								(i + (offsets[j](2) > 0 ? 1 : -1) + (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) < gridCells &&
								(i + (offsets[j](2) > 0 ? 1 : -1) + (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) >= 0)
							interpZ_X = m_simulation->getState()->m_signedDistance(i + (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) * 0.75 + m_simulation->getState()->m_signedDistance(i + (offsets[j](2) > 0 ? 1 : -1) + (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) * 0.25;

						if ((i + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) + (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) < gridCells &&
								(i + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) + (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) >= 0 &&
								(i + (offsets[j](2) > 0 ? 1 : -1) + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) < gridCells &&
								(i + (offsets[j](2) > 0 ? 1 : -1) + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) >= 0)
							interpZ_XY = m_simulation->getState()->m_signedDistance(i + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) + (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) * 0.75 + m_simulation->getState()->m_signedDistance(i + (offsets[j](2) > 0 ? 1 : -1) + (offsets[j](1) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (offsets[j](0) > 0 ? 1 : -1) * (m_simulation->getState()->m_dims(Z) + 1) * (m_simulation->getState()->m_dims(Y) + 1)) * 0.25;

						double Y0 = interpZ * 0.75 + interpZ_Y * 0.25;
						double Y1 = interpZ_X * 0.75 + interpZ_XY * 0.25;

						double dist = Y0 * 0.75 + Y1 * 0.25;

						if (dist <= 0)
						{
							//pos += Eigen::Vector3d::Random() * m_simulation->getState()->m_gridSizeHorizontal.minCoeff() / 10.0;
							points.push_back(point);
							colors.push_back(Eigen::Vector3d(0.8, 0.0, 0.0) * std::min(m_simulation->getState()->m_velocity.row(i).norm() / 3.0, 1.0) + Eigen::Vector3d(0.2, 0.2, 0.2));
						}
					}
				}
				else
				{
					if (m_simulation->getState()->m_signedDistance(i) <= 0)
					{
						//pos += Eigen::Vector3d::Random() * m_simulation->getState()->m_gridSizeHorizontal.minCoeff() / 10.0;
						points.push_back(pos);
						colors.push_back(Eigen::Vector3d(0.8, 0.0, 0.0) * std::min(m_simulation->getState()->m_velocity.row(i).norm() / 3.0, 1.0) + Eigen::Vector3d(0.2, 0.2, 0.2));
					}
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
			if (ImGui::Button("Step Once", ImVec2(-1, 0)))
			{
				stepOnce();
			}
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
			std::string strings[4] = { "Hydrostatic Volume", "Falling Sphere", "Suspended Column", "Dam Break"};
			for (int i = 0; i < 4; i++)
			{
				if (ImGui::RadioButton(strings[i].c_str(), m_scenario == i))
				{
					m_scenario = (Scenario)i;
				}
			}
		}
		if (ImGui::CollapsingHeader("Simulation Settings", ImGuiTreeNodeFlags_DefaultOpen))
		{
			bool gravity = true, pressure = true, recalcLevelSet = false;
			double cutoff = 3;
			if (m_simulation)
			{
				gravity = m_simulation->getGravityEnabled();
				pressure = m_simulation->getPressureEnabled();
				recalcLevelSet = m_simulation->getRecalculateLevelSet();
				cutoff = m_simulation->m_cutoffDistance;
			}
			if (ImGui::InputDouble("Cutoff Distance", &cutoff, 0.5, 0.5, "%.2f"))
			{
				m_simulation->m_cutoffDistance = cutoff;
			}

			if (ImGui::Checkbox("Gravity", &gravity))
			{
				if (m_simulation && gravity != m_simulation->getGravityEnabled())
				{
					m_simulation->toggleGravity();
				}
				//gravity = !gravity;
			}
			if (ImGui::Checkbox("Pressure", &pressure))
			{
				if (m_simulation && pressure != m_simulation->getPressureEnabled())
				{
					m_simulation->togglePressure();
				}
				//pressure = !pressure;
			}
			if (ImGui::Checkbox("Recalculate Level Set", &recalcLevelSet))
			{
				if (m_simulation && recalcLevelSet != m_simulation->getRecalculateLevelSet())
				{
					m_simulation->toggleRecalculateLevelSet();
				}
				//recalcLevelSet = !recalcLevelSet;
			}
		}
		if (ImGui::CollapsingHeader("Miscellaneous", ImGuiTreeNodeFlags_DefaultOpen))
		{
			bool logDebug = spdlog::default_logger()->level() == spdlog::level::debug;
			if (ImGui::Checkbox("Log Level = Debug", &logDebug))
			{
				//logDebug = !logDebug;
				if (logDebug)
					spdlog::set_level(spdlog::level::debug);
				else
					spdlog::set_level(spdlog::level::info);
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
			bool steponce = please_stepOnce;
			status_mutex.unlock();
			if (pausenow && !please_stepOnce)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
			}
			else if (please_stepOnce)
			{
				m_simulation->step(0.01);
				status_mutex.lock();
				please_stepOnce = false;
				status_mutex.unlock();
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