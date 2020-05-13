#include "ParticleVisualizer.h"
#include <igl/unproject.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <Eigen/Core>

namespace FluidVisualizer
{
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

	void ParticleVisualizer::resetSimulation()
	{
		if (!m_simulation)
			return;

		killSimThread();
		please_die = running = false;
		please_pause = true;
		//initSimulation();
		updateRenderGeometry();
		sim_thread = new std::thread(&ParticleVisualizer::runSimThread, this);
	}

	bool ParticleVisualizer::drawCallback(igl::opengl::glfw::Viewer& viewer)
	{
		if (!m_simulation)
			return false;

		//render(viewer);

		render_mutex.lock();
		
		viewer.data().clear();
		viewer.data().set_mesh(renderQ, renderF);
		viewer.data().add_points(renderP, Eigen::RowVector3d(0.19, 0.2, 1.0));
		viewer.data().set_colors(renderC);

		render_mutex.unlock();


		return false;
	}

	void ParticleVisualizer::updateRenderGeometry()
	{
		double radius = m_simulation->getState()->m_gridSizeHorizontal.minCoeff() / 2.0 / 2.0;
		int numcirclewedges = 10;

		std::vector<Eigen::Vector3d> verts;
		std::vector<Eigen::Vector3d> points;
		std::vector<Eigen::Vector3d> vertexColors;
		std::vector<Eigen::Vector3i> faces;

		int idx = 0;
		double eps = 1e-4;

		/* FLOOR */
		if (true)
		{
			for (int i = 0; i < 6; i++)
			{
				vertexColors.push_back(Eigen::Vector3d(0.3, 1.0, 0.3));
			}

			verts.push_back(Eigen::Vector3d(-1, -0.5, eps));
			verts.push_back(Eigen::Vector3d(1, -0.5, eps));
			verts.push_back(Eigen::Vector3d(-1, -1, eps));

			faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));

			verts.push_back(Eigen::Vector3d(-1, -1, eps));
			verts.push_back(Eigen::Vector3d(1, -0.5, eps));
			verts.push_back(Eigen::Vector3d(1, -1, eps));
			faces.push_back(Eigen::Vector3i(idx + 3, idx + 4, idx + 5));
			idx += 6;
		}

		int gridCells = m_simulation->getState()->getGridMatrixSize(false);


		for (int i = 0; i < gridCells; i++)
		{
			Eigen::Vector3d color(0, 0, 1.0);

			for (int scatterBalls = 0; scatterBalls < 6; scatterBalls++)
			{


				Eigen::Vector3d pos = m_simulation->getState()->getLocalCoordinatesOfElement(i, false);
				pos += Eigen::Vector3d::Random() / m_simulation->getState()->m_gridSizeHorizontal.minCoeff();

				points.push_back(Eigen::Vector3d(pos[0], pos[1], pos[2]));

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

		renderQ.resize(verts.size(), 3);
		renderC.resize(vertexColors.size(), 3);
		for (int i = 0; i < verts.size(); i++)
		{
			renderQ.row(i) = verts[i];
			renderC.row(i) = vertexColors[i];
		}
		renderF.resize(faces.size(), 3);
		for (int i = 0; i < faces.size(); i++)
			renderF.row(i) = faces[i];

		renderP.resize(points.size(), 3);
		for (int i = 0; i < points.size(); i++)
			renderP.row(i) = points[i];


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