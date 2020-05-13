#pragma once

#include <mutex>
#include <thread>
#include "../lib/EulerSimulation.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>



namespace FluidVisualizer
{

	class ParticleVisualizer
	{
	public:
		ParticleVisualizer(FluidSimulation::EulerSimulation* simulation) : m_simulation(simulation) { };
		virtual ~ParticleVisualizer()
		{
			killSimThread();
		}
		bool drawGUI(igl::opengl::glfw::imgui::ImGuiMenu& menu);
		void mouseClicked(double x, double y, int button);
		bool mouseScroll(igl::opengl::glfw::Viewer& viewer, float delta);
		bool mouseCallback(igl::opengl::glfw::Viewer& viewer, int button, int modifier);
		bool keyCallback(igl::opengl::glfw::Viewer& viewer, unsigned int key, int modifiers);
		bool drawCallback(igl::opengl::glfw::Viewer& viewer);
		void resetSimulation();
		void toggleSimulation();

		void updateRenderGeometry();

		virtual void tick() {};
		//void render(igl::opengl::glfw::Viewer& viewer);
		bool isPaused();
		void pause();
		void reset();
		void run();

	protected:
		Eigen::MatrixXd renderQ;
		Eigen::MatrixXi renderF;
		Eigen::MatrixXd renderC;
		Eigen::MatrixXd renderP;

		void runSimThread();
		void killSimThread();

		std::thread* sim_thread;
		bool please_pause;
		bool please_die;
		bool running;
		std::mutex render_mutex;
		std::mutex status_mutex;
		std::mutex message_mutex;

	private:
		FluidSimulation::EulerSimulation* m_simulation;
	};
}