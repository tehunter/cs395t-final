// TH_CS395T_final.cpp : Defines the entry point for the application.
//

#include "TH_CS395T_final.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <thread>
#include <Eigen/Sparse>
#include "lib/EulerSimulation.h"
#include "vis/ParticleVisualizer.h"



static FluidSimulation::EulerSimulation* simulation = NULL;
static FluidVisualizer::ParticleVisualizer* visualizer = NULL;

using namespace std;


int main()
{
	igl::opengl::glfw::Viewer viewer;

	Eigen::Vector3i dimSize(10, 10, 20);
	Eigen::Vector3d gridWidth(0.2, 0.2, 0.2);

	simulation = new FluidSimulation::EulerSimulation(dimSize, gridWidth);
	visualizer = new FluidVisualizer::ParticleVisualizer(simulation);

	viewer.data().show_lines = true;
	viewer.data().set_face_based(false);
	viewer.core().is_animating = true;
	viewer.callback_key_pressed = [](igl::opengl::glfw::Viewer& viewer, unsigned int key, int modifiers) { return visualizer->keyCallback(viewer, key, modifiers); };
	viewer.callback_pre_draw = [](igl::opengl::glfw::Viewer& viewer) { return visualizer->drawCallback(viewer); };
	viewer.callback_mouse_down = [](igl::opengl::glfw::Viewer& viewer, int button, int modifier) { return visualizer->mouseCallback(viewer, button, modifier); };
	//viewer.callback_mouse_scroll = &(visualizer.mouseScroll);

	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]() {visualizer->drawGUI(menu); };

	viewer.launch();


	Eigen::Vector3d a;
	a << 0.0, 1.0, 2.0;

	cout << a;

	cout << "Hello CMake." << endl;
	return 0;
}
