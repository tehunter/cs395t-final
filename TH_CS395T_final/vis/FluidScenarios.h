#pragma once

#include "../lib/EulerSimulation.h"

namespace FluidVisualizer
{
	enum Scenario { Hydrostatic, Sphere, SuspendedColumn, DamBreak };

	void setupScenario(FluidSimulation::EulerSimulation* m_simulation, Scenario m_scenario);
}