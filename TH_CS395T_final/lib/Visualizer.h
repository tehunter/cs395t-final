#pragma once

#include <Eigen/Sparse>
#include "EulerState.h"

namespace FluidSimulation
{
	class Visualizer
	{
	public:
		void getMesh(const EulerState& state);
		void getLevelSet(const EulerState& state);
	private:
		
	};
}
