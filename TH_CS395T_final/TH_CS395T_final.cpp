// TH_CS395T_final.cpp : Defines the entry point for the application.
//

#include "TH_CS395T_final.h"
#include <Eigen/Sparse>

using namespace std;

int main()
{
	Eigen::Vector3d a;
	a << 0.0, 1.0, 2.0;

	cout << a;

	cout << "Hello CMake." << endl;
	return 0;
}
