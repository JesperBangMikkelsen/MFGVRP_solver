#include "MFGVRP_Solver.h"
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
int main() {
	std::string path = "MFGVRP/";

	
	for (const auto& entry : fs::directory_iterator(path)) { 
		MFGVRP_Solver solver = MFGVRP_Solver();

		std::cout << entry.path().string() << "\n";
		solver.LoadData(entry.path().string(), true);

	}
	//MFGVRP_Solver solver = MFGVRP_Solver();
	//try {
	//	solver.LoadData("VRP Data - BIG - EVRP - N50.txt", true);
	//}
	//catch (IloException& e) {
	//	std::cerr << "Concert Exception: " << e << std::endl;
	//}
	//
	return 0;
}
