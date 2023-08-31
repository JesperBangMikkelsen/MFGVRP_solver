//Including header file:
#ifndef _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include <iostream>
#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include "valarray"
#include <numeric>
#include <list>
#include <set>


////JENS CVRPSEP:
#include "CVRPSEP/CNSTRMGR.H"
#include "CVRPSEP/BASEGRPH.H"
#include "CVRPSEP/BINPACK.H"
#include "CVRPSEP/blocks.h"
#include "CVRPSEP/brnching.h"
#include "CVRPSEP/CAPSEP.H"
#include "CVRPSEP/combsep.h"
#include "CVRPSEP/COMPCUTS.H"
#include "CVRPSEP/COMPRESS.H"
#include "CVRPSEP/CUTBASE.H"
#include "CVRPSEP/FCAPFIX.H"
#include "CVRPSEP/FCISEP.H"
#include "CVRPSEP/fcits.h"
#include "CVRPSEP/glmsep.h"
#include "CVRPSEP/Grsearch.h"
#include "CVRPSEP/HPMSTAR.H"
#include "CVRPSEP/HTOURSEP.H"
#include "CVRPSEP/intap.h"
#include "CVRPSEP/MSTARSEP.H"
#include "CVRPSEP/MXF.H"
#include "CVRPSEP/NEWHTOUR.H"
#include "CVRPSEP/SORT.H"
#include "CVRPSEP/STRCOMB.H"
#include "CVRPSEP/STRNGCMP.H"
#include "CVRPSEP/TWOMATCH.H"
#include "CVRPSEP/MEMMOD.H"



typedef IloArray<IloNumVarArray> IloVar2DMatrix;  //!< An IloArray of IloNumVarArrays
typedef IloArray<IloVar2DMatrix> IloVar3DMatrix;  //!< An IloArray of IloNumVar2DArrays
typedef IloArray<IloNumArray> IloNum2DMatrix;     //!< An IloArray of IloNums
typedef IloArray<IloNum2DMatrix> IloNum3DMatrix;     //!< An IloArray of Ilo2DNums


struct CutTypes {
	int FixedChargerSingle = 0;
	int FixedChargers = 0;
	int NoChargePath = 0;
	int NoChargeSet = 0;
	int FixedEdges = 0;
	int FixedSingleEdge = 0;
	int ChargerPath = 0;
	int CapCuts = 0;
	int MultiStar = 0;
	int FCI = 0;
	int TimeInfeasibleEVPath = 0;
	int TimeInfeasibleICEVPath = 0;
	int TimeInfeasibleSet = 0;
	int Comb = 0;
	int TimeInfeasibleEVSet = 0;
};


struct testStats {
	std::string instanceName="";											//!< Name of instance
	int n=0;																//!< Number of customers
	int r=0;																//!< Number of chargers
	double EVs=0;															//!< Number of vehicles
	double ICEVs = 0;														//!< Number of vehicles
	double bestUpperBound=0;												//!< Best upper bound. Equal to optimal solution if no optimality gap is left
	double bestLowerBound=0;												//!< Best lower bound on the instance
	bool isOptimal=false;													//!< Flags whether solution is optimal
	bool isInfeasible = false;												//!< Flags whether an instance was infeasible 
	double time=0;															//!< Total time
	double CapSepTime = 0;													//!< Total time spent on capacity cuts:
	double EnumerationTime = 0;												//!< Total time spent on enumerating paths
	double TotalSepTime=0;													//!< Total time for seperation
	double TotalLoadSepData = 0;											//!< Total time spent on fetching solution
	double AddCutsTime = 0;													//!< Total time spent on adding cuts
	double LBAtRootNode=0;													//!< Lower bound before branching
	double timeRootNode=0;													//!< Time spent at the root node
	double currentUB=INFINITY;												
	double ConvertVals = 0;													//!< Time spent on converting one vector setup to matrix
	double SetupSepVectors = 0;												//!< Time spent on setting up vectors for separation routines
	CutTypes *numberOfCutsAtRootNode=new CutTypes;							//!< Number of cuts added of each type at the root node
	CutTypes *numberOfBindingCutsAfterRootNode = new CutTypes;				//!< Number of binding cuts of each type after root node
	CutTypes *numberOfBasicCutsAfterRootNode = new CutTypes;				//!< Number of basic cuts after the root node of each type
	int treeSize=0;															//!< Number of nodes search in tree
	int treeDepth=0;														//!< Maximum depth of tree
	int treeUBDepth=-1;														//!< Depth where UB was found
	int treeFeasDepth=-1;													//!< Depth where first feasible solution was found
	CutTypes *totalNumberOfCuts = new CutTypes;								//!< Total umber of cuts added of each type 
	int TSPSolved = 0;										//!< Number of TSP solved for each type
	std::vector<int> TSPNodes;												//!< Vector containing information about the number of nodes for each TSP solved
	std::vector<double> TSPTime;											//!< Vector containing information about the solution time for each TSP solved
	std::vector<std::string> InequalityType;								//!< Vector containing information about what type of inequality we try to seperate during TSP
	double avgTSPTime=0;													//!< Average time for solving TSP
	double maxTSPTime=0;													//!< Longest time spend on solving a TSP
	double minTSPtime=0;													//!< Shortest time spend on solving a TSP
	double NumberofChargers=0;												//!< Number of chargers used
	double NumberOfRoutes=0;												//!< Number of routes

};


class MFGVRP_Solver {

private:
	

	//TSP LP
	IloModel TSPModel; //TSP model
	IloCplex TSPCplex; //TSP problem solver
	IloObjective TSPObj; //Objective of TSP problem on energy
	IloObjective TSPObjTime;
	IloObjective TSPObjAdjusted; //Addjusted objective function for TSP.
	IloRangeArray TSPCons; //Constraints for the TSP problem

	//SPP Problem:
	IloModel SPPModel; //SPP model
	IloCplex SPPCplex; //SPP problem solver
	IloObjective SPPObj; //Objective of SPP problem
	IloRangeArray SPPCons; //Constraints for the SPP problem
	IloNum2DMatrix cSPP; //Distance between node i and j
	IloNum2DMatrix SPPUBs; //Upper Bounds for variables in SPP
	IloIntArray rhs;

	//Parameters
	//IloNum3DMatrix c; //Cost for driving between customer i and j with mod p
	
	
	IloNum2DMatrix u; //Distance from charging station r to node j.
	IloInt n; //Number of customers 
	IloInt E; //Number of electric vehicles
	IloInt C; //Number of conventional vehicles 
	IloInt f; //Number of charging stations
	IloInt G; //Number of green zones
	IloNumArray d; //Contains demands at all customers
	IloNumArray s; //Contains service times for all customers
	IloInt Q; //Contains load capcity
	IloNum B; //Contains battery cap for vehicles
	IloInt T; // Max tour duration
	IloNum r; //Refueling rate (time spend per unit energy refuelled)
	IloNum g; //Consumption rate (energy used per distance unit travelled)
	IloInt index; //Index variable used accross functions.

	//Cuts
	std::vector<std::vector<int>> Routes;
	std::vector<int> Route;
	IloRangeArray CapCuts; //Holds all capacity cuts
	
	//Customer sets
	std::vector<int> TSPCustomers;
	std::vector<int> GreenZoneCustomers;

	//Solution
	IloNum3DMatrix xSol1;
	IloNum2DMatrix ySol1;
	IloNum load; // Load of a route in the current solution

	//Expressions
	IloExpr expr; //New expression 1
	IloExpr expr2; //New expression 2


	//Functions:
	void BuildModels(); //Used to build Cplex models
	void Run(); //Runs the algorithm
	void SolutionChecker(); //Checks whether solution is valid
	void GetSolution(); //Gets current solution
	void PrintSolution(); //Prints the optimal solution
	void FindCapCuts(); //Identifies cap cuts
	int NonDominatedChargingPaths();

	//Vectors for indentifying violated paths:
	bool* AddedToPath; //Vector of customers currently added to path
	int* LevelNode; //Vector of nodes at each level
	int* LevelIdx; //Index of customers currently searched at a given level
	double* Slack;
	double* SlackSet;
	double* Esum;
	double* Tsum;
	double* LHSSet;
	int* minTravelTime;



public:
	struct Graph {
		int V; // No. of vertices

		// Pointer to an array containing adjacency lists
		std::list<int>* adj;
		std::list<double>* edgeValues; // added pointer to store edge values
		std::vector<int>* ChargerTails;
		std::vector<int>* ChargerHeads;
		std::vector<double>* ChargerVals;
		double* sums;
		int* Customers;
		double* MinDists;
		int Charger1_Edge;
		int Charger2_Edge;
		int Charger1;
		int Charger2;
		int Customer1;
		int Customer2;
		int Customer1_Edge;
		int Customer2_Edge;
		double edge1;
		double edge2;
		// A function used by DFS
		void DFSUtil(int v, bool visited[],double dist,double service,double Cursum);
		void DFSUtilFC(int v, bool visited[],double dist, double Cursum);
		double sumEdgesToChargers(const std::vector<int>& components, bool tail);
		double sumEdgeValues(const std::vector<int>& components);
	public:
		Graph(int V,int f); // Constructor
		~Graph();
		void addEdge(int* v, int* w, double* edgevalue, int m,int start); // adjusted function to take in edgevalue as parameter
		
		void connectedComponents();
		void connectedComponentsFC();
		std::vector<int> ConnectedComponents;
		std::vector<int> LargestViolatedSet;
		double MaxSum;
		double ExtraSum;
		double CharSum1;
		double CharSum2;
		double GreedyObj;
		double MaxGreedy;
		double ServiceTime;
		double Cap = 0;
		double ExtraDist = 0;
		bool Energy;
		bool ICEV;
		bool FixedChargers;
		bool FCFound;
		bool SetUpdated;
		bool OneWay;
		bool FixedArcs;


		CnstrMgrPointer MyCuts;
		MFGVRP_Solver* MFGVRP;
	};

	double Greedy(int* Components, int nCities, bool Energy, bool ICEV);
	double GreedyFixedChargers(int* Components, int nCities, bool OneWay, bool FixedArc);

	//Used for determining intersection:
	struct Point
	{
		int x;
		int y;
	};
	int CutCounter;
	double TSPCutOff=0;
	//Point to statistics
	testStats* Stats;

	//Cplex features
	IloEnv env; //Cplex enviroment
	IloModel Model; //Model
	IloCplex Cplex; //Cplex solver
	IloObjective Obj; //Objective function

	IloConstraintArray user;
	IloConstraintArray lazy;
	//IloRangeArray infeas;
	//IloModel infeasModel; //Model
	//IloCplex infeasCplex; //Cplex solver

	IloCplex::MIPCallbackI::NodeId NodeID;

	//Checks
	bool IsFractional; //True if solution contains fractional values
	bool HasCleanedUp; //Checks if model needs to be cleaned up before running it.
	bool Feasible; // Checks if solution is feasible
	int DrawCounter; // Counts number of iterations to properly name draw file
	double ObjValTemporary = 0;
	bool RunCuts;

	//Vectors used for generating non dominated paths:
	//Non-dominated paths:
	std::vector<std::vector<std::vector<int>>> R; // Non-dominated charging paths.
	std::vector<std::vector<int>> Rj; //Non dominated charging paths from j
	std::vector<int> Rij; // Non-dominated charging paths between a given i,j;

	
	//Public parameters:
	IloNum2DMatrix c; //Distance between node i and j
	std::vector<double> MinCharDist; //Sorted vector of customers closest to a charging station

	//Decision variables
	IloVar3DMatrix y; //1 if charging arc r is used between customer i and j
	IloVar3DMatrix x; //1 if direct arc between customer i and j is used
	IloNumVarArray xDummy;
	IloNumVarArray yDummy;
	IloVar2DMatrix xTSP; //TSP variables i,j
	IloVar2DMatrix xSPP; //SPP variables i,j
	IloNumVarArray SOC; //SOC upon arrival at j
	IloVar2DMatrix e; //Energy charged on the way to j
	IloNumVarArray time; //Arrival time at j
	std::vector<int> xCoord; //x coordinates of the nodes
	std::vector<int> yCoord; //y coordinates of the nodes
	std::vector<std::vector<std::pair <int, int> >> GZPoints; //Coordinates for green zones
	std::vector<std::pair<int, int>> CustomersPoints;

	//Constraints
	IloRangeArray cons; //Holds all initial constraints
	//Fractional Cuts
	char IntegerAndFeasible;
	int MaxNoOfCuts, NoOfEdges;
	double EpsForIntegrality, MaxViolation;

	//Energy cuts:
	//void CreateEnergyCnstrMgr(EnergyCnstrMgrPointer* ECMP, int Dim);	//USED FOR TOURNAMENT INEQUALITIES:
	//void AddEnergyCnstr(EnergyCnstrMgrPointer ECMP, int ConsType,int IntListSize, int* IntList,double RHS);
	//void ExpandEnergyCnstrMgr (EnergyCnstrMgrPointer ECMP, int NewDim);
	void TournamentInequalities(int NoOfCustomers, int Battery, int NoOfEdges, int* EdgeTail, int* EdgeHead, double* EdgeX, int* ChargingStation, CnstrMgrPointer MyCuts); //Still to be exploited how to efficiently seperate energy cuts. Start with implementing x^e cuts.

	//PARAMETERS:
	//Edge information: NoOfEdges, EdgeTail, EdgeHead, EdgeX, ChargingStation, MinCharDist = Distance to nearest charging station
	//Checks on violated cuts: slack = Cumulative slack, esum = Cumulative energy consumption, tsum = Cumulative time consumption
	//Cut sets: res = Matrix of paths in edges, resP=Matrix of paths in customers, temp = One path in edges, P = Path in actual nodes
	//Indices: index = Index of the next edge to check, LatestInclusion = The edge that was last included
	//Checks: AddEdge = True of edge shoud be added to path, AddToP = True if path should be added to matrix of paths
	//Cuts: MyCuts = Holds pointer to overall set of identified cuts
	void FindCutSets(int NoOfEdges, int* EdgeTail, int* EdgeHead, double* EdgeX, int* ChargingStation, double slack, double esum, double tsum, std::vector<std::vector<bool>>& res, std::vector<bool> temp, std::vector<std::vector<int>>& resP, std::vector<int> P, int index, int LatestInclusion, bool AddEdge, bool AddToP, CnstrMgrPointer MyCuts);
	void FindViolatedPaths(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts,bool isRoot);
	void FindViolatedChargerPaths(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts);
	void FindViolatedDepotPaths(int NoOfEdges, std::vector<int> &EdgeTail, std::vector<int>& EdgeHead, std::vector<int>& Charger, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts);
	void CheckRouteFeasibility(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<int>& Charger, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts);
	double RunTSPProblem(IloIntArray rhs, int* P, int NoOfCustomers,const std::string &CutType, int* FixedEdge = nullptr, bool AlterObjective = false, bool OneWay = false, bool time=false);
	void TwoNodeInequalities(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<int>& Charger, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts);
	void FixedChargers(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts,bool isRoot);
	void FindTimeInfeasiblePaths(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<double>& EdgeX, std::vector<int>& Origin,bool ICEV, CnstrMgrPointer MyCuts, bool isRoot);
	void AgressiveSetSearch(int NoOfEdges, int* EdgeTail, int* EdgeHead, double* EdgeX, CnstrMgrPointer MyCuts, bool Energy,bool ICEV);
	void AgressiveSetSearch_FC(int NoOfEdges, int NoOfCharEdges, int* EdgeTail, int* EdgeHead, double* EdgeX, CnstrMgrPointer MyCuts);
	void ExportCplexModel();
	void SetupTSP();

	//Public functions:
	//Constructors and deconstructors:
	MFGVRP_Solver(); //Constructor used to initialize models and parameters

	//Use later (Needs to be filled in Solver_VRP.cpp)
	~MFGVRP_Solver(); //Deconstructor used to free memory
	void LoadData(const std::string &FileName, bool RunAllCuts = false); //Used to load in problem data
	

	void CheckFeasibility(); //Checks feasibility of the tour

	int getNumCust() { return n; } //Gets number of customers
	int getChargers() { return f; } //Gets number of charging stations
	


	IloEnv getEnv() { return env; } //Gets environment of cplex model

	void setTSPCustomers(int* P, int NoOfCustomers) {TSPCustomers.clear(); for (int i = 1; i < NoOfCustomers; i++) TSPCustomers.push_back(P[i]);}
	int getNoTSPCustomers() { return TSPCustomers.size(); }
	int getTSPCutomers(int i) { return TSPCustomers[i]; }
	double getTSPObjVal() { return TSPCplex.getObjValue(); }
	double	getDemand(int idx) { return d[idx]; } // Get demand
	double getCapacity() { return Q; } // Get capacity
	double getBatteryCap() { return B; } //Get battery capacity
	double getMaxTime() { return T; }
	double getBatteryConsumption(int i, int j) { return g * c[i][j]; } //Gets energy consumption between i and j
	double getTravelTime(int i, int j) { return c[i][j]; }
	double getServiceTime(int idx) { return s[idx]+minTravelTime[idx]; }
	int getChargingStation(int i, int j, int r) { return R[i][j][r]; }
	
	int calcMinTravelTime(int idx) {
		double minVal = 999999;  
		for (int i = 1; i < n+1; i++) {
			if (i == idx) continue;
			if (minVal > c[idx][i]) minVal = c[idx][i];
		}
		return int(minVal);
	}

	//Green zone functions:
	void SetupGreenZone();
	bool onSegment(Point p, Point q, Point r);
	int orientation(Point p, Point q, Point r);
	bool doIntersect(Point p1, Point q1, Point p2, Point q2);
	
	//Return 1 if point lies inside the polytope
	int pnpoly(std::vector<std::pair<int, int>>& poly, std::pair<int, int>& testPoint)
	{
		int nvert = poly.size();
		int i, j, c = 0;
		for (i = 0, j = nvert - 1; i < nvert; j = i++) {
			if (((poly[i].second > testPoint.second) != (poly[j].second > testPoint.second)) &&
				(testPoint.first < (poly[j].first - poly[i].first) * (testPoint.second - poly[i].second) / (poly[j].second - poly[i].second) + poly[i].first))
				c = !c;
		}
		return c;
	}

	void SetupSPP(int idx);
	void RunSPP();
	void ResetSPP();
	void AddResultsToFile();
	double calc_distances(int x, int y) {  return sqrt(pow(x, 2) + pow(y, 2)); }

	void FreeMemoryFromTestStats(testStats *Stats) {
		delete Stats->numberOfCutsAtRootNode;
		delete Stats->numberOfBindingCutsAfterRootNode;
		delete Stats->numberOfBasicCutsAfterRootNode;
		delete Stats->totalNumberOfCuts;

		std::vector<int>().swap(Stats->TSPNodes);											
		std::vector<double>().swap(Stats->TSPTime);
		std::vector<std::string>().swap(Stats->InequalityType);
	};

	//Slice a vector:
	template <typename T>
	std::vector<T> slicing(std::vector<T> const& v,
		int X, int Y)
	{

		// Begin and End iterator
		auto first = v.begin() + X;
		auto last = v.begin() + Y + 1;

		// Copy the element
		std::vector<T> vector(first, last);

		// Return the results
		return vector;
	}
};




#endif