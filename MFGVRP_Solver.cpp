#include "MFGVRP_Solver.h"

//FRACTIONALITY CUTS:


ILOUSERCUTCALLBACK1(FracCuts, MFGVRP_Solver&, solver) {
	if (getMIPRelativeGap() < 0.015) return;
	double start = 0;
	double end = 0;
	double startTot = solver.Cplex.getTime();
	double StartAddingCuts = 0;
	double endTot = 0;
	double gap = getMIPRelativeGap();
	IloInt i;
	IloEnv env = solver.getEnv();
	IloInt n = solver.getNumCust()+1;
	//std::string consName;
	//std::string consName = "";
	bool edgesNeedFix = false;
	double edgeVal;
	double edgeValX;
	double edgeValY;
	bool isRoot = false;

	IloExpr expr(env);
	IloExpr expr2(env);
	int* d = new int[n];
	int* s = new int[n];
	int Q = solver.getCapacity();
	int T = solver.getMaxTime();
	//std::vector<float> EnergyComps;
	//std::vector<float>TravelTime;
	std::vector<int> CutList;
	std::vector<int> ExtList;
	bool CutsAdded = false;
	int NNodes = getNnodes();
	int NoOfChargingEdges = 0;
	if ( NNodes == 0) isRoot = true;

	int LHS = 0;
	//int* EdgeHead1 = new int[ceil(n * n / 2)];
	//int* EdgeTail1 = new int[ceil(n * n / 2)];
	//double* EdgeX1 = new double[ceil(n * n / 2)];

	std::vector<int> EdgeTail; EdgeTail.reserve(n * 2); EdgeTail.push_back(0);
	std::vector<int> EdgeHead; EdgeHead.reserve(n * 2); EdgeHead.push_back(0);
	std::vector<double> EdgeX; EdgeX.reserve(n * 2); EdgeX.push_back(-1);
	

	int TotDemand = 0;
	//One vector setup
	//std::vector<double> temp_EdgeX;
	//std::vector<int> temp_EdgeTail;
	//std::vector<int> temp_EdgeHead;

	std::vector<int> EdgeTail_ICEV; EdgeTail_ICEV.reserve(n * 2); EdgeTail_ICEV.push_back(0);
	std::vector<int> EdgeHead_ICEV; EdgeHead_ICEV.reserve(n * 2); EdgeHead_ICEV.push_back(0);
	std::vector<double> EdgeX_ICEV; EdgeX_ICEV.reserve(n * 2); EdgeX_ICEV.push_back(-1);
	std::vector<int> EdgeTail_Charger; EdgeTail_Charger.reserve(n * 2); EdgeTail_Charger.push_back(0);
	std::vector<int> EdgeHead_Charger; EdgeHead_Charger.reserve(n * 2); EdgeHead_Charger.push_back(0);
	std::vector<double> EdgeX_Charger; EdgeX_Charger.reserve(n * 2); EdgeX_Charger.push_back(-1);
	std::vector<int> EdgeTail_Direct; EdgeTail_Direct.reserve(n * 2); EdgeTail_Direct.push_back(0);
	std::vector<int> EdgeHead_Direct; EdgeHead_Direct.reserve(n * 2); EdgeHead_Direct.push_back(0);
	std::vector<double> EdgeX_Direct; EdgeX_Direct.reserve(n * 2); EdgeX_Direct.push_back(-1);
	std::vector<int> Origin; Origin.reserve(n * 2); Origin.push_back(0);
	//For FCI:
	int* Label = new int[n + 1];
	int MaxIdx;
	int MinIdx;
	int k;


	//For Multistar:
	int* NList = new int[n + 1];
	int* TList = new int[n + 1];
	int* CList = new int[n + 1];
	int intA;
	int intB;
	int intL;

	//For Comb:
	int NoOfTeeth = 0;

	int OriginID = 1;
	int Qmin = 0;


	std::vector<int> reduced_EdgeTail;
	std::vector<int> reduced_EdgeHead;
	std::vector<double> reduced_EdgeX;

	//CVRPSEP parameters:
	CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;
	CMGR_CreateCMgr(&MyCutsCMP, 100);
	CMGR_CreateCMgr(&MyOldCutsCMP, 100);
	


	IloNum3DMatrix xTemp(env, n);
	IloNum3DMatrix yTemp(env, n);


	IloNumArray xTemp1d(env, solver.xDummy.getSize());
	IloNumArray yTemp1d(env, solver.yDummy.getSize());

	start= solver.Cplex.getTime();
	getValues(xTemp1d, solver.xDummy);
	getValues(yTemp1d, solver.yDummy);
	end = solver.Cplex.getTime();
	solver.Stats->TotalLoadSepData += end - start;
	start = solver.Cplex.getTime();
	int cntX = 0;
	int cntY = 0;
	for (int i = 0; i < n; i++)
	{
		xTemp[i] = IloNum2DMatrix(env, n);
		yTemp[i] = IloNum2DMatrix(env, n);
		for (int j = 0; j < n; j++)
		{
			xTemp[i][j] = IloNumArray(env, 2);
			for (int k = 0; k < 2; k++)
			{
				xTemp[i][j][k] = xTemp1d[cntX];
				cntX++;
			}
			yTemp[i][j] = IloNumArray(env, solver.R[i][j].size());
			for (int r = 0; r < solver.R[i][j].size(); r++)
			{
				yTemp[i][j][r] = yTemp1d[cntY];
				cntY++;
			}
		}
	}

	xTemp1d.end();
	yTemp1d.end();
	
	
	end = solver.Cplex.getTime();
	solver.Stats->ConvertVals += end - start;
	start = solver.Cplex.getTime();
	for (IloInt i = 0; i < n; i++)
	{
		if (i > 0) { d[i] = solver.getDemand(i); TotDemand += d[i]; s[i] = solver.getServiceTime(i); }
		Qmin += d[i];
		for (IloInt j = i + 1; j < n; j++)
		{
			edgeVal = 0;
			expr.clear();
		
			for (IloInt r = 0; r < solver.R[i][j].size(); r++)
			{
				edgeValY = yTemp[i][j][r] + yTemp[j][i][r];
				edgeVal += edgeValY;
				expr += solver.y[i][j][r] + solver.y[j][i][r];
				if (edgeValY > 0.001)
				{
					//One vector setup:
					if (yTemp[i][j][r] > 0.001)
					{
						//EdgeTail_Charger[ChargingArcs]=i;
						//EdgeHead_Charger[ChargingArcs]=solver.R[i][j][r];
						//EdgeX_Charger[ChargingArcs]=yTemp[i][j][r];
						//ChargingArcs++;
						//EdgeTail_Charger[ChargingArcs]=solver.R[i][j][r];
						//EdgeHead_Charger[ChargingArcs] = j;
						//EdgeX_Charger[ChargingArcs] = yTemp[i][j][r];
						//ChargingArcs++;
						EdgeTail_Charger.push_back(i);
						EdgeHead_Charger.push_back(solver.R[i][j][r]);
						EdgeX_Charger.push_back(yTemp[i][j][r]);
						Origin.push_back(OriginID);
						EdgeTail_Charger.push_back(solver.R[i][j][r]);
						EdgeHead_Charger.push_back(j);
						EdgeX_Charger.push_back(yTemp[i][j][r]);
						Origin.push_back(OriginID);
						OriginID++;
					}
					if (yTemp[j][i][r] > 0.001)
					{

						//EdgeTail_Charger[ChargingArcs] = j;
						//EdgeHead_Charger[ChargingArcs] = solver.R[j][i][r];
						//EdgeX_Charger[ChargingArcs] = yTemp[j][i][r];
						//ChargingArcs++;
						//EdgeTail_Charger[ChargingArcs] = solver.R[j][i][r];
						//EdgeHead_Charger[ChargingArcs] = i;
						//EdgeX_Charger[ChargingArcs] = yTemp[j][i][r];
						//ChargingArcs++;
						EdgeTail_Charger.push_back(j);
						EdgeHead_Charger.push_back(solver.R[j][i][r]);
						EdgeX_Charger.push_back(yTemp[j][i][r]);
						Origin.push_back(OriginID);
						EdgeTail_Charger.push_back(solver.R[j][i][r]);
						EdgeHead_Charger.push_back(i);
						EdgeX_Charger.push_back(yTemp[j][i][r]);
						Origin.push_back(OriginID);
						OriginID++;
					}
				}
			}
			edgeValX = xTemp[i][j][1] + xTemp[j][i][1] + xTemp[i][j][0] + xTemp[j][i][0];
			edgeVal +=edgeValX;
			expr += solver.x[i][j][0] + solver.x[j][i][0] + solver.x[i][j][1] + solver.x[j][i][1];


			if (edgeValX > 0.001) {
				//One vector setup
				if (xTemp[i][j][0] > 0.001) {
					if (i == 0) {

						EdgeTail_Charger.push_back(i);
						EdgeHead_Charger.push_back(j);
						EdgeX_Charger.push_back(xTemp[i][j][0]);
						Origin.push_back(i);
						//EdgeTail_Charger[ChargingArcs] = i;
						//EdgeHead_Charger[ChargingArcs] =j;
						//EdgeX_Charger[ChargingArcs] = xTemp[i][j][0];
						//ChargingArcs++;


					}
					else {
						EdgeTail_Direct.push_back(i);
						EdgeHead_Direct.push_back(j);
						EdgeX_Direct.push_back(xTemp[i][j][0]);

						//EdgeTail_Direct[DirectArcs] = i;
						//EdgeHead_Direct[DirectArcs] = j;
						//EdgeX_Direct[DirectArcs] = xTemp[i][j][0];
						//DirectArcs++;
					}

				}
				if (xTemp[j][i][0] > 0.001) {
					if (i == 0) {
						//EdgeTail_Charger[ChargingArcs] = j;
						//EdgeHead_Charger[ChargingArcs] = i;
						//EdgeX_Charger[ChargingArcs] = xTemp[j][i][0];
						//ChargingArcs++;

						EdgeTail_Charger.push_back(j);
						EdgeHead_Charger.push_back(i);
						EdgeX_Charger.push_back(xTemp[j][i][0]);
						Origin.push_back(i);
					}
					else {
						//EdgeTail_Direct[DirectArcs] = j;
						//EdgeHead_Direct[DirectArcs] = i;
						//EdgeX_Direct[DirectArcs] = xTemp[j][i][0];
						//DirectArcs++;

						EdgeTail_Direct.push_back(j);
						EdgeHead_Direct.push_back(i);
						EdgeX_Direct.push_back(xTemp[j][i][0]);
					}


				}
				if (xTemp[i][j][1] > 0.001) {

					EdgeTail_ICEV.push_back(i);
					EdgeHead_ICEV.push_back(j);
					EdgeX_ICEV.push_back(xTemp[i][j][1]);
					//EdgeTail_ICEV[ICEVArcs]=i;
					//EdgeHead_ICEV[ICEVArcs] = j;
					//EdgeX_ICEV[ICEVArcs] = xTemp[i][j][1];
					//ICEVArcs++;
					
				}
				if (xTemp[j][i][1] > 0.001) {
					EdgeTail_ICEV.push_back(j);
					EdgeHead_ICEV.push_back(i);
					EdgeX_ICEV.push_back(xTemp[j][i][1]);
					//EdgeTail_ICEV[ICEVArcs] = j;
					//EdgeHead_ICEV[ICEVArcs] = i;
					//EdgeX_ICEV[ICEVArcs] = xTemp[j][i][1];
					//ICEVArcs++;
				}
			}
			if (edgeVal > 0)
			{
				EdgeTail.push_back(i == 0 ? n : i);
				EdgeHead.push_back(j);
				EdgeX.push_back(edgeVal);
				if (edgeVal > 1.001 && i!=0)
				{
					edgesNeedFix = true;
					add(expr <= 1);
				}
			}
		}
	}
	Qmin -= (n - 2) * Q;
	NoOfChargingEdges = EdgeTail_Charger.size();
	//temp_EdgeTail.insert(temp_EdgeTail.begin(), 0);
	//temp_EdgeHead.insert(temp_EdgeHead.begin(), 0);
	//temp_EdgeX.insert(temp_EdgeX.begin(), 0);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			xTemp[i][j].end();
			yTemp[i][j].end();
		}
		xTemp[i].end();
		yTemp[i].end();
	}
	xTemp.end();
	yTemp.end();

	end = solver.Cplex.getTime();
	solver.Stats->SetupSepVectors += end - start;





	//FOR GETTING THE ROOT NODE:
	//if (getNnodes()==1)
	//{

	//}


	//printf("\nDepot arcs: %d\tDirect arcs: %d\tCharging arcs: %d\tICEV arcs:%d\n\n", DepotArcs, DirectArcs, ChargingArcs,ICEVArcs);
	//for (int i = 0; i < temp_EdgeHead.size(); i++)
	//{
	//	printf("%d\t%d\t%f\n", temp_EdgeTail[i], temp_EdgeHead[i],temp_EdgeX[i]);
	//}
	//printf("Stop!");
	//printf("\nAll arcs:");
	////for (int i = 0; i < EdgeHead.size(); i++)
	////{
	////	printf("\n%d\t%d\t%f", EdgeTail[i], EdgeHead[i], EdgeX[i]);
	////}

	//std::vector<int> p1 = solver.slicing(temp_EdgeTail, DepotArcs, DepotArcs+DirectArcs-1);
	//std::vector<int> p2 = solver.slicing(temp_EdgeHead, DepotArcs, DepotArcs + DirectArcs-1);
	//std::vector<double> p3 = solver.slicing(temp_EdgeX, DepotArcs, DepotArcs + DirectArcs-1);
	//printf("\n\nReduced vector\n");

	//for (int i = 0; i < p1.size(); i++)
	//{
	//	printf("%d\t%d\t%f\n", p1[i], p2[i], p3[i]);
	//}

	//Input for Graph:
	//std::cout << "\n\nMy draw values:\n";
	//for (IloInt i = 0; i < EdgeTailDraw.size(); i++)
	//{
	//	std::cout << "\n" << EdgeTailDraw[i] << "\t" << EdgeHeadDraw[i] << "\t" << EdgeDraw[i] << "\t" << EnergyComps[i] << "\t" << TravelTime[i];
	//}
	//std::cout << "\n\nMy draw values:\n";
	//std::cout << "\n\nMy load cuts values:\n";
	//for (IloInt i = 1; i <= solver.NoOfEdges; i++)
	//{
	//	std::cout << "\n" << EdgeTail[i] << " " << EdgeHead[i] << " " << EdgeX[i];
	//}
	

	//std::cout << "\n\nMy energy cuts values:\n";
	//for (IloInt i = 0; i < EdgeTailY.size(); i++)
	//{
	//	std::cout << "\n" << EdgeTailY[i] << " " << EdgeHeadY[i] << " " << EdgeY[i] << " " << EdgeChargStation[i];
	//}
	//printf("\n");
	//std::cout << "\n\nMy draw values:\n";
	//for (IloInt i = 0; i < EdgeTailDraw.size(); i++)
	//{
	//	std::cout << "\n" << EdgeTailDraw[i] << "\t" << EdgeHeadDraw[i] << "\t" << EdgeDraw[i];
	//}
	//std::cout << "\n\nx_ICEV:\n";
	//for (IloInt i = 0; i < n; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n; j++)
	//	{
	//		std::cout << getValue(solver.x[i][j][1])<<" ";
	//	}
	//}

	//std::cout << "\n\nx_EV:\n";
	//for (IloInt i = 0; i < n; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n; j++)
	//	{
	//		std::cout << getValue(solver.x[i][j][0]) << " ";
	//	}
	//}

	//std::cout << "\n\ny_EV:\n";
	//for (IloInt i = 0; i < n; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n; j++)
	//	{
	//		std::cout << "{";
	//		for (IloInt r = 0; r < solver.R[i][j].size(); r++)
	//		{
	//			std::cout << getValue(solver.y[i][j][r]) << ";";
	//		}
	//		std::cout << "} ";
	//	}
	//}

	////For multistar:
	//std::vector<int> NList;
	//std::vector<int> TList;
	//std::vector<int> CList;
	//double sigma = 0;
	//double lambda = 0;
	if (!edgesNeedFix) //If no edges needs to be fixed, we can initiate our seperation problem:
	{	
		start = solver.Cplex.getTime();
		int cntCut = MyCutsCMP->Size;
		//if (TotDemand > Q)
		//{
		
			CAPSEP_SeparateCapCuts(n - 1,
				d,
				Q,
				EdgeTail.size()-1,
				EdgeTail.data(),
				EdgeHead.data(),
				EdgeX.data(),
				MyOldCutsCMP,
				solver.MaxNoOfCuts,
				solver.EpsForIntegrality,
				&solver.IntegerAndFeasible,
				&solver.MaxViolation,
				MyCutsCMP);
		//}
			if (NNodes==0) {
				//FCISEP_SeparateFCIs(n - 1,
				//	d,
				//	Q,
				//	solver.NoOfEdges,
				//	EdgeTail.data(),
				//	EdgeHead.data(),
				//	EdgeX.data(),
				//	MyOldCutsCMP,
				//	20000,
				//	10,
				//	&solver.MaxViolation,
				//	MyCutsCMP
				//);
				if (getNiterations()%2==0)
				{
					if (solver.MaxViolation < 0.05)
					{
						////NEEDS TO BE CHANGED TO MULTISTAR CUTS
						MSTARSEP_SeparateMultiStarCuts(n - 1,
							d,
							Q,
							EdgeTail.size() - 1,
							EdgeTail.data(),
							EdgeHead.data(),
							EdgeX.data(),
							MyOldCutsCMP,
							solver.MaxNoOfCuts * 2 - MyCutsCMP->Size,
							&solver.MaxViolation,
							MyCutsCMP);
					}
					if (solver.MaxViolation < 0.1)
					{
						COMBSEP_SeparateCombs(
							n - 1,
							d,
							Q,
							Qmin,
							EdgeTail.size() - 1,
							EdgeTail.data(),
							EdgeHead.data(),
							EdgeX.data(),
							10,
							&solver.MaxViolation,
							MyCutsCMP
						);
					}
				}
				else
				{
					if (solver.MaxViolation < 0.1)
					{
						COMBSEP_SeparateCombs(
							n - 1,
							d,
							Q,
							Qmin,
							EdgeTail.size() - 1,
							EdgeTail.data(),
							EdgeHead.data(),
							EdgeX.data(),
							10,
							&solver.MaxViolation,
							MyCutsCMP
						);
					}
					if (solver.MaxViolation < 0.05)
					{
						MSTARSEP_SeparateMultiStarCuts(n - 1,
							d,
							Q,
							EdgeTail.size() - 1,
							EdgeTail.data(),
							EdgeHead.data(),
							EdgeX.data(),
							MyOldCutsCMP,
							solver.MaxNoOfCuts * 2 - MyCutsCMP->Size,
							&solver.MaxViolation,
							MyCutsCMP);
					}
				}
				
			}
		
		end = solver.Cplex.getTime();
		solver.Stats->CapSepTime += end - start;
		//}else {
		//	//FCISEP_SeparateFCIs(n - 1,
		//	//	d,
		//	//	Q,
		//	//	solver.NoOfEdges,
		//	//	EdgeTail.data(),
		//	//	EdgeHead.data(),
		//	//	EdgeX.data(),
		//	//	MyOldCutsCMP,
		//	//	100,
		//	//	10,
		//	//	&solver.MaxViolation,
		//	//	MyCutsCMP
		//	//);

		//}


		if (MyCutsCMP->Size == 0) {
			CAPSEP_SeparateCapCuts(n - 1,
				s,
				T,
				EdgeTail.size() - 1,
				EdgeTail.data(),
				EdgeHead.data(),
				EdgeX.data(),
				MyOldCutsCMP,
				solver.MaxNoOfCuts,
				solver.EpsForIntegrality,
				&solver.IntegerAndFeasible,
				&solver.MaxViolation,
				MyCutsCMP);
		}
		
		//printf("\nUser Cap");

		/*solver.DrawCounter++;*/
		/*std::string strFilename = "Grafer/Elbils_løsning_n21_f9.tex";
		char* FileName;
		FileName = &strFilename[0];
		char GraphName[] = "Elbils loesning (sparse)";
		DRAWGRAPH_DrawCVRPSparseGraph(n - 1, f, EdgeTailDraw.data(), 0, solver.xCoord.data(), solver.yCoord.data(), EdgeDraw.size(), EdgeTailDraw.data(), EdgeHeadDraw.data(), EnergyComps.data(), TravelTime.data(), d, s, EdgeDraw.data(), 0.7, GraphName, FileName, false);*/
		



		if (MyCutsCMP->Size==0) {
			//if (NNodes>0)
			//{
				
				start = solver.Cplex.getTime();

				solver.CutCounter = 0;

				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);

				//for (int i = DepotArcs + 1; i <= DepotArcs + DirectArcs; i++)
				//{
				//	//reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	//reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	//reduced_EdgeX.push_back(temp_EdgeX[i]);
				//	printf("\n%d\t%d\t%f", temp_EdgeTail[i], temp_EdgeHead[i], temp_EdgeX[i]);
				//}
				//
				//printf("\nNew");
				//for (int i = 1; i <EdgeHead_Direct.size(); i++)
				//{
				//	printf("\n%d\t%d\t%f", EdgeTail_Direct[i], EdgeHead_Direct[i], EdgeX_Direct[i]);
				//}

				//printf("\nNo charge");
				cntCut = MyCutsCMP->Size;
				////printf("\n start violated paths");
				//printf("\We stated no-charge user");
				//solver.FindViolatedPaths(reduced_EdgeTail.size() - 1, reduced_EdgeTail, reduced_EdgeHead, reduced_EdgeX, MyCutsCMP, isRoot);
				solver.FindViolatedPaths(EdgeTail_Direct.size() - 1, EdgeTail_Direct, EdgeHead_Direct, EdgeX_Direct, MyCutsCMP, isRoot);
				//printf("\nUser No charge");
				if (cntCut != MyCutsCMP->Size) solver.NodeID = getNodeId();
				//printf("\We finished no-charge user");

				//printf("\n end violated paths");
				if (MyCutsCMP->Size == 0)
				{
					cntCut = MyCutsCMP->Size;
					/*temp_EdgeTail.insert(temp_EdgeTail.begin(), 0);
					temp_EdgeHead.insert(temp_EdgeHead.begin(), 0);
					temp_EdgeX.insert(temp_EdgeX.begin(), 0);*/
					//printf("\n start violated Charger Paths");
					//printf("\nWe started Fixed user");



					//solver.FixedChargers(temp_EdgeTail.size() - ICEVArcs, temp_EdgeTail, temp_EdgeHead, temp_EdgeX, MyCutsCMP, isRoot);
					//printf("\nfixed");
					//EdgeTail_Charger.reserve(EdgeTail_Charger.size()+ EdgeTail_Direct.size() - 1);
					EdgeTail_Charger.insert(EdgeTail_Charger.end(), EdgeTail_Direct.begin()+1, EdgeTail_Direct.end());
					//EdgeHead_Charger.reserve(EdgeHead_Charger.size() + EdgeHead_Direct.size()-1);
					EdgeHead_Charger.insert(EdgeHead_Charger.end(), EdgeHead_Direct.begin()+1, EdgeHead_Direct.end());
					//EdgeX_Charger.reserve( EdgeX_Charger.size() + EdgeX_Direct.size()-1);
					EdgeX_Charger.insert(EdgeX_Charger.end(), EdgeX_Direct.begin()+1, EdgeX_Direct.end());




					solver.FixedChargers(EdgeTail_Charger.size(), EdgeTail_Charger, EdgeHead_Charger, EdgeX_Charger, MyCutsCMP, isRoot);
					//printf("\nUser Fixed");
					//printf("\n end violated Charger Paths");
				}
			//}
			
			//Check EVs:
			if (MyCutsCMP->Size == 0)
			{


				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);

				//for (int i = DepotArcs + 1; i <= DepotArcs + DirectArcs + ChargingArcs; i++)
				//{
				//	//reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	//reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	//reduced_EdgeX.push_back(temp_EdgeX[i]);
				//	printf("\n%d\t%d\t%f", temp_EdgeTail[i], temp_EdgeHead[i], temp_EdgeX[i]);
				//}


				//printf("\nNew");
				//for (int i = 1; i <EdgeHead_Charger.size(); i++)
				//{
				//	printf("\n%d\t%d\t%f", EdgeTail_Charger[i], EdgeHead_Charger[i], EdgeX_Charger[i]);
				//}
				//printf("\nOrigins");
				//for (int i = 1; i < Origin.size(); i++)
				//{
				//	printf("\n%d", Origin[i]);
				//}
				//printf("\ntime ev");
				//cntCut = MyCutsCMP->Size;
				//solver.FindTimeInfeasiblePaths(reduced_EdgeTail.size() - 1, reduced_EdgeTail, reduced_EdgeHead, reduced_EdgeX, false, MyCutsCMP,isRoot);
				solver.FindTimeInfeasiblePaths(EdgeTail_Charger.size() - 1, EdgeTail_Charger, EdgeHead_Charger, EdgeX_Charger,Origin, false, MyCutsCMP, isRoot);
				//printf("\nend");
				if (cntCut != MyCutsCMP->Size) solver.NodeID = getNodeId();
				//printf("\nUser Time paths");
			}
			//Check ICEVS
			if( MyCutsCMP->Size == 0)
			{
				//cntCut = MyCutsCMP->Size;
				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);

				//for (int i = DepotArcs + DirectArcs + ChargingArcs + 1; i < temp_EdgeTail.size(); i++)
				//{
				//	reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	reduced_EdgeX.push_back(temp_EdgeX[i]);
				//}
				cntCut = MyCutsCMP->Size;

				//solver.FindTimeInfeasiblePaths(reduced_EdgeTail.size() - 1, reduced_EdgeTail, reduced_EdgeHead, reduced_EdgeX, true, MyCutsCMP,isRoot);
				//printf("\ntime icev");

				solver.FindTimeInfeasiblePaths(EdgeTail_ICEV.size() - 1, EdgeTail_ICEV, EdgeHead_ICEV, EdgeX_ICEV,Origin, true, MyCutsCMP, isRoot);
				solver.NodeID = getNodeId();
				//printf("\nUser ICEV");
			}
			if (NNodes == 0 && MyCutsCMP->Size == 0)
			{

				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);
				//cntCut = MyCutsCMP->Size;
				//for (int i = DepotArcs + 1; i <= DepotArcs + DirectArcs; i++)
				//{
				//	reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	reduced_EdgeX.push_back(temp_EdgeX[i]);
				//	printf("\n%d\t%d\t%f", temp_EdgeTail[i], temp_EdgeHead[i], temp_EdgeX[i]);
				//}


				//printf("\nNew");
				//for (int i = 1; i <EdgeHead_Direct.size(); i++)
				//{
				//	printf("\n%d\t%d\t%f", EdgeTail_Direct[i], EdgeHead_Direct[i], EdgeX_Direct[i]);
				//}


				//solver.AgressiveSetSearch(reduced_EdgeTail.size(), reduced_EdgeTail.data(), reduced_EdgeHead.data(), reduced_EdgeX.data(), MyCutsCMP, true, false);
				//printf("\nagg no charge");
				solver.AgressiveSetSearch(EdgeTail_Direct.size(), EdgeTail_Direct.data(), EdgeHead_Direct.data(), EdgeX_Direct.data(), MyCutsCMP, true, false);
				//printf("\nUser agg charge");
			}


			if (NNodes == 0 && MyCutsCMP->Size == 0)
			{
				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);

				//for (int i = DepotArcs + 1; i <= DepotArcs + DirectArcs + ChargingArcs; i++)
				//{
				//	reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	reduced_EdgeX.push_back(temp_EdgeX[i]);
				//}
				//printf("\nagg time ev");

				//cntCut = MyCutsCMP->Size;
				//solver.AgressiveSetSearch(reduced_EdgeTail.size(), reduced_EdgeTail.data(), reduced_EdgeHead.data(), reduced_EdgeX.data(), MyCutsCMP, false, false);
				solver.AgressiveSetSearch(EdgeTail_Charger.size(), EdgeTail_Charger.data(), EdgeHead_Charger.data(), EdgeX_Charger.data(), MyCutsCMP, false, false);
				//if (cntCut != MyCutsCMP->Size) solver.NodeID = getNodeId();
				//printf("\nUser agg Time EV");
			}
			if (NNodes == 0 && MyCutsCMP->Size == 0)
			{

				//solver.NodeID = getNodeId();
				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);

				//for (int i = DepotArcs + DirectArcs + ChargingArcs + ICEVDepotArcs + 1; i < temp_EdgeTail.size(); i++)
				//{
				//	reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	reduced_EdgeX.push_back(temp_EdgeX[i]);
				//	/*printf("\n%d\t%d\t%f", temp_EdgeTail[i], temp_EdgeHead[i], temp_EdgeX[i]);*/
				//}
				//printf("\nagg time ICEV");
				//solver.AgressiveSetSearch(reduced_EdgeTail.size(), reduced_EdgeTail.data(), reduced_EdgeHead.data(), reduced_EdgeX.data(), MyCutsCMP, false, true);
				solver.AgressiveSetSearch(EdgeTail_ICEV.size(), EdgeTail_ICEV.data(), EdgeHead_ICEV.data(), EdgeX_ICEV.data(), MyCutsCMP, false, true);
				//printf("\nUser agg ICEV");
			}
			//if (NNodes == 0 && MyCutsCMP->Size == 0)
			//{
			//	//printf("\nStart User Agg FIX");
			//	solver.AgressiveSetSearch_FC(EdgeTail_Charger.size(), NoOfChargingEdges, EdgeTail_Charger.data(), EdgeHead_Charger.data(), EdgeX_Charger.data(),MyCutsCMP);
			//	printf("\nUser Agg FIX");
			//}
			solver.NodeID = getNodeId();
			
			end = solver.Cplex.getTime();
			solver.Stats->EnumerationTime += end-start;
		}//, EnergyCuts);

		
		//std::cout << "\n\nNumber of violated cuts: " << MyCutsCMP->Size << "\n\nCuts are:";
		if (MyCutsCMP->Size != 0)
		{
			StartAddingCuts = solver.Cplex.getTime();
			CutsAdded = true;
			for (IloInt cons = 0; cons < MyCutsCMP->Size; cons++)
			{
				//consName = "";
				//	std::cout << "\n";
				//	for (IloInt j = 1; j <= MyCutsCMP->CPL[c]->IntListSize; j++)
				//	{

				//		std::cout << MyCutsCMP->CPL[c]->IntList[j] << " ";
				//	}

				//	std::cout << " with RHS = " << MyCutsCMP->CPL[c]->RHS<<"\n\nList is: ";

				CutList.clear();
				ExtList.clear();
				if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_CAP)
				{

					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->CapCuts++;
						solver.Stats->numberOfCutsAtRootNode->CapCuts++;
					}
					else solver.Stats->totalNumberOfCuts->CapCuts++;
					
					//printf("\n");
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j] == n ? 0 : MyCutsCMP->CPL[cons]->IntList[j]);
						//printf("%d ", CutList.back());
					}
					
					expr.clear();
					for (IloInt i = 0; i < CutList.size(); i++)
					{
						//std::cout << CutList[i]<<" ";
					}


					for (int i = 0; i < CutList.size(); i++)
					{
						for (int j = 0 ; j < CutList.size(); j++)
						{
							expr += solver.x[CutList[i]][CutList[j]][0] + solver.x[CutList[i]][CutList[j]][1];
							for (int r = 0; r < solver.R[CutList[i]][CutList[j]].size(); r++)
							{
								expr += solver.y[CutList[i]][CutList[j]][r];
							}
						}
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//consName = "Cap_Cut";
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_FCI)
				{

					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->FCI++;
						solver.Stats->numberOfCutsAtRootNode->FCI++;
					}
					else solver.Stats->totalNumberOfCuts->FCI++;
					
					for (int j = 0; j < n; j++) Label[j] = 0;

					MaxIdx = 0;
					for (int SubsetNr = 1;
						SubsetNr <= MyCutsCMP->CPL[cons]->ExtListSize;
						SubsetNr++)
					{

						MinIdx = MaxIdx + 1;
						MaxIdx = MinIdx + MyCutsCMP->CPL[cons]->ExtList[SubsetNr] - 1;

						for (int j = MinIdx; j <= MaxIdx; j++)
						{
							k = MyCutsCMP->CPL[cons]->IntList[j];
							Label[k] = SubsetNr;
						}
					}

					expr.clear();
					for (int i = 0; i < n; i++)
					{
						for (int j = i + 1; j < n; j++)
						{
							if (Label[i] == 0 && Label[j] != 0)
							{
								expr += solver.x[i][j][0] + solver.x[j][i][0] + solver.x[i][j][1] + solver.x[j][i][1];
								for (int r = 0; r < solver.R[i][j].size(); r++) expr += solver.y[i][j][r] + solver.y[j][i][r];
							}
							else if (Label[i] != 0 && Label[j] == 0) {
								expr += solver.x[i][j][0] + solver.x[j][i][0] + solver.x[i][j][1] + solver.x[j][i][1];
								for (int r = 0; r < solver.R[i][j].size(); r++) expr += solver.y[i][j][r] + solver.y[j][i][r];
							}
						}
					}
					for (int p = 1; p <= MyCutsCMP->CPL[cons]->ExtListSize; p++)
					{
						for (int i = 0; i < n; i++)
						{
							for (int j = i + 1; j < n; j++)
							{
								if (Label[i] != p && Label[j] == p)
								{
									expr += solver.x[i][j][0] + solver.x[j][i][0] + solver.x[i][j][1] + solver.x[j][i][1];
									for (int r = 0; r < solver.R[i][j].size(); r++) expr += solver.y[i][j][r] + solver.y[j][i][r];
								}
								else if (Label[i] == p && Label[j] != p)
								{
									expr += solver.x[i][j][0] + solver.x[j][i][0] + solver.x[i][j][1] + solver.x[j][i][1];
									for (int r = 0; r < solver.R[i][j].size(); r++) expr += solver.y[i][j][r] + solver.y[j][i][r];
								}
							}
						}
					}

					add(expr >= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr >= MyCutsCMP->CPL[cons]->RHS);
					//consName = "FCI_Cut";
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_STR_COMB)
				{

					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->Comb++;
						solver.Stats->numberOfCutsAtRootNode->Comb++;
					}
					else solver.Stats->totalNumberOfCuts->Comb++;

					// Book memory:
					NoOfTeeth = MyCutsCMP->CPL[cons]->Key;
					int j = 0;
					char** InTooth = new char* [n + 1];
					for (int i = 0; i < n + 1; i++) InTooth[i] = new char[NoOfTeeth + 1];

					for (int Node = 0; Node <= n; Node++)
						for (int Tooth = 0; Tooth <= NoOfTeeth; Tooth++)
							InTooth[Node][Tooth] = 0;


					// Read cut:
					for (int k = 1; k <= MyCutsCMP->CPL[cons]->IntListSize; k++)
					{
						j = MyCutsCMP->CPL[cons]->IntList[k];
						if (j == n) j = 0;
						InTooth[j][0] = 1; /* Node j is in the handle */
					}
					for (int t = 1; t <= NoOfTeeth; t++)
					{
						//printf("\n");
						MinIdx = MyCutsCMP->CPL[cons]->ExtList[t];
						if (t == NoOfTeeth)
							MaxIdx = MyCutsCMP->CPL[cons]->ExtListSize;
						else
							MaxIdx = MyCutsCMP->CPL[cons]->ExtList[t + 1] - 1;
						for (int k = MinIdx; k <= MaxIdx; k++)
						{
							j = MyCutsCMP->CPL[cons]->ExtList[k];
							if (j == n) j = 0;
							InTooth[j][t] = 1; /* Node j is in tooth t */
							//printf("%d", InTooth[j][t]);
						}
					}

					//Add cut:
					expr.clear();
					for (int i = 0; i < n; i++)
					{
						for (int j = i + 1; j < n; j++)
						{
							if (InTooth[i][0] != InTooth[j][0]) {
								expr += solver.x[i][j][0] + solver.x[j][i][0] + solver.x[i][j][1] + solver.x[j][i][1];
								for (int r = 0; r < solver.R[i][j].size(); r++) expr += solver.y[i][j][r] + solver.y[j][i][r];
							}
							for (int t = 1; t <= NoOfTeeth; t++)
							{
								if (InTooth[i][t] != InTooth[j][t]) {
									expr += solver.x[i][j][0] + solver.x[j][i][0] + solver.x[i][j][1] + solver.x[j][i][1];
									for (int r = 0; r < solver.R[i][j].size(); r++) expr += solver.y[i][j][r] + solver.y[j][i][r];
								}
							}
						}

					}

					add(expr >= MyCutsCMP->CPL[cons]->RHS);

					//solver.infeas.add(expr >= MyCutsCMP->CPL[cons]->RHS);
					//consName = "COMB_Cut";
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
					//Release memomry:
					for (int i = 0; i < n + 1; i++) delete[] InTooth[i];
					delete[] InTooth;

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_MSTAR) {
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->MultiStar++;
						solver.Stats->numberOfCutsAtRootNode->MultiStar++;
					}
					else solver.Stats->totalNumberOfCuts->MultiStar++;

					/* Nucleus: */
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; j++)
						NList[j] = MyCutsCMP->CPL[cons]->IntList[j];
					/* Satellites: */
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->ExtListSize; j++)
						TList[j] = MyCutsCMP->CPL[cons]->ExtList[j];
					/* Connectors: */
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->CListSize; j++)
						CList[j] = MyCutsCMP->CPL[cons]->CList[j];

					/* Coefficients of the cut: */
					intA = MyCutsCMP->CPL[cons]->A;
					intB = MyCutsCMP->CPL[cons]->B;
					intL = MyCutsCMP->CPL[cons]->L;
					/* Lambda=L/B, Sigma=A/B */

					/*Add the cut to the LP*/
					expr.clear();
					for (int i = 1; i <= MyCutsCMP->CPL[cons]->IntListSize; i++)
					{
						for (int j = i + 1; j <= MyCutsCMP->CPL[cons]->IntListSize; j++)
						{
							expr += 2 * intB * (solver.x[NList[i]][NList[j]][0] + solver.x[NList[j]][NList[i]][0] + solver.x[NList[i]][NList[j]][1] + solver.x[NList[j]][NList[i]][1]);
							for (int r = 0; r < solver.R[NList[i]][NList[j]].size(); r++) expr += 2 * intB * (solver.y[NList[i]][NList[j]][r] + solver.y[NList[j]][NList[i]][r]);
						}
					}
					for (int i = 1; i <= MyCutsCMP->CPL[cons]->CListSize; i++)
					{
						for (int j = 1; j <= MyCutsCMP->CPL[cons]->ExtListSize; j++)
						{
							expr += intA * (solver.x[CList[i]][TList[j]][0] + solver.x[TList[j]][CList[i]][0] + solver.x[CList[i]][TList[j]][1] + solver.x[TList[j]][CList[i]][1]);
							for (int r = 0; r < solver.R[CList[i]][TList[j]].size(); r++) expr += intA * (solver.y[CList[i]][TList[j]][r] + solver.y[TList[j]][CList[i]][r]);
						}
					}
					add(expr <= 2 * intB * MyCutsCMP->CPL[cons]->IntListSize - intL);
					//solver.infeas.add(expr <= 2 * intB * MyCutsCMP->CPL[cons]->IntListSize - intL);
					//consName = "MTSTAR_Cut";
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_PATH)
				{

					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->NoChargePath++;
						solver.Stats->numberOfCutsAtRootNode->NoChargePath++;
					}
					else solver.Stats->totalNumberOfCuts->NoChargePath++;
					//std::cout << "\nNCP: ";
					//consName = "NCP -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
						//std::cout << CutList[j-1] << " ";
					}
					//std::cout << "\nPath: ";
					//for (int j = 0; j < CutList.size(); ++j)
					//{
					//	std::cout << CutList[j] << " ";
					//}
					expr.clear();

					for (IloInt i = 0; i < CutList.size() - 1; i++)
					{
						expr += solver.x[CutList[i]][CutList[i+1]][0] + solver.x[CutList[i + 1]][CutList[i]][0];
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					expr2.clear();
					expr.clear();

					for (int i = 0; i < CutList.size()-1; i++)
					{
						for (int j = i; j < CutList.size(); j++)
						{
							expr += solver.x[CutList[i]][CutList[j]][0];
							expr2 += solver.x[CutList[j]][CutList[i]][0];
						}
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					add(expr2 <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_SET)
				{
					//consName = "NCS -> ";
				/*printf("\nEnergy problem");*/
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->NoChargeSet++;
						solver.Stats->numberOfCutsAtRootNode->NoChargeSet++;
					}
					else solver.Stats->totalNumberOfCuts->NoChargeSet++;

					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
					}

					//std::cout << "\nSet: ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						//std::cout << CutList[j - 1] << " ";
					}
					expr.clear();

					for (IloInt i = 0; i < CutList.size(); i++)
					{
						for (IloInt j = 0; j < CutList.size(); j++)
						{
							if (i < j)
							{
								expr += solver.x[CutList[i]][CutList[j]][0] + solver.x[CutList[j]][CutList[i]][0];
							}
						}
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_CHARGER_PATH) {
					//consName = "CP -> ";
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->ChargerPath++;
						solver.Stats->numberOfCutsAtRootNode->ChargerPath++;
					}
					else solver.Stats->totalNumberOfCuts->ChargerPath++;
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
					}

					expr.clear();
					LHS = 0;
					//Find arcs to add
					for (int i = 0; i < CutList.size() - 1; i++)
					{
						if (CutList[i] >= n) {
							for (int j = 0; j < n; j++)
							{
								if (std::find(CutList.begin() + 1, CutList.end() - 1, j) == CutList.end() - 1)
								{
									for (int r = 0; r < solver.R[CutList[i + 1]][j].size(); r++)
									{
										if (solver.R[CutList[i + 1]][j][r] == CutList[i])
										{
											expr += solver.y[j][CutList[i + 1]][r] + solver.y[CutList[i + 1]][j][r];
										}
									}
								}
								

							}
							LHS++;
						}
						else if (CutList[i + 1] >= n)
						{
							for (int j = 0; j < n; j++)
							{
								if (std::find(CutList.begin() + 1, CutList.end() - 1, j) == CutList.end() - 1)
								{
									for (int r = 0; r < solver.R[CutList[i]][j].size(); r++)
									{
										if (solver.R[CutList[i]][j][r] == CutList[i + 1])
										{
											expr += solver.y[CutList[i]][j][r] + solver.y[j][CutList[i]][r];
										}
									}
								}

							}
							LHS++;
						}
						else {
							expr += solver.x[CutList[i]][CutList[i + 1]][0] + solver.x[CutList[i + 1]][CutList[i]][0];
							LHS++;
						}
					}
					add(expr <=LHS-1);
					//solver.infeas.add(expr <= LHS-1);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_CHARGERS) { //ADD DEPOT CONSTRAINTS
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->FixedChargers++;
						solver.Stats->numberOfCutsAtRootNode->FixedChargers++;
					}
					else solver.Stats->totalNumberOfCuts->FixedChargers++;
					//consName = "FC->";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
						//printf("%d ", CutList[j - 1]);
					}

					expr.clear();

					//Add arcs connecting S:
					for (int i = 1; i < CutList.size() - 1; i++)
						for (int j = 1; j < CutList.size() - 1; j++)
							expr += 3 * solver.x[CutList[i]][CutList[j]][0]+ 3 * solver.x[CutList[i]][CutList[j]][1];

					//Add connection from r1 to S
					if (CutList[0] > 0) {
						for (int i = 0; i < n; i++)
							for (int j = 1; j < CutList.size() - 1; j++)
								for (int r = 0; r < solver.R[i][CutList[j]].size(); r++)
									if (CutList[0] == solver.R[i][CutList[j]][r]) expr += solver.y[i][CutList[j]][r];
					}
					else
					{
						for (int j = 1; j < CutList.size() - 1; j++)
							expr += solver.x[0][CutList[j]][0];
					}

					//Add connection from S to r2
					if (CutList.back() > 0) {
						for (int i = 0; i < n; i++)
							for (int j = 1; j < CutList.size() - 1; j++)
								for (int r = 0; r < solver.R[i][CutList[j]].size(); r++)
									if (CutList.back() == solver.R[i][CutList[j]][r]) expr += solver.y[CutList[j]][i][r];
					}
					else
					{
						for (int j = 1; j < CutList.size() - 1; j++)
							expr += solver.x[CutList[j]][0][0];
					}


					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_CHARGER_SINGLE) {
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->FixedChargerSingle++;
						solver.Stats->numberOfCutsAtRootNode->FixedChargerSingle++;
					}
					else solver.Stats->totalNumberOfCuts->FixedChargerSingle++;
					//consName = "FCS->";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
					}

					//Add arcs connecting S
					expr.clear();
					for (int i = 1; i < CutList.size(); i++)
						for (int j = 1; j < CutList.size(); j++)
							expr += 2 * solver.x[CutList[i]][CutList[j]][0] + 2 * solver.x[CutList[i]][CutList[j]][1];

					
					if (CutList[0] > 0) { //If the charging option is an external charger
						for (int i = 0; i < n; i++)
							for (int j = 1; j < CutList.size(); j++)
								for (int r = 0; r < solver.R[i][CutList[j]].size(); r++)
									if (solver.R[i][CutList[j]][r] == CutList[0]) expr += solver.y[i][CutList[j]][r];
					}
					else { //If the charging option is the depot
						for (int i = 1; i < CutList.size(); i++)
						{
							expr += solver.x[0][CutList[i]][0];
						}
					}

					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_EDGES) { //ADD FIXED EDGES CONSTRAINTS
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->FixedEdges++;
						solver.Stats->numberOfCutsAtRootNode->FixedEdges++;
					}
					else solver.Stats->totalNumberOfCuts->FixedEdges++;
					//std::cout << "\nList Tail: ";
						//consName = "FE -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";

						//std::cout << CutList[j-1] << " ";
					}
					//std::cout << "\nPath: ";
					//for (int j = 0; j < List.size(); ++j)
					//{
					//	std::cout << List[j] << " ";
					//}
					double r1_dist = solver.c[CutList[0]][CutList[1]];
					double r2_dist = solver.c[CutList[CutList.size() - 2]][CutList[CutList.size() - 1]];

					expr.clear();
					//Add all arcs directly connecting S
					for (IloInt i = 1; i < CutList.size() - 1; i++)
						for (IloInt j = 1; j < CutList.size() - 1; j++)
							expr += solver.x[CutList[i]][CutList[j]][0];


					for (int i = 0; i < n; i++)
					{

						if (std::find(CutList.begin() + 1, CutList.end() - 1, i) != CutList.end() - 1) continue;
						if (solver.c[i][CutList[1]] + solver.MinCharDist[i] >= r1_dist) expr += solver.x[i][CutList[1]][0];// +solver.x[CutList[1]][i][0];
						if (solver.c[CutList[CutList.size() - 2]][i] + solver.MinCharDist[i] >= r2_dist) expr += solver.x[CutList[CutList.size() - 2]][i][0];// +solver.x[i][CutList[CutList.size() - 2]][0];

						for (int r = 0; r < solver.R[i][CutList[1]].size(); r++)
						{
							if (solver.c[solver.R[i][CutList[1]][r]][CutList[1]] >= r1_dist) expr += solver.y[i][CutList[1]][r];// + solver.y[CutList[1]][i][r];
						}
						for (int r = 0; r < solver.R[CutList[CutList.size() - 2]][i].size(); r++)
						{
							if (solver.c[CutList[CutList.size() - 2]][solver.R[CutList[CutList.size() - 2]][i][r]] >= r2_dist) expr += solver.y[CutList[CutList.size() - 2]][i][r];// +solver.y[i][CutList[CutList.size() - 2]][r];
						}
					}

					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_SINGLEEDGE) { //ADD FIXED SINGLE EDGE CONSTRAINTS

					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->FixedSingleEdge++;
						solver.Stats->numberOfCutsAtRootNode->FixedSingleEdge++;
					}
					else solver.Stats->totalNumberOfCuts->FixedSingleEdge++;
					//consName = "FES -> ";
				//std::cout << "\nList Tail: ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";

						//std::cout << ListTail[j-1] << " ";
					}
					//std::cout << "\nPath: ";
					//for (int j = 0; j < List.size(); ++j)
					//{
					//	std::cout << List[j] << " ";
					//}
					double r1_dist = solver.c[CutList[0]][CutList[1]];
					expr.clear();
					for (IloInt i = 1; i < CutList.size(); i++)
					{
						for (IloInt j = 1; j < CutList.size(); j++)
						{
							expr += solver.x[CutList[i]][CutList[j]][0];
						}


					}

					for (int i = 0; i < n; i++)
					{
						if (std::find(CutList.begin() + 1, CutList.end(), i) != CutList.end()) continue;
						if (solver.c[i][CutList[1]] + solver.MinCharDist[i] >= r1_dist) expr += solver.x[i][CutList[1]][0];// +solver.x[CutList[1]][i][0];
						for (int r = 0; r < solver.R[i][CutList[1]].size(); r++)
						{
							if (solver.c[solver.R[i][CutList[1]][r]][CutList[1]] >= r1_dist) expr += solver.y[i][CutList[1]][r];// +solver.y[CutList[1]][i][r];
						}
					}


					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());



				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_INFEASIBLE_ICEV_PATH)
				{
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->TimeInfeasibleICEVPath++;
						solver.Stats->numberOfCutsAtRootNode->TimeInfeasibleICEVPath++;
					}
					else solver.Stats->totalNumberOfCuts->TimeInfeasibleICEVPath++;
					//std::cout << "\nList Tail: ";
					//consName = "Time_ICEV -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";

						//std::cout << ListTail[j-1] << " ";
					}
					//std::cout << "\nPath: ";
					//for (int j = 0; j < CutList.size(); ++j)
					//{
					//	std::cout << CutList[j] << " ";
					//}
					expr.clear();

					for (IloInt i = 0; i < CutList.size() - 1; i++)
					{

						expr += solver.x[CutList[i]][CutList[i + 1]][0] + solver.x[CutList[i]][CutList[i + 1]][1] + solver.x[CutList[i + 1]][CutList[i]][0] + solver.x[CutList[i + 1]][CutList[i]][1];
						for (int r = 0; r < solver.R[CutList[i]][CutList[i + 1]].size(); r++) expr += solver.y[CutList[i]][CutList[i + 1]][r] + solver.y[CutList[i + 1]][CutList[i]][r];


					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_TIME_INFEASIBLE_SET)
				{
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->TimeInfeasibleSet++;
						solver.Stats->numberOfCutsAtRootNode->TimeInfeasibleSet++;
					}
					else solver.Stats->totalNumberOfCuts->TimeInfeasibleSet++;
					//consName = "Time_set -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";

					}

					//std::cout << "\nSet: ";
					//for (int j = 1; j <= MyCutsCMP->CPL[c]->IntListSize; ++j)
					//{
					//	std::cout << CutList[j - 1] << " ";
					//}
					expr.clear();

					for (IloInt i = 0; i < CutList.size(); i++)
					{
						for (IloInt j = 0; j < CutList.size(); j++)
						{
							if (i < j)
							{
								expr += solver.x[CutList[i]][CutList[j]][0] + solver.x[CutList[j]][CutList[i]][0]+ solver.x[CutList[i]][CutList[j]][1] + solver.x[CutList[j]][CutList[i]][1];
								for (int r = 0; r < solver.R[CutList[i]][CutList[j]].size(); r++) expr += solver.y[CutList[i]][CutList[j]][r] + solver.y[CutList[j]][CutList[i]][r];
							}
						}
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());


				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_INFEASIBLE_EV_PATH) {

					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->TimeInfeasibleEVPath++;
						solver.Stats->numberOfCutsAtRootNode->TimeInfeasibleEVPath++;
					}
					else solver.Stats->totalNumberOfCuts->TimeInfeasibleEVPath++;
					//consName = "Time_EV ->";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
						
						//printf("%d ", CutList[j - 1]);
					}

					expr.clear();
					LHS = 0;
					//Find arcs to add
					for (int i = 0; i < CutList.size() - 1; i++)
					{
						if (CutList[i] >= n && i == 0)
						{
							for (int j = 0; j < n; j++)
							{
								for (int r = 0; r < solver.R[j][CutList[i + 1]].size(); r++)
								{
									if (solver.R[j][CutList[i + 1]][r] == CutList[i]) expr += solver.y[j][CutList[i + 1]][r] + solver.y[CutList[i + 1]][j][r];
								}
							}
							LHS++;
						}
						else if (CutList[i + 1] >= n && i + 1 == CutList.size() - 1) {
							for (int j = 0; j < n; j++)
							{
								for (int r = 0; r < solver.R[CutList[i]][j].size(); r++)
								{
									if (solver.R[CutList[i]][j][r] == CutList[i + 1]) expr += solver.y[CutList[i]][j][r] + solver.y[j][CutList[i]][r];
								}
							}
							LHS++;
						}
						else if (CutList[i + 1] >= n)
						{
							for (int r = 0; r < solver.R[CutList[i]][CutList[i + 2]].size(); r++)
							{
								if (solver.R[CutList[i]][CutList[i + 2]][r] == CutList[i + 1]) { expr += solver.y[CutList[i]][CutList[i + 2]][r] + solver.y[CutList[i + 2]][CutList[i]][r]; break; }
							}
							i++;
							LHS++;
						}
						else {
							expr += solver.x[CutList[i]][CutList[i + 1]][0] + solver.x[CutList[i + 1]][CutList[i]][0];
							LHS++;
						}
					}
					add(expr <= LHS - 1);
					//solver.infeas.add(expr <= LHS - 1);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_TIME_INFEASIBLE_EV_SET)
				{
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->TimeInfeasibleEVSet++;
						solver.Stats->numberOfCutsAtRootNode->TimeInfeasibleEVSet++;
					}
					else solver.Stats->totalNumberOfCuts->TimeInfeasibleEVSet++;
					//consName = "Time_EV_Set";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";

						//printf("%d ", CutList[j - 1]);
					}
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->ExtListSize; ++j)
					{
						ExtList.push_back(MyCutsCMP->CPL[cons]->ExtList[j]);
						//consName = consName + std::to_string(ExtList[j - 1]) + " ";

						//printf("%d ", CutList[j - 1]);
					}
					
					//Arcs connecting nodes in S
					expr.clear();
					for (int i = 0; i < CutList.size(); i++)
					{
						for (int j = 0; j < CutList.size(); j++)
						{
							expr += solver.x[CutList[i]][CutList[j]][0];
							for (int r = 0; r < solver.R[CutList[i]][CutList[j]].size(); r++)
							{
								if (std::find(ExtList.begin(),ExtList.end(),solver.R[CutList[i]][CutList[j]][r]) != ExtList.end())
								{
									expr += 2 * solver.y[CutList[i]][CutList[j]][r];
								}
							}
							
						}
					}

					for (int i = 0; i < n; i++)
					{
						if (std::find(CutList.begin(),CutList.end(),i)==CutList.end())
						{
							for (int j = 0; j < CutList.size(); j++)
							{
								for (int r = 0; r < solver.R[i][CutList[j]].size(); r++)
								{
									if (std::find(ExtList.begin(), ExtList.end(), solver.R[i][CutList[j]][r]) != ExtList.end())
									{
										expr += solver.y[i][CutList[j]][r];
									}
								}
								for (int r = 0; r < solver.R[CutList[j]][i].size(); r++)
								{
									if (std::find(ExtList.begin(), ExtList.end(), solver.R[CutList[j]][i][r]) != ExtList.end())
									{
										expr += solver.y[CutList[j]][i][r];
									}
								}
							}
						}
						
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
				}
			}


			for (i = 0; i < MyCutsCMP->Size; i++)
			{
				CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
			}

			MyCutsCMP->Size = 0;
			end = solver.Cplex.getTime();
			solver.Stats->AddCutsTime += end - StartAddingCuts;
			//if (solver.DrawCounter % 10 == 0)
			//{
			//	std::string strFilename = "Grafer/Tid/Elbils_løsning_n21_f9_tid_" + std::to_string(solver.DrawCounter) + ".tex";
			//	char* FileName;
			//	FileName = &strFilename[0];
			//	std::string strGraphname = "Elbils loesning (tidsforbrug) eksempel " + std::to_string(solver.DrawCounter);
			//	char* GraphName = &strGraphname[0];
			//	DRAWGRAPH_DrawCVRPSparseGraph(n - 1, f, EdgeTailDraw.data(), 0, solver.xCoord.data(), solver.yCoord.data(), EdgeDraw.size(), EdgeTailDraw.data(), EdgeHeadDraw.data(), EnergyComps.data(), TravelTime.data(), d, s, EdgeDraw.data(), 0.7, GraphName, FileName, true);
			//	strFilename = "Grafer/Energi/Elbils_løsning_n21_f9_energi_" + std::to_string(solver.DrawCounter) + ".tex";
			//	FileName = &strFilename[0];
			//	strGraphname = "Elbils loesning (energiforbrug) eksempel " + std::to_string(solver.DrawCounter);
			//	GraphName = &strGraphname[0];
			//	DRAWGRAPH_DrawCVRPSparseGraph(n - 1, f, EdgeTailDraw.data(), 0, solver.xCoord.data(), solver.yCoord.data(), EdgeDraw.size(), EdgeTailDraw.data(), EdgeHeadDraw.data(), EnergyComps.data(), TravelTime.data(), d, s, EdgeDraw.data(), 0.7, GraphName, FileName, false);
			//}
		}
	}
	


	delete[] CList;
	delete[] TList;
	delete[] NList;
	delete[] Label;
	delete[] d;
	delete[] s;

	if (NNodes == 0) {
		solver.Stats->LBAtRootNode = getBestObjValue();
		solver.Stats->timeRootNode = solver.Cplex.getTime();
	}

	if (solver.Stats->treeDepth < getCurrentNodeDepth()) solver.Stats->treeDepth = getCurrentNodeDepth();
	endTot = solver.Cplex.getTime();
	solver.Stats->TotalSepTime += endTot - startTot;
}
// LAZY CUTS::
//RIGHT NOW CAP CUTS ADDED TOGETHER WITH BATTERY CUTS. CAP CUTS WORKS AS INTENDED, BATTERY CUTS NEEDS TO BE CHECKED FOR WHETHER THEY ALWAYS REACH CHARGING STATION BEFORE RUNNING OUT OF BATTERY:
ILOLAZYCONSTRAINTCALLBACK1(LazyCapCuts, MFGVRP_Solver&, solver) {
	double start = 0;
	double end = 0;
	double startTot = solver.Cplex.getTime();
	double endTot = 0;
	double startAddingCuts = 0;
	IloInt i;
	IloEnv env = solver.getEnv();
	IloInt n = solver.x.getSize();
	bool edgesNeedFix = false;
	bool isRoot = false;
	//std::string consName;
	double edgeVal;
	double edgeValX;
	double edgeValY;
	IloExpr expr(env);
	IloExpr expr2(env);
	//std::string consName;
	int* d = new int[n];
	int Q = solver.getCapacity();
	int* s = new int[n];
	int T = solver.getMaxTime();

	bool CutsAdded = false;
	//std::vector<float> EnergyComps;
	//std::vector<float>TravelTime;
	std::vector<int> CutList;
	std::vector<int> ExtList;
	int NNodes = getNnodes(); 
	//int NNodes = 0;
	if (NNodes == 0)isRoot = true;
	
	int LHS = 0;
	//int* EdgeHead1 = new int[ceil(n * n / 2)];
	//int* EdgeTail1 = new int[ceil(n * n / 2)];
	//double* EdgeX1 = new double[ceil(n * n / 2)];
	int NoOfChargingEdges = 0;
	//For FCI:
	//int* Label = new int[n + 1];
	//int MaxIdx;
	//int MinIdx;
	//int k;


	//For Multistar:
	//int* NList = new int[n + 1];
	//int* TList = new int[n + 1];
	//int* CList = new int[n + 1];
	//int intA;
	//int intB;
	//int intL;

	std::vector<int> EdgeTail; EdgeTail.reserve(n * 5); EdgeTail.push_back(0);
	std::vector<int> EdgeHead; EdgeHead.reserve(n * 5); EdgeHead.push_back(0);
	std::vector<double> EdgeX; EdgeX.reserve(n * 5); EdgeX.push_back(-1);
	//std::vector<int> EdgeTailY;
	//std::vector<int> EdgeHeadY;
	//std::vector<double> EdgeY;
	////std::vector<double> EdgeDraw;
	////std::vector<int> EdgeTailDraw;
	////std::vector<int> EdgeHeadDraw;
	//std::vector<int> EdgeChargStation;
	//EdgeTail.push_back(0);
	//EdgeHead.push_back(0);
	//EdgeX.push_back(0);

	////One vector setup
	//std::vector<double> temp_EdgeX;
	//std::vector<int> temp_EdgeTail;
	//std::vector<int> temp_EdgeHead;

	std::vector<int> EdgeTail_ICEV; EdgeTail_ICEV.reserve(n * 5); EdgeTail_ICEV.push_back(0);
	std::vector<int> EdgeHead_ICEV; EdgeHead_ICEV.reserve(n * 5); EdgeHead_ICEV.push_back(0);
	std::vector<double> EdgeX_ICEV; EdgeX_ICEV.reserve(n * 5); EdgeX_ICEV.push_back(-1);
	std::vector<int> EdgeTail_Charger; EdgeTail_Charger.reserve(n * 5); EdgeTail_Charger.push_back(0);
	std::vector<int> EdgeHead_Charger; EdgeHead_Charger.reserve(n * 5); EdgeHead_Charger.push_back(0);
	std::vector<double> EdgeX_Charger; EdgeX_Charger.reserve(n * 5); EdgeX_Charger.push_back(-1);
	std::vector<int> EdgeTail_Direct; EdgeTail_Direct.reserve(n * 5); EdgeTail_Direct.push_back(0);
	std::vector<int> EdgeHead_Direct; EdgeHead_Direct.reserve(n * 5); EdgeHead_Direct.push_back(0);
	std::vector<double> EdgeX_Direct; EdgeX_Direct.reserve(n * 5); EdgeX_Direct.push_back(-1);
	std::vector<int> Origin; Origin.reserve(n * 5); Origin.push_back(0);

	int ChargingArcs = 0;
	int	DepotArcs = 0;
	int DirectArcs = 0;
	int ICEVArcs = 0;
	int ICEVDepotArcs = 0;
	int OriginID = 1;
	std::vector<int> reduced_EdgeTail;
	std::vector<int> reduced_EdgeHead;
	std::vector<double> reduced_EdgeX;


	//CVRPSEP parameters:
	int ChargingStation;
	CnstrMgrPointer MyCutsCMP, MyOldCutsCMP;

	//solver.Create_Cut_Mgr(&CutSet, 100);
	CMGR_CreateCMgr(&MyCutsCMP, 100);
	CMGR_CreateCMgr(&MyOldCutsCMP, 100);

	IloNum3DMatrix xTemp(env, n);
	IloNum3DMatrix yTemp(env, n);


	IloNumArray xTemp1d(env, solver.xDummy.getSize());
	IloNumArray yTemp1d(env, solver.yDummy.getSize());

	start = solver.Cplex.getTime();
	getValues(xTemp1d, solver.xDummy);
	getValues(yTemp1d, solver.yDummy);
	end = solver.Cplex.getTime();
	solver.Stats->TotalLoadSepData += end - start;
	start = solver.Cplex.getTime();
	int cntX = 0;
	int cntY = 0;
	for (int i = 0; i < n; i++)
	{
		xTemp[i] = IloNum2DMatrix(env, n);
		yTemp[i] = IloNum2DMatrix(env, n);
		for (int j = 0; j < n; j++)
		{
			xTemp[i][j] = IloNumArray(env, 2);
			for (int k = 0; k < 2; k++)
			{
				xTemp[i][j][k] = xTemp1d[cntX];
				cntX++;
			}
			yTemp[i][j] = IloNumArray(env, solver.R[i][j].size());
			for (int r = 0; r < solver.R[i][j].size(); r++)
			{
				yTemp[i][j][r] = yTemp1d[cntY];
				cntY++;
			}
		}
	}
	xTemp1d.end();
	yTemp1d.end();
	end = solver.Cplex.getTime();
	solver.Stats->ConvertVals += end - start;
	start = solver.Cplex.getTime();
	for (IloInt i = 0; i < n; i++)
	{
		if (i > 0) { d[i] = solver.getDemand(i);  s[i] = solver.getServiceTime(i); }
		
		for (IloInt j = i + 1; j < n; j++)
		{
			edgeVal = 0;
			expr.clear();

			for (IloInt r = 0; r < solver.R[i][j].size(); r++)
			{
				edgeValY = yTemp[i][j][r] + yTemp[j][i][r];
				edgeVal += edgeValY;
				expr += solver.y[i][j][r] + solver.y[j][i][r];
				if (edgeValY > 0.001)
				{
					//One vector setup:
					if (yTemp[i][j][r] > 0.001)
					{
						EdgeTail_Charger.push_back(i);
						EdgeHead_Charger.push_back(solver.R[i][j][r]);
						EdgeX_Charger.push_back(yTemp[i][j][r]);
						Origin.push_back(OriginID);
						EdgeTail_Charger.push_back(solver.R[i][j][r]);
						EdgeHead_Charger.push_back(j);
						EdgeX_Charger.push_back(yTemp[i][j][r]);
						Origin.push_back(OriginID);
						OriginID++;
					}
					if (yTemp[j][i][r] > 0.001)
					{
						EdgeTail_Charger.push_back(j);
						EdgeHead_Charger.push_back(solver.R[j][i][r]);
						EdgeX_Charger.push_back(yTemp[j][i][r]);
						Origin.push_back(OriginID);
						EdgeTail_Charger.push_back(solver.R[j][i][r]);
						EdgeHead_Charger.push_back(i);
						EdgeX_Charger.push_back(yTemp[j][i][r]);
						Origin.push_back(OriginID);
						OriginID++;
					}
				}
			}
			edgeValX = xTemp[i][j][1] + xTemp[j][i][1] + xTemp[i][j][0] + xTemp[j][i][0];
			edgeVal += edgeValX;
			expr += solver.x[i][j][0] + solver.x[j][i][0] + solver.x[i][j][1] + solver.x[j][i][1];


			if (edgeValX > 0.001) {
				//One vector setup
				if (xTemp[i][j][0] > 0.001) {
					if (i == 0) {
						EdgeTail_Charger.push_back(i);
						EdgeHead_Charger.push_back(j);
						EdgeX_Charger.push_back(xTemp[i][j][0]);
						Origin.push_back(i);
					}
					else {
						EdgeTail_Direct.push_back(i);
						EdgeHead_Direct.push_back(j);
						EdgeX_Direct.push_back(xTemp[i][j][0]);
					}

				}
				if (xTemp[j][i][0] > 0.001) {
					if (i == 0) {
						EdgeTail_Charger.push_back(j);
						EdgeHead_Charger.push_back(i);
						EdgeX_Charger.push_back(xTemp[j][i][0]);
						Origin.push_back(i);

					}
					else {
						EdgeTail_Direct.push_back(j);
						EdgeHead_Direct.push_back(i);
						EdgeX_Direct.push_back(xTemp[j][i][0]);
					}


				}
				if (xTemp[i][j][1] > 0.001) {
					EdgeTail_ICEV.push_back(i);
					EdgeHead_ICEV.push_back(j);
					EdgeX_ICEV.push_back(xTemp[i][j][1]);
				}
				if (xTemp[j][i][1] > 0.001) {
					EdgeTail_ICEV.push_back(j);
					EdgeHead_ICEV.push_back(i);
					EdgeX_ICEV.push_back(xTemp[j][i][1]);
				}
			}

			if (edgeVal > 0)
			{

				EdgeTail.push_back(i == 0 ? n : i);
				EdgeHead.push_back(j);
				EdgeX.push_back(edgeVal);
				if (edgeVal > 1.001 && i != 0)
				{
					edgesNeedFix = true;
					add(expr <= 1);
					//solver.infeas.add(expr <= 1);
				}
			}
		}
	}

	//for (int i = 0; i < n; i++)
	//{
	//	for (int j = 0; j < n; j++)
	//	{
	//		xTemp[i][j].end();
	//		yTemp[i][j].end();
	//	}
	//	xTemp[i].end();
	//	yTemp[i].end();
	//}
	//xTemp.end();
	//yTemp.end();
	NoOfChargingEdges = EdgeTail_Charger.size();
	end = solver.Cplex.getTime();
	solver.Stats->SetupSepVectors += end - start;


	//std::vector<int> p1 = solver.slicing(temp_EdgeTail, DepotArcs, DepotArcs + DirectArcs - 1);
	//std::vector<int> p2 = solver.slicing(temp_EdgeHead, DepotArcs, DepotArcs + DirectArcs - 1);
	//std::vector<double> p3 = solver.slicing(temp_EdgeX, DepotArcs, DepotArcs + DirectArcs - 1);
	//printf("\n\nReduced vector\n");
	//printf("\nDepot arcs: %d\tDirect arcs: %d\tCharging arcs: %d\tICEV arcs:%d\n\n", DepotArcs, DirectArcs, ChargingArcs,ICEVArcs);
	//for (int i = 0; i < temp_EdgeHead.size(); i++)
	//{
	//	printf("%d\t%d\t%f\n", temp_EdgeTail[i], temp_EdgeHead[i],temp_EdgeX[i]);
	//}

	//for (int i = 0; i < p1.size(); i++)
	//{
	//	printf("%d\t%d\t%f\n", p1[i], p2[i], p3[i]);
	//}
	//Input for Graph:
	//std::cout << "\n\nMy draw values:\n";
	//for (IloInt i = 0; i < EdgeTailDraw.size(); i++)
	//{
	//	std::cout << "\n" << EdgeTailDraw[i] << "\t" << EdgeHeadDraw[i] << "\t" << EdgeDraw[i] << "\t" << EnergyComps[i] << "\t" << TravelTime[i];
	//}
	//std::cout << "\n\nMy draw values:\n";
	//std::cout << "\n\nMy load cuts values:\n";
	//for (IloInt i = 1; i <= solver.NoOfEdges; i++)
	//{
	//	std::cout << "\n" << EdgeTail[i] << " " << EdgeHead[i] << " " << EdgeX[i];
	//}
	//

	//std::cout << "\n\nMy energy cuts values:\n";
	//for (IloInt i = 0; i < EdgeTailY.size(); i++)
	//{
	//	std::cout << "\n" << EdgeTailY[i] << " " << EdgeHeadY[i] << " " << EdgeY[i] << " " << EdgeChargStation[i];
	//}

	//std::cout << "\n\nMy draw values:\n";
	//for (IloInt i = 0; i < EdgeTailDraw.size(); i++)
	//{
	//	std::cout << "\n" << EdgeTailDraw[i] << "\t" << EdgeHeadDraw[i] << "\t" << EdgeDraw[i];
	//}
	//std::cout << "\n\nx_ICEV:\n";
	//for (IloInt i = 0; i < n; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n; j++)
	//	{
	//		std::cout << getValue(solver.x[i][j][1])<<" ";
	//	}
	//}

	//std::cout << "\n\nx_EV:\n";
	//for (IloInt i = 0; i < n; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n; j++)
	//	{
	//		std::cout << getValue(solver.x[i][j][0]) << " ";
	//	}
	//}

	//std::cout << "\n\ny_EV:\n";
	//for (IloInt i = 0; i < n; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n; j++)
	//	{
	//		std::cout << "{";
	//		for (IloInt r = 0; r < solver.R[i][j].size(); r++)
	//		{
	//			std::cout << getValue(solver.y[i][j][r]) << ";";
	//		}
	//		std::cout << "} ";
	//	}
	//}

	if (!edgesNeedFix) //If no edges needs to be fixed, we can initiate our seperation problem:
	{
		start = solver.Cplex.getTime();
		int cntCut = MyCutsCMP->Size;
		CAPSEP_SeparateCapCuts(n - 1,
			d,
			Q,
			EdgeTail.size()-1,
			EdgeTail.data(),
			EdgeHead.data(),
			EdgeX.data(),
			MyOldCutsCMP,
			solver.MaxNoOfCuts,
			solver.EpsForIntegrality,
			&solver.IntegerAndFeasible,
			&solver.MaxViolation,
			MyCutsCMP);

		if (MyCutsCMP->Size==0)
		{
			CAPSEP_SeparateCapCuts(n - 1,
				s,
				T,
				EdgeTail.size() - 1,
				EdgeTail.data(),
				EdgeHead.data(),
				EdgeX.data(),
				MyOldCutsCMP,
				solver.MaxNoOfCuts,
				solver.EpsForIntegrality,
				&solver.IntegerAndFeasible,
				&solver.MaxViolation,
				MyCutsCMP);
		}
		end = solver.Cplex.getTime();
		solver.Stats->CapSepTime += end - start;
		/*solver.ObjValTemporary = getObjValue();*/

		/////*solver.DrawCounter++;*/
		//if (solver.ObjValTemporary > 840.94 && solver.ObjValTemporary < 840.95)
		//{
		//	std::string strFilename = "Grafer/EVRP_SOL.tex";
		//	char* FileName;
		//	FileName = &strFilename[0];
		//	char GraphName[] = "Elbils loesning (sparse)";
		//	DRAWGRAPH_DrawCVRPSparseGraph(n - 1, f, EdgeTailDraw.data(), 0, solver.xCoord.data(), solver.yCoord.data(), EdgeDraw.size(), EdgeTailDraw.data(), EdgeHeadDraw.data(), EnergyComps.data(), TravelTime.data(), d, s, EdgeDraw.data(), 0.7, GraphName, FileName, false);
		
		//}

		if (MyCutsCMP->Size == 0) {
			//if (NNodes>0)
			//{

			start = solver.Cplex.getTime();

			solver.CutCounter = 0;

			//reduced_EdgeTail.clear();
			//reduced_EdgeHead.clear();
			//reduced_EdgeX.clear();
			//reduced_EdgeTail.push_back(0);
			//reduced_EdgeHead.push_back(0);
			//reduced_EdgeX.push_back(0);

			//for (int i = DepotArcs + 1; i <= DepotArcs + DirectArcs; i++)
			//{
			//	//reduced_EdgeTail.push_back(temp_EdgeTail[i]);
			//	//reduced_EdgeHead.push_back(temp_EdgeHead[i]);
			//	//reduced_EdgeX.push_back(temp_EdgeX[i]);
			//	printf("\n%d\t%d\t%f", temp_EdgeTail[i], temp_EdgeHead[i], temp_EdgeX[i]);
			//}
			//
			//printf("\nNew");
			//for (int i = 1; i <EdgeHead_Direct.size(); i++)
			//{
			//	printf("\n%d\t%d\t%f", EdgeTail_Direct[i], EdgeHead_Direct[i], EdgeX_Direct[i]);
			//}

			//printf("\nNo charge");
			cntCut = MyCutsCMP->Size;
			////printf("\n start violated paths");
			////printf("\We stated no-charge user");
			//solver.FindViolatedPaths(reduced_EdgeTail.size() - 1, reduced_EdgeTail, reduced_EdgeHead, reduced_EdgeX, MyCutsCMP, isRoot);
			solver.FindViolatedPaths(EdgeTail_Direct.size() - 1, EdgeTail_Direct, EdgeHead_Direct, EdgeX_Direct, MyCutsCMP, isRoot);

			if (cntCut != MyCutsCMP->Size) solver.NodeID = getNodeId();
			//printf("\We finished no-charge user");
			if (MyCutsCMP->Size==0)
			{
				EdgeTail_Charger.insert(EdgeTail_Charger.end(), EdgeTail_Direct.begin() + 1, EdgeTail_Direct.end());
				EdgeHead_Charger.insert(EdgeHead_Charger.end(), EdgeHead_Direct.begin() + 1, EdgeHead_Direct.end());
				EdgeX_Charger.insert(EdgeX_Charger.end(), EdgeX_Direct.begin() + 1, EdgeX_Direct.end());
			}
			//printf("\n end violated paths");			
			//if (NNodes == 0 && MyCutsCMP->Size == 0)
			//{
			//	
			//	solver.AgressiveSetSearch_FC(EdgeTail_Charger.size(), NoOfChargingEdges, EdgeTail_Charger.data(), EdgeHead_Charger.data(), EdgeX_Charger.data(), MyCutsCMP);
			//}
			if (MyCutsCMP->Size == 0)
			{
				cntCut = MyCutsCMP->Size;
				/*temp_EdgeTail.insert(temp_EdgeTail.begin(), 0);
				temp_EdgeHead.insert(temp_EdgeHead.begin(), 0);
				temp_EdgeX.insert(temp_EdgeX.begin(), 0);*/
				//printf("\n start violated Charger Paths");
				//printf("\nWe started Fixed user");

				//printf("\nSep Charger");
				//for (int i = 0; i < EdgeTail_Charger.size(); i++)
				//{
				//	printf("\n%d\t%d\t%f", EdgeTail_Charger[i], EdgeHead_Charger[i], EdgeX_Charger[i]);
				//}
				//printf("\nSep Direct");
				//for (int i = 0; i < EdgeTail_Direct.size(); i++)
				//{
				//	printf("\n%d\t%d\t%f", EdgeTail_Direct[i], EdgeHead_Direct[i], EdgeX_Direct[i]);
				//}
				//printf("\nSep Consolidated");
				
				//solver.FixedChargers(temp_EdgeTail.size() - ICEVArcs, temp_EdgeTail, temp_EdgeHead, temp_EdgeX, MyCutsCMP, isRoot);
				//printf("\nfixed");
				//EdgeTail_Charger.reserve(EdgeTail_Charger.size() + EdgeTail_Direct.size() - 1);

				solver.FixedChargers(EdgeTail_Charger.size(), EdgeTail_Charger, EdgeHead_Charger, EdgeX_Charger, MyCutsCMP, isRoot);
				//for (int i = 0; i < EdgeTail_Charger.size(); i++)
				//{
				//	printf("\n%d\t%d\t%f", EdgeTail_Charger[i], EdgeHead_Charger[i], EdgeX_Charger[i]);
				//}
				//printf("\n end violated Charger Paths");
			}
			//}

			//Check EVs:
			if (MyCutsCMP->Size == 0)
			{


				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);

				//for (int i = DepotArcs + 1; i <= DepotArcs + DirectArcs + ChargingArcs; i++)
				//{
				//	//reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	//reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	//reduced_EdgeX.push_back(temp_EdgeX[i]);
				//	printf("\n%d\t%d\t%f", temp_EdgeTail[i], temp_EdgeHead[i], temp_EdgeX[i]);
				//}


				//printf("\nNew");
				//for (int i = 1; i <EdgeHead_Charger.size(); i++)
				//{
				//	printf("\n%d\t%d\t%f", EdgeTail_Charger[i], EdgeHead_Charger[i], EdgeX_Charger[i]);
				//}

				//printf("\ntime ev lazy");
				//cntCut = MyCutsCMP->Size;
				//if (getObjValue()>= 356.52 && getObjValue()<= 356.53)
				//{
				//	cntCut = cntCut;
				//}
				//solver.FindTimeInfeasiblePaths(reduced_EdgeTail.size() - 1, reduced_EdgeTail, reduced_EdgeHead, reduced_EdgeX, false, MyCutsCMP,isRoot);
				
				//if (MyCutsCMP->Size==0)
				//{
				//	printf("\n%f\n", getObjValue());

					//for (int i = 0; i < EdgeTail_Charger.size(); i++)
					//{
					//	printf("\n%d\t%d\t%f", EdgeTail_Charger[i], EdgeHead_Charger[i], EdgeX_Charger[i]);
					//}
					//int CurNode = 0;
					//int nRoutes = 0;
					//printf("\nRoutes:");
					//for (int i = 0; i < n; i++)
					//{
					//	if (xTemp[CurNode][i][0] > 0.9) {
					//		std::cout << "\n" << CurNode << " " << i;
					//		CurNode = i;
					//		for (int j = 0; j < n; j++)
					//		{
					//			if (xTemp[CurNode][j][0] > 0.9) {
					//				std::cout << " " << j;
					//				CurNode = j;
					//				j = -1;
					//			}
					//			if (j != -1)
					//			{
					//				for (int r = 0; r < solver.R[CurNode][j].size(); r++)
					//				{
					//					if (yTemp[CurNode][j][r] > 0.9)
					//					{
					//						std::cout << " " << solver.R[CurNode][j][r] << " " << j;
					//						CurNode = j;
					//						j = -1;
					//						break;
					//					}
					//				}
					//			}

					//			if (CurNode == 0)
					//			{
					//				//std::cout << "\n";
					//				nRoutes++;
					//				break;
					//			}
					//		}

					//	}
					//	for (int r1 = 0; r1 < solver.R[CurNode][i].size(); r1++)
					//	{
					//		if (yTemp[CurNode][i][r1] > 0.9) {
					//			std::cout << "\n" << CurNode << " " << solver.R[CurNode][i][r1] << " " << i;
					//			CurNode = i;
					//			for (int j = 0; j < n; j++)
					//			{
					//				if (xTemp[CurNode][j][0] > 0.9) {
					//					std::cout << " " << j;
					//					CurNode = j;
					//					j = -1;
					//				}
					//				if (j != -1) {
					//					for (int r = 0; r < solver.R[CurNode][j].size(); r++)
					//					{
					//						if (yTemp[CurNode][j][r] > 0.9)
					//						{
					//							std::cout << " " << solver.R[CurNode][j][r] << " " << j;
					//							CurNode = j;
					//							j = -1;
					//							break;
					//						}
					//					}
					//				}

					//				if (CurNode == 0)
					//				{
					//					nRoutes++;
					//					//std::cout << "\n";
					//					break;
					//				}

					//			}
					//		}
					//	}


				
					//}

				solver.FindTimeInfeasiblePaths(EdgeTail_Charger.size() - 1, EdgeTail_Charger, EdgeHead_Charger, EdgeX_Charger,Origin, false, MyCutsCMP, isRoot);
				//	if (nRoutes==1)
				//	{
				//if (MyCutsCMP->Size==0)
				//{
				//	solver.FixedChargers(EdgeTail_Charger.size(), EdgeTail_Charger, EdgeHead_Charger, EdgeX_Charger, MyCutsCMP, isRoot);
				//}
				////		solver.FindTimeInfeasiblePaths(EdgeTail_Charger.size() - 1, EdgeTail_Charger, EdgeHead_Charger, EdgeX_Charger, Origin, false, MyCutsCMP, isRoot);
				//	}
				//	
				//}
				//printf("\nend");
				bool CutWasAdded = false;
				if (0 != MyCutsCMP->Size) CutWasAdded=true;
				
				//if (MyCutsCMP->Size==0)
				//{
				//	std::cout << "\n" << getObjValue();
				//	solver.NodeID = getNodeId();
				//}

			}
			//Check ICEVS
			if (MyCutsCMP->Size == 0)
			{
				//cntCut = MyCutsCMP->Size;
				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);

				//for (int i = DepotArcs + DirectArcs + ChargingArcs + 1; i < temp_EdgeTail.size(); i++)
				//{
				//	reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	reduced_EdgeX.push_back(temp_EdgeX[i]);
				//}
				cntCut = MyCutsCMP->Size;

				//solver.FindTimeInfeasiblePaths(reduced_EdgeTail.size() - 1, reduced_EdgeTail, reduced_EdgeHead, reduced_EdgeX, true, MyCutsCMP,isRoot);
				//printf("\ntime icev");

				solver.FindTimeInfeasiblePaths(EdgeTail_ICEV.size() - 1, EdgeTail_ICEV, EdgeHead_ICEV, EdgeX_ICEV,Origin, true, MyCutsCMP, isRoot);
				solver.NodeID = getNodeId();

			}
			if (NNodes == 0 && MyCutsCMP->Size == 0)
			{

				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);
				//cntCut = MyCutsCMP->Size;
				//for (int i = DepotArcs + 1; i <= DepotArcs + DirectArcs; i++)
				//{
				//	reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	reduced_EdgeX.push_back(temp_EdgeX[i]);
				//	printf("\n%d\t%d\t%f", temp_EdgeTail[i], temp_EdgeHead[i], temp_EdgeX[i]);
				//}


				//printf("\nNew");
				//for (int i = 1; i <EdgeHead_Direct.size(); i++)
				//{
				//	printf("\n%d\t%d\t%f", EdgeTail_Direct[i], EdgeHead_Direct[i], EdgeX_Direct[i]);
				//}


				solver.AgressiveSetSearch(reduced_EdgeTail.size(), reduced_EdgeTail.data(), reduced_EdgeHead.data(), reduced_EdgeX.data(), MyCutsCMP, true, false);
				//printf("\nagg no charge");
				//solver.AgressiveSetSearch(EdgeTail_Direct.size(), EdgeTail_Direct.data(), EdgeHead_Direct.data(), EdgeX_Direct.data(), MyCutsCMP, true, false);

			}


			if (NNodes == 0 && MyCutsCMP->Size == 0)
			{
				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);

				//for (int i = DepotArcs + 1; i <= DepotArcs + DirectArcs + ChargingArcs; i++)
				//{
				//	reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	reduced_EdgeX.push_back(temp_EdgeX[i]);
				//}
				////printf("\nagg time ev");

				//cntCut = MyCutsCMP->Size;
				//solver.AgressiveSetSearch(reduced_EdgeTail.size(), reduced_EdgeTail.data(), reduced_EdgeHead.data(), reduced_EdgeX.data(), MyCutsCMP, false, false);
				solver.AgressiveSetSearch(EdgeTail_Charger.size(), EdgeTail_Charger.data(), EdgeHead_Charger.data(), EdgeX_Charger.data(), MyCutsCMP, false, false);
				//if (cntCut != MyCutsCMP->Size) solver.NodeID = getNodeId();
			}
			if (NNodes == 0 && MyCutsCMP->Size == 0)
			{

				//solver.NodeID = getNodeId();
				//reduced_EdgeTail.clear();
				//reduced_EdgeHead.clear();
				//reduced_EdgeX.clear();
				//reduced_EdgeTail.push_back(0);
				//reduced_EdgeHead.push_back(0);
				//reduced_EdgeX.push_back(0);

				//for (int i = DepotArcs + DirectArcs + ChargingArcs + ICEVDepotArcs + 1; i < temp_EdgeTail.size(); i++)
				//{
				//	reduced_EdgeTail.push_back(temp_EdgeTail[i]);
				//	reduced_EdgeHead.push_back(temp_EdgeHead[i]);
				//	reduced_EdgeX.push_back(temp_EdgeX[i]);
				//	/*printf("\n%d\t%d\t%f", temp_EdgeTail[i], temp_EdgeHead[i], temp_EdgeX[i]);*/
				//}
				//printf("\nagg time ICEV");
				//solver.AgressiveSetSearch(reduced_EdgeTail.size(), reduced_EdgeTail.data(), reduced_EdgeHead.data(), reduced_EdgeX.data(), MyCutsCMP, false, true);
				solver.AgressiveSetSearch(EdgeTail_ICEV.size(), EdgeTail_ICEV.data(), EdgeHead_ICEV.data(), EdgeX_ICEV.data(), MyCutsCMP, false, true);

			}
			end = solver.Cplex.getTime();
			solver.Stats->EnumerationTime += end - start;
		}//, EnergyCuts);




		//std::cout << "\n\nNumber of violated cuts: " << MyCutsCMP->Size << "\n\nCuts are:";
		if (MyCutsCMP->Size != 0)
		{
			startAddingCuts = solver.Cplex.getTime();
			//consName = "";
			CutsAdded = true;
			for (IloInt cons = 0; cons < MyCutsCMP->Size; cons++)
			{
				//	std::cout << "\n";
				//	for (IloInt j = 1; j <= MyCutsCMP->CPL[c]->IntListSize; j++)
				//	{

				//		std::cout << MyCutsCMP->CPL[c]->IntList[j] << " ";
				//	}

				//	std::cout << " with RHS = " << MyCutsCMP->CPL[c]->RHS<<"\n\nList is: ";

				CutList.clear();
				ExtList.clear();
				if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_CAP)
				{
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->CapCuts++;
						solver.Stats->numberOfCutsAtRootNode->CapCuts++;
					}
					else solver.Stats->totalNumberOfCuts->CapCuts++;
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j] == n ? 0 : MyCutsCMP->CPL[cons]->IntList[j]);
					}
					expr.clear();

					//for (IloInt i = 0; i < CutList.size(); i++)
					//{
					//	std::cout << CutList[i]<<" ";
					//}


					for (IloInt i = 0; i < CutList.size(); i++)
					{
						for (IloInt j = 0; j < CutList.size(); j++)
						{
							if (i < j)
							{
								expr += solver.x[CutList[i]][CutList[j]][0] + solver.x[CutList[i]][CutList[j]][1] + solver.x[CutList[j]][CutList[i]][0] + solver.x[CutList[j]][CutList[i]][1];
								for (IloInt r = 0; r < solver.R[CutList[i]][CutList[j]].size(); r++)
								{
									expr += solver.y[CutList[i]][CutList[j]][r] + solver.y[CutList[j]][CutList[i]][r];
								}
							}
						}
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//consName = "Cap_cut";
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

					//printf("We add cut CAP");
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_CT_MSTAR) {
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->MultiStar++;
						solver.Stats->numberOfCutsAtRootNode->MultiStar++;
					}
					else solver.Stats->totalNumberOfCuts->MultiStar++;
					std::vector<int> NList;
					std::vector<int> TList;
					std::vector<int> CList;
					NList.reserve(MyCutsCMP->CPL[cons]->IntListSize);
					TList.reserve(MyCutsCMP->CPL[cons]->IntListSize);
					CList.reserve(MyCutsCMP->CPL[cons]->IntListSize);
					//Nucleus:
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						NList.push_back(MyCutsCMP->CPL[cons]->IntList[j] == n ? 0 : MyCutsCMP->CPL[cons]->IntList[j]);
					}
					//Satellites
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->ExtListSize; ++j)
					{
						TList.push_back(MyCutsCMP->CPL[cons]->ExtList[j] == n ? 0 : MyCutsCMP->CPL[cons]->ExtList[j]);
					}
					//Connectors
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->CListSize; ++j)
					{
						CList.push_back(MyCutsCMP->CPL[cons]->CList[j] == n ? 0 : MyCutsCMP->CPL[cons]->CList[j]);
					}

					std::cout << "\nInt list: ";
					for (int i = 0; i < NList.size(); i++)
					{
						std::cout << NList[i] << " ";
					}
					std::cout << "\nExt list: ";
					for (int i = 0; i < TList.size(); i++)
					{
						std::cout << TList[i] << " ";
					}
					std::cout << "\nC list: ";
					for (int i = 0; i < CList.size(); i++)
					{
						std::cout << CList[i] << " ";
					}
					//consName = "MTSTAR_cut";
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_PATH)
				{
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->NoChargePath++;
						solver.Stats->numberOfCutsAtRootNode->NoChargePath++;
					}
					else solver.Stats->totalNumberOfCuts->NoChargePath++;
					//consName = "NCP -> ";
					//std::cout << "\nList Tail: ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//std::cout << ListTail[j-1] << " ";
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
					}
					//std::cout << "\nPath: ";
					//for (int j = 0; j < CutList.size(); ++j)
					//{
					//	//std::cout << CutList[j] << " ";
					//}
					//
					expr.clear();

					for (IloInt i = 0; i < CutList.size() - 1; i++)
					{
						expr += solver.x[CutList[i]][CutList[i + 1]][0] + solver.x[CutList[i + 1]][CutList[i]][0];
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					expr.clear();
					expr2.clear();
					for (int i = 0; i < CutList.size() - 1; i++)
					{
						for (int j = i; j < CutList.size(); j++)
						{
							expr += solver.x[CutList[i]][CutList[j]][0];
							expr2 += solver.x[CutList[j]][CutList[i]][0];
						}
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					add(expr2 <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);

					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_NO_CHARGE_SET)
				{
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->NoChargeSet++;
						solver.Stats->numberOfCutsAtRootNode->NoChargeSet++;
					}
					else solver.Stats->totalNumberOfCuts->NoChargeSet++;

					//consName = "NCS -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j] == n ? 0 : MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";

					}
					//std::cout << "\nSet: ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						//std::cout << CutList[j - 1] << " ";
					}
					expr.clear();

					for (IloInt i = 0; i < CutList.size(); i++)
					{
						for (IloInt j = i + 1; j < CutList.size(); j++)
						{
							if (i < j)
							{
								expr += solver.x[CutList[i]][CutList[j]][0] + solver.x[CutList[j]][CutList[i]][0];
							}
						}
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

					//printf("We add cut no charge set");
					//Draw graph of current solution:
					//DRAWGRAPH_DrawCVRPSparseGraph(n - 1, f, CutList.data() ,CutList.size(), solver.xCoord.data(), solver.yCoord.data(), EdgeDraw.size(), EdgeTailDraw.data(), EdgeHeadDraw.data(),EdgeTailDraw.data(),EdgeHeadDraw.data(),EdgeDraw.data(),0, EdgeDraw.data(), 0.7, GraphName, FileName);
										//for (int i = n; i <n+f ; i++)
					//{
				}
				//else if (MyCutsCMP->CPL[cons]->CType == CMGR_ENERGY_DEPOT) { //ADD DEPOT CONSTRAINTS


				//for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
				//{
				//	CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j] == n ? 0 : MyCutsCMP->CPL[cons]->IntList[j]);
				//	//std::cout << "\n" << CutList[j - 1];
				//}

				//expr.clear();

				//for (IloInt i = 0; i < CutList.size(); i++)
				//{
				//	for (IloInt j = 0; j < n; j++)
				//	{
				//		if (!std::count(CutList.begin(), CutList.end(), j))
				//		{
				//			if (i == 0 && j == 0) continue;
				//			else if (i == CutList.size() - 1 && j == 0) continue;
				//			else expr += solver.x[CutList[i]][j][0]+ solver.x[j][CutList[i]][0];
				//		}
				//		for (int r = 0; r < solver.R[CutList[i]][j].size(); r++)
				//		{
				//			expr += solver.y[CutList[i]][j][r]+ solver.y[j][CutList[i]][r];
				//		}
				//	}
				//}

				//expr -= (solver.x[0][CutList[0]][0] + solver.x[CutList[0]][0][0] + solver.x[0][CutList[CutList.size() - 1]][0] + solver.x[CutList[CutList.size() - 1]][0][0]);
				//add(expr >= 0);
				//}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_CHARGER_PATH) {
					/*printf("Start Charger Path\n\n");*/

					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->ChargerPath++;
						solver.Stats->numberOfCutsAtRootNode->ChargerPath++;
					}
					else solver.Stats->totalNumberOfCuts->ChargerPath++;
					//consName = "CP -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";

						//printf("%d ", CutList[j - 1]);
					}

					expr.clear();
					LHS = 0;
					//Find arcs to add
					for (int i = 0; i < CutList.size() - 1; i++)
					{
						if (CutList[i] >= n) {
							for (int j = 0; j < n; j++)
							{
								if (std::find(CutList.begin() + 1, CutList.end() - 1, j) == CutList.end() - 1)
								{
									for (int r = 0; r < solver.R[CutList[i + 1]][j].size(); r++)
									{
										if (solver.R[CutList[i + 1]][j][r] == CutList[i])
										{
											expr += solver.y[j][CutList[i + 1]][r] + solver.y[CutList[i + 1]][j][r];
										}
									}
								}

							}
							LHS++;
						}
						else if (CutList[i + 1] >= n)
						{
							for (int j = 0; j < n; j++)
							{
								if (std::find(CutList.begin() + 1, CutList.end() - 1, j) == CutList.end() - 1)
								{
									for (int r = 0; r < solver.R[CutList[i]][j].size(); r++)
									{
										if (solver.R[CutList[i]][j][r] == CutList[i + 1])
										{
											expr += solver.y[CutList[i]][j][r] + solver.y[j][CutList[i]][r];
										}
									}
								}

							}
							LHS++;
						}
						else {
							expr += solver.x[CutList[i]][CutList[i + 1]][0] + solver.x[CutList[i + 1]][CutList[i]][0];
							LHS++;
						}
					}
					add(expr <= LHS-1);
					//solver.infeas.add(expr <= LHS-1);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

					//printf("\nEnd Charger Path");
					//printf("We add cut Charger path");
				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_CHARGERS) {
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->FixedChargers++;
						solver.Stats->numberOfCutsAtRootNode->FixedChargers++;
					}
					else solver.Stats->totalNumberOfCuts->FixedChargers++;
					//consName = "FC -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
						//printf("%d ", CutList[j - 1]);
					}

					expr.clear();

					//Add arcs connecting S:
					for (int i = 1; i < CutList.size() - 1; i++)
						for (int j = 1; j < CutList.size() - 1; j++)
							expr += 3 * solver.x[CutList[i]][CutList[j]][0] + 3 * solver.x[CutList[i]][CutList[j]][1];

					//Add connections from r1 to S
					if (CutList[0] > 0) {
						for (int i = 0; i < n; i++)
							for (int j = 1; j < CutList.size() - 1; j++)
								for (int r = 0; r < solver.R[i][CutList[j]].size(); r++)
									if (CutList[0] == solver.R[i][CutList[j]][r]) expr += solver.y[i][CutList[j]][r];
					}
					else
					{
						for (int j = 1; j < CutList.size() - 1; j++)
							expr += solver.x[0][CutList[j]][0];
					}

					//Add connections from S to r2
					if (CutList.back() > 0) {
						for (int i = 0; i < n; i++)
							for (int j = 1; j < CutList.size() - 1; j++)
								for (int r = 0; r < solver.R[i][CutList[j]].size(); r++)
									if (CutList.back() == solver.R[i][CutList[j]][r]) expr += solver.y[CutList[j]][i][r];
					}
					else
					{
						for (int j = 1; j < CutList.size() - 1; j++)
							expr += solver.x[CutList[j]][0][0];
					}


					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_CHARGER_SINGLE) {
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->FixedChargerSingle++;
						solver.Stats->numberOfCutsAtRootNode->FixedChargerSingle++;
					}
					else solver.Stats->totalNumberOfCuts->FixedChargerSingle++;
					//consName = "FCS -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
					}

					//Add arcs connecting S
					expr.clear();
					for (int i = 1; i < CutList.size(); i++)
						for (int j = 1; j < CutList.size(); j++)
							expr += 2 * solver.x[CutList[i]][CutList[j]][0] + 2 * solver.x[CutList[i]][CutList[j]][1];


					if (CutList[0] > 0) { //If the charging option is an external charger
						for (int i = 0; i < n; i++)
							for (int j = 1; j < CutList.size(); j++)
								for (int r = 0; r < solver.R[i][CutList[j]].size(); r++)
									if (solver.R[i][CutList[j]][r] == CutList[0]) expr += solver.y[i][CutList[j]][r];
					}
					else { //If the charging option is the depot
						for (int i = 1; i < CutList.size(); i++)
						{
							expr += solver.x[0][CutList[i]][0];
						}
					}

					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_EDGES) { //ADD FIXED EDGES CONSTRAINTS
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->FixedEdges++;
						solver.Stats->numberOfCutsAtRootNode->FixedEdges++;
					}
					else solver.Stats->totalNumberOfCuts->FixedEdges++;
					//std::cout << "\nList Tail: ";
						//consName = "FE -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";

						//std::cout << CutList[j-1] << " ";
					}
					//std::cout << "\nPath: ";
					//for (int j = 0; j < List.size(); ++j)
					//{
					//	std::cout << List[j] << " ";
					//}
					double r1_dist = solver.c[CutList[0]][CutList[1]];
					double r2_dist = solver.c[CutList[CutList.size() - 2]][CutList[CutList.size() - 1]];

					expr.clear();
					//Add all arcs directly connecting S
					for (IloInt i = 1; i < CutList.size() - 1; i++)
						for (IloInt j = 1; j < CutList.size() - 1; j++)
							expr += solver.x[CutList[i]][CutList[j]][0];


					for (int i = 0; i < n; i++)
					{

						if (std::find(CutList.begin() + 1, CutList.end() - 1, i) != CutList.end() - 1) continue;
						if (solver.c[i][CutList[1]] + solver.MinCharDist[i] >= r1_dist) expr += solver.x[i][CutList[1]][0];// + solver.x[CutList[1]][i][0];
							if (solver.c[CutList[CutList.size() - 2]][i] + solver.MinCharDist[i] >= r2_dist) expr += solver.x[CutList[CutList.size() - 2]][i][0];// + solver.x[i][CutList[CutList.size() - 2]][0];

						for (int r = 0; r < solver.R[i][CutList[1]].size(); r++)
						{
							if (solver.c[solver.R[i][CutList[1]][r]][CutList[1]] >= r1_dist) expr += solver.y[i][CutList[1]][r];// + solver.y[CutList[1]][i][r];
						}
						for (int r = 0; r < solver.R[CutList[CutList.size() - 2]][i].size(); r++)
						{
							if (solver.c[CutList[CutList.size() - 2]][solver.R[CutList[CutList.size() - 2]][i][r]] >= r2_dist) expr += solver.y[CutList[CutList.size() - 2]][i][r];// + solver.y[i][CutList[CutList.size() - 2]][r];
						}
					}

					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_FIXED_SINGLEEDGE) { //ADD FIXED SINGLE EDGE CONSTRAINTS


					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->FixedSingleEdge++;
						solver.Stats->numberOfCutsAtRootNode->FixedSingleEdge++;
					}
					else solver.Stats->totalNumberOfCuts->FixedSingleEdge++;
					//consName = "FES -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";

						//std::cout << ListTail[j-1] << " ";
					}
					//std::cout << "\nPath: ";
					//for (int j = 0; j < List.size(); ++j)
					//{
					//	std::cout << List[j] << " ";
					//}
					double r1_dist = solver.c[CutList[0]][CutList[1]];
					expr.clear();
					for (IloInt i = 1; i < CutList.size(); i++)
					{
						for (IloInt j = 1; j < CutList.size(); j++)
						{
							/*if (j <= i)continue;*/
							expr += solver.x[CutList[i]][CutList[j]][0];
						}


					}



					for (int i = 0; i < n; i++)
					{
						if (std::find(CutList.begin() + 1, CutList.end(), i) != CutList.end()) continue;
						if (solver.c[i][CutList[1]] + solver.MinCharDist[i] >= r1_dist) expr += solver.x[i][CutList[1]][0];// +solver.x[CutList[1]][i][0];
						for (int r = 0; r < solver.R[i][CutList[1]].size(); r++)
						{
							if (solver.c[solver.R[i][CutList[1]][r]][CutList[1]] >= r1_dist) expr += solver.y[i][CutList[1]][r];// +solver.y[CutList[1]][i][r];
						}
					}


					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());



				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_INFEASIBLE_ICEV_PATH)
				{
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->TimeInfeasibleICEVPath++;
						solver.Stats->numberOfCutsAtRootNode->TimeInfeasibleICEVPath++;
					}
					else solver.Stats->totalNumberOfCuts->TimeInfeasibleICEVPath++;
					//std::cout << "\nList Tail: ";
					//consName = "Time_ICEV -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
						//std::cout << ListTail[j-1] << " ";
					}
					//std::cout << "\nPath: ";
					//for (int j = 0; j < CutList.size(); ++j)
					//{
					//	std::cout << CutList[j] << " ";
					//}
					expr.clear();

					for (IloInt i = 0; i < CutList.size() - 1; i++)
					{
						
						expr += solver.x[CutList[i]][CutList[i+1]][0] + solver.x[CutList[i]][CutList[i+1]][1] + solver.x[CutList[i + 1]][CutList[i]][0] + solver.x[CutList[i + 1]][CutList[i]][1];
						for (int r = 0; r < solver.R[CutList[i]][CutList[i+1]].size(); r++) expr += solver.y[CutList[i]][CutList[i+1]][r]+ solver.y[CutList[i+1]][CutList[i]][r];


					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);

					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_TIME_INFEASIBLE_SET)
				{

					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->TimeInfeasibleSet++;
						solver.Stats->numberOfCutsAtRootNode->TimeInfeasibleSet++;
					}
					else solver.Stats->totalNumberOfCuts->TimeInfeasibleSet++;
					//consName = "Time_set -> ";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
					}

					//std::cout << "\nSet: ";
					//for (int j = 1; j <= MyCutsCMP->CPL[c]->IntListSize; ++j)
					//{
					//	std::cout << CutList[j - 1] << " ";
					//}
					expr.clear();

					for (IloInt i = 0; i < CutList.size(); i++)
					{
						for (IloInt j = 0; j < CutList.size(); j++)
						{
							if (i < j)
							{
								expr += solver.x[CutList[i]][CutList[j]][0] + solver.x[CutList[j]][CutList[i]][0]+ solver.x[CutList[i]][CutList[j]][1] + solver.x[CutList[j]][CutList[i]][1];
								for (int r = 0; r < solver.R[CutList[i]][CutList[j]].size(); r++) expr += solver.y[CutList[i]][CutList[j]][r] + solver.y[CutList[j]][CutList[i]][r];
							}
						}
					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_INFEASIBLE_EV_PATH) {
					

					//for (int i = 1; i < EdgeHead_Charger.size(); i++)
					//{
					//	printf("\n%d\t%d\t%f", EdgeTail_Charger[i], EdgeHead_Charger[i], EdgeX_Charger[i]);
					//}
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->TimeInfeasibleEVPath++;
						solver.Stats->numberOfCutsAtRootNode->TimeInfeasibleEVPath++;
					}
					else solver.Stats->totalNumberOfCuts->TimeInfeasibleEVPath++;
					//consName = "Time_EV -> ";
					//printf("\nSet:\n");
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";
						//printf("%d ", CutList[j - 1]);
						//std::cout << MyCutsCMP->CPL[cons]->IntList[j] << "\n";
					}

					expr.clear();
					LHS = 0;
					//Find arcs to add
					for (int i = 0; i < CutList.size() - 1; i++)
					{
						if (CutList[i] >= n && i == 0)
						{
							for (int j = 0; j < n; j++)
							{
								for (int r = 0; r < solver.R[j][CutList[i + 1]].size(); r++)
								{
									if (solver.R[j][CutList[i + 1]][r] == CutList[i]) expr += solver.y[j][CutList[i + 1]][r] + solver.y[CutList[i + 1]][j][r];
								}
							}
							LHS++;
						}
						else if (CutList[i + 1] >= n && i + 1 == CutList.size() - 1) {
							for (int j = 0; j < n; j++)
							{
								for (int r = 0; r < solver.R[CutList[i]][j].size(); r++)
								{
									if (solver.R[CutList[i]][j][r] == CutList[i + 1]) expr += solver.y[CutList[i]][j][r] + solver.y[j][CutList[i]][r];
								}
							}
							LHS++;
						}
						else if (CutList[i + 1] >= n)
						{
							for (int r = 0; r < solver.R[CutList[i]][CutList[i + 2]].size(); r++)
							{
								if (solver.R[CutList[i]][CutList[i + 2]][r] == CutList[i + 1]) { expr += solver.y[CutList[i]][CutList[i + 2]][r] + solver.y[CutList[i + 2]][CutList[i]][r]; break; }
							}
							i++;
							LHS++;
						}
						else {
							expr += solver.x[CutList[i]][CutList[i + 1]][0] + solver.x[CutList[i + 1]][CutList[i]][0];
							LHS++;
						}
					}
					add(expr <= LHS - 1);

					solver.cons.add(expr <= LHS - 1);
					
					//std::cout << "Variables: " << std::endl;
					//for (IloExpr::LinearIterator it(expr.getLinearIterator()); it.ok(); ++it) {
					//	std::cout << "  " << it.getVar().getName() << std::endl;
					//}
					//LHS = LHS;
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());

				}
				else if (MyCutsCMP->CPL[cons]->CType == CMGR_TIME_INFEASIBLE_EV_SET)
				{
					if (NNodes == 0) {
						solver.Stats->totalNumberOfCuts->TimeInfeasibleEVSet++;
						solver.Stats->numberOfCutsAtRootNode->TimeInfeasibleEVSet++;
					}
					else solver.Stats->totalNumberOfCuts->TimeInfeasibleEVSet++;
					//consName = "Time_EV_set";
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->IntListSize; ++j)
					{
						CutList.push_back(MyCutsCMP->CPL[cons]->IntList[j]);
						//consName = consName + std::to_string(CutList[j - 1]) + " ";

						//printf("%d ", CutList[j - 1]);
					}
					for (int j = 1; j <= MyCutsCMP->CPL[cons]->ExtListSize; ++j)
					{
						ExtList.push_back(MyCutsCMP->CPL[cons]->ExtList[j]);
						//consName = consName + std::to_string(ExtList[j - 1]) + " ";

						//printf("%d ", CutList[j - 1]);
					}

					//Arcs connecting nodes in S
					expr.clear();
					for (int i = 0; i < CutList.size(); i++)
					{
						for (int j = 0; j < CutList.size(); j++)
						{
							expr += solver.x[CutList[i]][CutList[j]][0];
							for (int r = 0; r < solver.R[CutList[i]][CutList[j]].size(); r++)
							{
								if (std::find(ExtList.begin(), ExtList.end(), solver.R[CutList[i]][CutList[j]][r]) != ExtList.end())
								{
									expr += 2 * solver.y[CutList[i]][CutList[j]][r];
								}
							}
						}
					}

					for (int i = 0; i < n; i++)
					{
						if (std::find(CutList.begin(), CutList.end(), i) == CutList.end())
						{
							for (int j = 0; j < CutList.size(); j++)
							{
								for (int r = 0; r < solver.R[i][CutList[j]].size(); r++)
								{
									if (std::find(ExtList.begin(), ExtList.end(), solver.R[i][CutList[j]][r]) != ExtList.end())
									{
										expr += solver.y[i][CutList[j]][r];
									}
								}
								for (int r = 0; r < solver.R[CutList[j]][i].size(); r++)
								{
									if (std::find(ExtList.begin(), ExtList.end(), solver.R[CutList[j]][i][r]) != ExtList.end())
									{
										expr += solver.y[CutList[j]][i][r];
									}
								}
							}
						}

					}
					add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas.add(expr <= MyCutsCMP->CPL[cons]->RHS);
					//solver.infeas[solver.infeas.getSize() - 1].setName(consName.c_str());
				}
			}


			for (i = 0; i < MyCutsCMP->Size; i++)
			{
				CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
			}

			MyCutsCMP->Size = 0;
			end = solver.Cplex.getTime();
			solver.Stats->AddCutsTime += end - startAddingCuts;
			//if (solver.DrawCounter % 10 == 0)
			//{
			//	std::string strFilename = "Grafer/Tid/Elbils_løsning_n21_f9_tid_" + std::to_string(solver.DrawCounter) + ".tex";
			//	char* FileName;
			//	FileName = &strFilename[0];
			//	std::string strGraphname = "Elbils loesning (tidsforbrug) eksempel " + std::to_string(solver.DrawCounter);
			//	char* GraphName = &strGraphname[0];
			//	DRAWGRAPH_DrawCVRPSparseGraph(n - 1, f, EdgeTailDraw.data(), 0, solver.xCoord.data(), solver.yCoord.data(), EdgeDraw.size(), EdgeTailDraw.data(), EdgeHeadDraw.data(), EnergyComps.data(), TravelTime.data(), d, s, EdgeDraw.data(), 0.7, GraphName, FileName, true);
			//	strFilename = "Grafer/Energi/Elbils_løsning_n21_f9_energi_" + std::to_string(solver.DrawCounter) + ".tex";
			//	FileName = &strFilename[0];
			//	strGraphname = "Elbils loesning (energiforbrug) eksempel " + std::to_string(solver.DrawCounter);
			//	GraphName = &strGraphname[0];
			//	DRAWGRAPH_DrawCVRPSparseGraph(n - 1, f, EdgeTailDraw.data(), 0, solver.xCoord.data(), solver.yCoord.data(), EdgeDraw.size(), EdgeTailDraw.data(), EdgeHeadDraw.data(), EnergyComps.data(), TravelTime.data(), d, s, EdgeDraw.data(), 0.7, GraphName, FileName, false);
			//}
		}
	}

	//delete[] Label;
	//delete[] NList;
	//delete[] TList;
	//delete[] CList;
	delete[] d;
	delete[] s;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			xTemp[i][j].end();
			yTemp[i][j].end();
		}
		xTemp[i].end();
		yTemp[i].end();
	}
	xTemp.end();
	yTemp.end();


	if (getNnodes() == 0) {
		solver.Stats->LBAtRootNode = getBestObjValue();
		solver.Stats->timeRootNode = solver.Cplex.getTime();
	}

	if (solver.Stats->treeDepth < getCurrentNodeDepth()) solver.Stats->treeDepth = getCurrentNodeDepth();

	if (!CutsAdded) {
		if (solver.Stats->treeFeasDepth == -1) solver.Stats->treeFeasDepth = getCurrentNodeDepth();
		if (solver.Stats->currentUB > getObjValue()) {
			solver.Stats->currentUB = getObjValue();
			solver.Stats->treeUBDepth = getCurrentNodeDepth();
		}
	}
	endTot = solver.Cplex.getTime();
	solver.Stats->TotalSepTime += endTot - startTot;
}

//TSP cuts LAZY:
ILOLAZYCONSTRAINTCALLBACK1(LazyTSPCuts, MFGVRP_Solver&, solver) {
	IloEnv env = solver.getEnv();
	IloExpr expr(env);
	


	double ObjVal = getObjValue();
	double B = solver.getBatteryCap();
	double T = solver.getMaxTime();

	int NoOFCustomers = solver.getNoTSPCustomers();
	int* TSPCustomers = new int[NoOFCustomers];
	for (int i = 0; i < NoOFCustomers; i++) TSPCustomers[i] = solver.getTSPCutomers(i);
	int idx = TSPCustomers[0];
	int cnt = 0;
	std::vector<int> SubTour;
	if (getBestObjValue()<=solver.TSPCutOff+0.03)
	{
		int f = solver.getChargers();
		int n = solver.getNumCust();
		IloNum2DMatrix xVals(env, n + f + 1);


		for (IloInt i = 0; i < n + f + 1; i++) {
			xVals[i] = IloNumArray(env, n + f + 1);
			getValues(xVals[i], solver.xTSP[i]);
		}
		while (true)
		{
			for (int i = 0; i < NoOFCustomers; i++)
			{
				if (xVals[idx][TSPCustomers[i]] >= 0.99) {

					idx = TSPCustomers[i];
					SubTour.push_back(idx);
					if (idx != TSPCustomers[0]) { ++cnt; }
					break;
				}
			}
			if (idx == TSPCustomers[0]) {
				if (cnt < NoOFCustomers - 1)
				{
					expr.clear();
					for (int i = 0; i < SubTour.size(); i++) for (int j = 0; j < SubTour.size(); j++) expr += solver.xTSP[SubTour[i]][SubTour[j]];
					add(expr <= cnt);
					break;
				}
				else break;
			}
		}
		for (int i = 0; i < n + f + 1; i++)
			xVals[i].end();
		xVals.end();
	}


	delete[] TSPCustomers;

}
//Define constructor:
MFGVRP_Solver::MFGVRP_Solver() {
	//Initializing CPLEX models:
	Model = IloModel(env);

	Cplex = IloCplex(Model);
	//Cplex.setParam(IloCplex::WorkDir, "c:/cplex/");
	//Cplex.setParam(IloCplex::NodeFileInd, 3);
	//Cplex.setParam(IloCplex::Param::MIP::Tolerances::UpperCutoff, 309.5);
	//Cplex.setParam(IloCplex::WorkMem, 500.0);
	//Cplex.setParam(IloCplex::SolnPoolCapacity, 1);
	//Cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 0); //Only process root node
	//Cplex.setParam(IloCplex::Threads, 1);
	//Cplex.setParam(IloCplex::ParallelMode, 1);
	//Cplex.setParam(IloCplex::EpRHS, 1e-9);
	//Cplex.setParam(IloCplex::EpAGap, 0.0);
	Cplex.setParam(IloCplex::VarSel, 4);
	Cplex.setParam(IloCplex::MemoryEmphasis, true);
	//Variables
	x = IloVar3DMatrix(env);
	y = IloVar3DMatrix(env);
	//SOC = IloNumVarArray(env);
	//e = IloVar2DMatrix(env);
	//time = IloNumVarArray(env);

	//Parameters:
	c = IloNum2DMatrix(env); //Cost for driving between customer i and j with charging mode p
	d = IloNumArray(env); //Customer demands
	s = IloNumArray(env); //Contains service times all customers
	DrawCounter = 0;

	lazy = IloRangeArray(env);
	user = IloRangeArray(env);
	//infeas = IloRangeArray(env);
	//Constraints:
	cons = IloRangeArray(env); //Holds all initial constraints
	//infeas = IloRangeArray(env);
	CapCuts = IloRangeArray(env); //Holds all capacity cuts
	//Expressions
	expr = IloExpr(env); //New expression 1
	expr2 = IloExpr(env); //New expression 2
	//infeas = IloRangeArray(env);

	//Fractional cuts:
	EpsForIntegrality = 0.0001;
	MaxViolation = 0;

	//Statistics:
	Stats = new testStats;
}

//---------------------------------------------------------------------------------------------------------------------------//
MFGVRP_Solver::~MFGVRP_Solver() {
	if (!HasCleanedUp) {
		//Freeing memory occupied by variables:
		if (x.getImpl()) {
			if (n > 0) for (IloInt i = 0; i < x.getSize(); i++) x[i].end();
			x.end();
		}
		if (y.getImpl()) {
			for (int i = 0; i < y.getSize(); i++)
			{
				for (int j = 0; j < y[i].getSize(); j++)
				{
					y[i][j].end();
				}
				y[i].end();
			}
			y.end();
		}

		if (c.getImpl()) {
			for (int i = 0; i < c.getSize(); i++)
			{
				c[i].end();
			}
			c.end();
		}

		if (d.getImpl()) d.end();
		if (s.getImpl()) s.end();
		if (SOC.getImpl()) SOC.end();
		if (time.getImpl()) time.end();
		if (e.getImpl())
		{
			for (int i = 0; i < e.getSize(); i++)
			{
				e[i].end();
			}
			e.end();

		}
		if (env.getImpl()) env.end();
		//if (Cplex.getImpl()) Cplex.end();
		//if (cons.getImpl()) cons.end();
		//if (Model.getImpl()) Model.end();
		//
		//if (TSPCplex.getImpl()) TSPModel.end();
		//if (TSPCons.getImpl()) TSPCons.end();
		//if (TSPModel.getImpl()) TSPModel.end();

		

		//Freeing memory occupied by statistics:
		FreeMemoryFromTestStats(Stats);
		delete Stats;
		HasCleanedUp = true;
	}
}

//---------------------------------------------------------------------------------------------------------------------------//
//Define load data
void MFGVRP_Solver::LoadData(const std::string &FileName, bool RunAllCuts) {
	RunCuts = RunAllCuts;
	double SomeNumber; //Used to load in data to cplex arrays
	int x_in, y_in;
	std::vector<std::pair<int, int>> GZ;
	std::ifstream in(FileName, std::ios::app);
	if (!in) {
		std::cout << "Cannot open file.\n";
	}

	in >> n; //Number of customers 
	in >> f; //Number of charging stations
	in >> E; //Number of electric vehicles
	in >> C; //Number of conventional vehicles
	in >> G; //Number of green zones
	in >> T; //Max tour duration
	in >> B; //Battery capacity 
	in >> g; //consumption rate
	in >> r; //Refueling rate
	in >> Q; //Load capacity


	//Load in demands:
	for (IloInt i = 0; i < n + 1; i++)
	{
		in >> SomeNumber;
		d.add(SomeNumber);
	}

	//Load in service times:
	for (IloInt i = 0; i < n + 1; i++)
	{
		in >> SomeNumber;
		s.add(SomeNumber);
	}

	//std::cout << "\n\nC_DIR_ij: ";
	//for (IloInt i = 0; i < n + f + 1; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n + 1; j++)
	//	{
	//		in >> SomeNumber;
	//		std::cout << SomeNumber << " ";
	//	}
	//}

	//Load coordinates
	for (IloInt i = 0; i < n + f + 1; i++)
	{
		in >> SomeNumber;
		xCoord.push_back(SomeNumber);
		in >> SomeNumber;
		yCoord.push_back(SomeNumber);
	}

	//std::cout << "\n\nC_DIR_ij: ";
	//Load direct distance matrix
	c = IloNum2DMatrix(env, n + f + 1);
	minTravelTime = new int[n + 1];
	for (IloInt i = 0; i < n + f + 1; i++)
	{
		c[i] = IloNumArray(env, n + f + 1); //Add dimension
		//std::cout << "\n";
		for (IloInt j = 0; j < n + f + 1; j++)
		{
			c[i][j] = calc_distances(xCoord[i] - xCoord[j], yCoord[i] - yCoord[j]);
			//std::cout << c[i][j] << " ";
		}
		if (i<=n)
		{
			minTravelTime[i] = calcMinTravelTime(i);
		}
	}


	//Save customer points:
	for (IloInt i = 1; i < n + 1; i++) CustomersPoints.push_back(std::make_pair(xCoord[i], yCoord[i]));

	//Load green zone points:
	for (int i = 0; i < G; i++)
	{
		in >> SomeNumber;
		GZ.clear();
		for (int j = 0; j < SomeNumber; j++)
		{
			in >> x_in;
			in >> y_in;
			GZ.push_back(std::make_pair(x_in, y_in));
		}
		GZPoints.push_back(GZ);
	}

	//std::cout << "\nDistance to charging stations: ";
	////Load distances to charging nodes
	//for (IloInt i = 0; i < f; i++)
	//{
	//	u.add(IloNumArray(env, n + 1));
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n + 1; j++)
	//	{
	//		in >> c[i][j];
	//		std::cout << c[i][j] << " ";
	//	}
	//}
	//std::cout << "\nc_REC_ij: ";
	////Find minimum charging paths
	//for (IloInt i = 0; i < n + 1; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n + 1; j++)
	//	{
	//		minDist = 1000;
	//		for (IloInt p = 0; p < f; p++)
	//		{
	//			if (minDist > u[p][i] + u[p][j])
	//			{
	//				minDist = u[p][i] + u[p][j];
	//				MinDistCharg= p;
	//			}
	//		}
	//		if (i != j) {c[i][j][1] = minDist;}
	//		else { c[i][j][1] = 999;}
	//		std::cout << c[i][j][1] << " ";
	//		
	//	}
	//}

	in.close();
	SetupGreenZone();
	int numY = NonDominatedChargingPaths();

	//Setup variables:
	y = IloVar3DMatrix(env, n + 1); //Setting 1st dimensions for y variables
	x = IloVar3DMatrix(env, n + 1); //Setting 1st dimensions for x variables
	xDummy = IloNumVarArray(env, (n + 1) * (n + 1) *2, 0, 1, ILOINT);
	yDummy = IloNumVarArray(env, numY, 0, 1, ILOINT);

	for (IloInt i = 0; i < n + 1; i++)
	{
		x[i] = IloVar2DMatrix(env, n + 1); //Setting 2nd dimensions for x variables
		y[i] = IloVar2DMatrix(env, n + 1); //Setting 2nd dimensions for y variables
		for (IloInt j = 0; j < n + 1; j++)
		{
			x[i][j] = IloNumVarArray(env, 2, 0, 1, ILOINT); 
			y[i][j] = IloNumVarArray(env, R[i][j].size(), 0, 1, ILOINT); 
			//xSol1[i][j] = IloNumArray(env, 2);
			for (IloInt k = 0; k < 2; k++)
			{
				std::string xName = "x_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
				x[i][j][k].setName(xName.c_str());
			}
			for (IloInt r = 0; r < R[i][j].size(); r++)
			{
				std::string yName = "y_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(r);
				y[i][j][r].setName(yName.c_str());
			}
		}
	}

	//Initiate vectors for identifying paths:
	AddedToPath = new bool[n + f + 1];
	LevelNode = new int[n + f + 1];
	LevelIdx = new int[n + f + 1];
	Slack = new double[n + f + 1];
	SlackSet = new double[n + f + 1];
	Esum = new double[n + f + 1];
	Tsum = new double[n + f + 1];
	rhs=IloIntArray (env, n * 2 + 2*f + 2);

	//Set number of cap cuts:
	if(n>50)MaxNoOfCuts = 10; else MaxNoOfCuts = 5;

	Stats->instanceName = FileName;
	Stats->n = n;
	Stats->r = f;

	BuildModels();
}

//---------------------------------------------------------------------------------------------------------------------------//
void MFGVRP_Solver::BuildModels() {

	int count = 0;
	int consindex = 0;
	std::string consName = "";
	//PROBLEM SETUP
	//Set tolerance:
	Cplex.setParam(IloCplex::EpGap, 0);

	//Build objective function: (1)
	expr.clear();
	for (IloInt i = 0; i < n + 1; i++)
	{
		for (IloInt j = 0; j < n + 1; j++)
		{
			for (IloInt k = 0; k < 2; k++)
			{
				expr+= c[i][j]* x[i][j][k];
			}
			for (IloInt r = 0; r < R[i][j].size(); r++)
			{
				expr+= (c[R[i][j][r]][i] + c[R[i][j][r]][j])*y[i][j][r];
			}

		}
	}
	Obj = IloObjective(env);
	Obj = IloMinimize(env,expr);
	Model.add(Obj);
	//Cons 2: All customers must be visited by either ICEV or EV (2)
	for (IloInt j = 1; j < n + 1; j++)
	{
		expr.clear();
		//expr2.clear();
		for (IloInt i = 0; i < n + 1; i++)
		{
			for (IloInt k = 0; k < 2; k++)
			{
				expr += x[i][j][k];
				//expr2 += x[j][i][k];
			}
			for (IloInt r = 0; r < R[i][j].size(); r++)
			{
				expr += y[i][j][r];
				//expr2 += y[j][i][r];
			}
		}
		cons.add(expr == 1);
		//infeas.add(expr == 1);
		consName = "Must_enter_customer_" + std::to_string(j);
		cons[count].setName(consName.c_str());
		count++;
		//cons.add(expr2 == 1);
		//consName = "Must_leave_customer_" + std::to_string(j);
		//cons[count].setName(consName.c_str());
		//count++;
	}

	//LB on Vehicles:
	
	expr.clear();
	expr2.clear();
	for (IloInt j = 1; j < n + 1; j++)//(5)
	{
		expr += x[0][j][0];// +x[0][j][1];
		expr2 += x[0][j][1];// +x[0][j][1];
		for (IloInt r = 0; r < R[0][j].size(); r++)
		{
			expr += y[0][j][r];
		}

	}

	cons.add(expr == E);
	//infeas.add(expr >= std::ceil(DemandSum / Q));
	consName = "LB_on_vehicles";
	cons[count].setName(consName.c_str());
	count++;
	cons.add(expr2 == C);
	//infeas.add(expr >= std::ceil(DemandSum / Q));
	consName = "LB_on_vehicles";
	cons[count].setName(consName.c_str());
	count++;

	////DUMMY CONSTRAINT TO LEAVE OUT ICEV'S
	//for (IloInt j = 1; j < n + 1; j++)
	//{
	//	expr.clear();
	//	for (IloInt i = 0; i < n + 1; i++)
	//	{
	//		expr += x[i][j][1];
	//	}
	//	cons.add(expr <= 0);
	//	//infeas.add(expr <= C);
	//	count++;
	//}

	////DUMMY CONSTRAINT TO LEAVE OUT ICEV'S
	//for (IloInt j = 1; j < n + 1; j++)
	//{
	//	expr.clear();
	//	for (IloInt i = 0; i < n + 1; i++)
	//	{
	//		expr += x[j][i][1];
	//	}
	//	cons.add(expr <= C);
	//	infeas.add(expr <= C);
	//	count++;
	//}
	//Don't drive depot to depot and from/to the same customer through a charger:

	for (int i = 0; i < n + 1; i++)
	{
		for (int k = 0; k < 2; k++) x[i][i][k].setUB(0);
		for (int r = 0; r < R[i][i].size(); r++) y[i][i][r].setUB(0);
		for (int r = 0; r < R[0][i].size(); r++) if (R[0][i][r] == n + 1) { y[0][i][r].setUB(0); break; }
		for (int r = 0; r < R[i][0].size(); r++) if (R[i][0][r] == n + 1) { y[i][0][r].setUB(0); break; }
	}

	//cons.add(x[0][0][0] + x[0][0][0] == 0);
	//count++;

	// Cons 3+4: If a vehicle enter a node, it must leave it again 
	for (IloInt j = 0; j < n + 1; j++)
	{
		expr.clear(); //(3)
		for (IloInt i = 0; i < n + 1; i++)
		{
			expr += x[i][j][0] - x[j][i][0];
			for (IloInt r = 0; r < R[i][j].size(); r++)
			{
				expr += y[i][j][r] - y[j][i][r];
			}
		}
		cons.add(expr == 0);
		//infeas.add(expr == 0);
		consName = "Must_leave_customer_" + std::to_string(j) + "_if_visited_by_EV";
		cons[count].setName(consName.c_str());
		count++;
		expr.clear(); //(4)
		for (IloInt i = 0; i < n + 1; i++)
		{
			expr += x[i][j][1] - x[j][i][1];
		}
		cons.add(expr == 0);
		//infeas.add(expr == 0);
		consName = "Must_leave_customer_" + std::to_string(j) + "_if_visited_by_ICEV";
		cons[count].setName(consName.c_str());
		count++;

	}
	//Cons 5: We cannot use more EVs than available:
	for (IloInt j = 1; j < n + 1; j++)//(5)
	{
		expr += x[0][j][0];
		for (IloInt r = 0; r < R[0][j].size(); r++)
		{
			expr += y[0][j][r];
		}

	}
	cons.add(expr <= E);
	//infeas.add(expr <= E);
	consName = "We_cannot_use more_EVs_than_available";
	cons[count].setName(consName.c_str());
	count++;
	expr.clear();
	////Cons 6:  We cannot use more ICEVs than available:
	//for (IloInt j = 1; j < n + 1; j++)//(6)
	//{
	//	expr += x[0][j][1];
	//}
	//cons.add(expr == C);
	////infeas.add(expr <= C);
	//consName = "We_cannot_use more_ICEVs_than_available";
	//cons[count].setName(consName.c_str());
	//count++;


	//Setup dummy variables for solution fetching 
	int cnt = 0;
	
	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				cons.add(x[i][j][k] - xDummy[cnt] == 0);
				cnt++;
			}
		}
	}


	cnt = 0;

	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			for (int r = 0; r < R[i][j].size(); r++)
			{
				cons.add(y[i][j][r] - yDummy[cnt] == 0);
				cnt++;
			}
		}
	}
	// Cons 8: All EVs leave depot with full charge (8)
	//SOC[0].setLB(B);

	//// Cons 9: We can only recharge if we use recharging path (9)
	//for (IloInt j = 0; j < n + 1; j++)
	//{
	//	for (IloInt i = 0; i < n + 1; i++)
	//	{
	//		expr.clear();
	//		for (IloInt r = 0; r < R[i][j].size(); r++)
	//		{
	//			expr += y[i][j][r];
	//		}
	//		cons.add(e[i][j] - B * expr <= 0);
	//		consName = "Only_recharge_if_y_" + std::to_string(i) + "_" + std::to_string(j) + "_is_selected";
	//		cons[count].setName(consName.c_str());
	//		count++;
	//	}
	//}

	//// Cons TEST CONSTRAINT: WE NEED SUFFICIENT BATTERY
	//expr.clear();
	//for (IloInt j = 0; j < n + 1; j++)
	//{
	//	expr -= B * x[0][j][0];
	//	for (IloInt i = 0; i < n + 1; i++)
	//	{
	//		expr += g * c[i][j] * x[i][j][0];
	//		for (IloInt r = 0; r < R[i][j].size(); r++)
	//		{
	//			expr += g * (c[R[i][j][r]][i] + c[R[i][j][r]][i]) * y[i][j][r] - e[i][j];
	//		}
	//	}
	//}
	//cons.add(expr <= 0);
	//consName = "Overall_battery_available";
	//cons[count].setName(consName.c_str());
	//count++;


	//// Cons TEST CONSTRAINT: GZ Customers
	//expr.clear();
	//for (IloInt i = 0; i < n + 1; i++)
	//{
	//	
	//	for (IloInt j = 0; j < GreenZoneCustomers.size(); j++)
	//	{
	//		expr += x[i][GreenZoneCustomers[j]][1];
	//	}
	//}
	//cons.add(expr == 0);
	//consName = "GREEN_ZONE_CONSTRAINT";
	//cons[count].setName(consName.c_str());
	//count++;
	//// Cons 10: We cannot recharge more than the battery capacity (10)
	//for (IloInt j = 0; j < n + 1; j++)
	//{
	//	for (IloInt i = 0; i < n + 1; i++)
	//	{
	//		expr.clear();
	//		expr2.clear();
	//		for (IloInt r = 0; r < R[i][j].size(); r++)
	//		{
	//			expr += y[i][j][r];
	//			expr2 += c[R[i][j][r]][i] *y[i][j][r];
	//		}
	//		cons.add(e[i][j] + SOC[i] - B * (1 - expr) - expr2 <= B);
	//		consName = "Max_recharge_for_" + std::to_string(i) + "_" + std::to_string(j);
	//		cons[count].setName(consName.c_str());
	//		count++;
	//	}
	//}

	//// Cons 11: We must regulate SOC at each node as we drive a route: (11)
	//for (IloInt i = 0; i < n + 1; i++)
	//{
	//	for (IloInt j = 1; j < n + 1; j++)
	//	{
	//		//if (i != j) {
	//		for (IloInt r = 0; r < R[i][j].size(); r++)
	//		{
	//			cons.add(SOC[j] - SOC[i] - e[i][j] + (B + g * c[i][j]) * x[i][j][0] + (B + g *(c[R[i][j][r]][i]+ c[R[i][j][r]][j])) * y[i][j][r] <= B);
	//			consName = "Regulate_charge_for_" + std::to_string(i) + "_" + std::to_string(j);
	//			cons[count].setName(consName.c_str());
	//			count++;
	//		}
	//		//}
	//	}
	//}

	//// Cons 12: We must ensure that we dont run out of battery and are able to reach a charging node (12)
	//for (IloInt i = 1; i < n + 1; i++)
	//{
	//	expr.clear();
	//	for (IloInt j = 0; j < n+1; j++)
	//	{
	//		for (IloInt r = 0; r < R[i][j].size(); r++)
	//		{
	//			cons.add(SOC[i] - g * c[R[i][j][r]][i] * y[i][j][r] >= 0);
	//			consName = "Must_reach_charging_station_before_depletion_" + std::to_string(i) + "_" + std::to_string(j);
	//			cons[count].setName(consName.c_str());
	//			count++;
	//		}
	//	}
	//}

	//// Cons 13: We must ensure that we can reach the depot with the available charge (13)
	//for (IloInt i = 1; i < n + 1; i++)
	//{
	//	expr.clear();
	//	for (IloInt r = 0; r < R[i][0].size(); r++)
	//	{
	//		//expr += y[i][0][r];
	//		cons.add(SOC[i] - g * c[i][0] * x[i][0][0] - g * (c[R[i][0][r]][i] + c[R[i][0][r]][0]) * y[i][0][r] + e[i][0] >= 0);
	//		consName = "Enough_battery_to_return_to_depot_from_" + std::to_string(i);
	//		cons[count].setName(consName.c_str());
	//		count++;
	//	}

	//}

	//// Cons 14: Specify arrival time of first customer (14)
	//for (IloInt j = 1; j < n + 1; j++)
	//{
	//	expr.clear();
	//	for (IloInt k = 0; k < 2; k++)
	//	{
	//		expr += x[0][j][k];
	//	}
	//	expr2.clear();
	//	for (IloInt r = 0; r < R[0][j].size(); r++)
	//	{
	//		expr2 += (c[R[0][j][r]][0] + c[R[0][j][r]][j]) * y[0][j][r];
	//	}
	//	cons.add(time[j] - c[0][j] * expr - expr2 - r * e[0][j] >= 0);
	//	consName = "Arrival_time_of_" + std::to_string(j) + "_if_first_customer_in_route";
	//	cons[count].setName(consName.c_str());
	//	count++;
	//}
	////TAKEN OUT AS WE DO NOT NEED EXACT TIME BUT ONLY ENSURE WE DO NOT EXCEED WORKDAY DURATION:
	//// //Cons 10: We must initiate time variable with the trip from the depot (21)
	////for (IloInt j = 1; j < n + 1; j++)
	////{
	////	expr.clear();
	////	for (IloInt k = 0; k < 2; k++)
	////	{
	////		expr += x[0][j][k];
	////	}
	////	cons.add(time[j] + (T-c[0][j][0])*expr + (T - c[0][j][1])*y[0][j] - e[0][j] <= T);
	////	cons[count].setName("Initiate time variable");
	////	count++;
	////}

	//// Cons 15: Regulates time during route (15)
	//for (IloInt i = 1; i < n + 1; i++)
	//{
	//	for (IloInt j = 1; j < n + 1; j++)
	//	{
	//		expr.clear();
	//		for (IloInt k = 0; k < 2; k++)
	//		{
	//			expr += x[i][j][k];
	//		}
	//		expr2.clear();
	//		for (IloInt r = 0; r < R[i][j].size(); r++)
	//		{
	//			expr2 += (T + s[i] + c[R[i][j][r]][i] + c[R[i][j][r]][j]) * y[i][j][r];
	//		}

	//		cons.add(time[i] - time[j] + (T + s[i] + c[i][j]) * expr + expr2 + r * e[i][j] <= T);
	//		consName = "Regulate_time_for_" + std::to_string(i) + "_" + std::to_string(j);
	//		cons[count].setName(consName.c_str());
	//		count++;
	//	}
	//}

	//// Cons 16: We must return to depot before maximum day duration (16)
	//for (IloInt i = 1; i < n + 1; i++)
	//{
	//	expr.clear();
	//	for (IloInt k = 0; k < 2; k++)
	//	{
	//		expr += x[i][0][k];
	//	}
	//	expr2.clear();
	//	for (IloInt r = 0; r < R[i][0].size(); r++)
	//	{
	//		expr2 += (s[i] + c[R[i][0][r]][i] + c[R[i][0][r]][0]) * y[i][0][r];
	//	}
	//	cons.add(time[i] + (s[i] + c[i][0]) * expr + expr2 + r * e[i][0] <= T);
	//	consName = "Return_from_" + std::to_string(i) + "_before_max_time";
	//	cons[count].setName(consName.c_str());
	//	count++;
	//}

	Model.add(cons);
	Cplex.setParam(IloCplex::TiLim, 3600);


	//Cplex.extract(Model);
	//Cplex.exportModel("MY_MODEL.lp");
	Run();

}

//---------------------------------------------------------------------------------------------------------------------------//
void MFGVRP_Solver::Run() {
	Cplex.use(LazyCapCuts(env, *this));
	Cplex.use(FracCuts(env, *this));
	double startTime = Cplex.getTime();
	SetupTSP();
	try
	{
		if (Cplex.solve()) {
			Stats->bestUpperBound = Cplex.getObjValue();
			for (int i = 0; i < n + 1; i++)
			{
				Stats->EVs += Cplex.getValue(x[0][i][0]);
				Stats->ICEVs += Cplex.getValue(x[0][i][1]);
				for (int r = 0; r < R[0][i].size(); r++)
				{
					Stats->EVs += Cplex.getValue(y[0][i][r]);
				}
			}
		}


		
		//infeas.add(cons);



		//infeas.add(x[0][14][0]==1);
		//infeas.add(x[14][6][0] == 1);
		//infeas.add(x[6][22][0] == 1);
		//infeas.add(x[22][19][0] == 1);
		//infeas.add(x[20][12][0] == 1);
		//infeas.add(y[19][20][2]==1);
		//infeas.add(x[12][0][0] == 1);
		//infeas.add(x[0][23][0] == 1);
		//infeas.add(x[23][0][0] == 1);

		//infeas.add(x[0][10][1] == 1);
		//infeas.add(x[10][17][1] == 1);
		//infeas.add(x[17][3][1] == 1);
		//infeas.add(x[3][1][1] == 1);
		//infeas.add(x[1][11][1] == 1);
		//infeas.add(x[11][8][1] == 1);
		//infeas.add(x[8][7][1] == 1);
		//infeas.add(x[7][5][1] == 1);
		//infeas.add(x[5][18][1] == 1);
		//infeas.add(x[18][0][1] == 1);

		//infeas.add(x[0][24][1] == 1);
		//infeas.add(x[24][15][1] == 1);
		//infeas.add(x[15][25][1] == 1);
		//infeas.add(x[25][21][1] == 1);
		//infeas.add(x[21][4][1] == 1);
		//infeas.add(x[4][13][1] == 1);
		//infeas.add(x[13][16][1] == 1);
		//infeas.add(x[16][2][1] == 1);
		//infeas.add(x[2][9][1] == 1);
		//infeas.add(x[9][0][1] == 1);



		
		///*if (infeasCplex.getStatus() == IloAlgorithm::Infeasible){*/
		//	IloModel FeasModel(env); //TSP model
		//	IloCplex FeasCplex(FeasModel); //TSP problem solver
		//	
		//	FeasModel.add(Obj);
		//	FeasModel.add(infeas);

		//	FeasCplex.solve();
		//	std::cout << "\n" << "Cplex stats:\n" << FeasCplex.getStatus() << "\n"<<FeasCplex.getObjValue()<<"\n";
		//	IloNumArray preferences(env);
		//	IloConstraintArray Infeasible(env);
		//	Infeasible.add(infeas);
		//	for (IloInt i = 0; i < Infeasible.getSize(); i++) {
		//		preferences.add(1.0);  // user may wish to assign unique preferences
		//	}

		//	if (FeasCplex.refineConflict(Infeasible, preferences)) {
		//		IloCplex::ConflictStatusArray conflict = FeasCplex.getConflict(Infeasible);
		//		env.getImpl()->useDetailedDisplay(IloTrue);
		//		std::cout << "Conflict :" << std::endl;
		//		for (IloInt i = 0; i < Infeasible.getSize(); i++) {
		//			if (conflict[i] == IloCplex::ConflictMember)
		//				std:: cout << "Proved  : " << Infeasible[i] << std::endl;
		//			if (conflict[i] == IloCplex::ConflictPossibleMember)
		//				std::cout << "Possible: " << Infeasible[i] << std::endl;
		//		}
		//	}
		//	else
		//		std::cout << "Conflict could not be refined" << std::endl;
		//	std::cout << std::endl;
		////}

		Stats->time =  Cplex.getTime() - startTime;
		Stats->treeSize = Cplex.getNnodes();

		if (Cplex.getStatus() == IloAlgorithm::Status::Optimal)
			Stats->isOptimal = true;
		else
			Stats->isOptimal = false;
		if (Cplex.getStatus() == IloAlgorithm::Status::Infeasible)
			Stats->isInfeasible = true;
		else
			Stats->isInfeasible = false;

		Stats->bestLowerBound = Cplex.getBestObjValue();
		if (!Stats->TSPNodes.empty()) {
			auto it = std::minmax_element(Stats->TSPTime.begin(), Stats->TSPTime.end());
			Stats->minTSPtime = *it.first;
			Stats->maxTSPTime = *it.second;
			Stats->avgTSPTime = std::accumulate(Stats->TSPTime.begin(), Stats->TSPTime.end(), 0.0) / Stats->TSPTime.size();
		}
		Stats->TSPSolved = Stats->TSPNodes.size();
		AddResultsToFile();
	} catch (IloException& e) {
	std::cerr << "Concert Exception: " << e << std::endl;
	}
	
	//printf("\n\nType\tNodes\tTime");
	//for (int i = 0; i < Stats->TSPNodes.size(); i++)
	//{
	//	std::cout << "\n" << Stats->InequalityType[i] << "\t" << Stats->TSPNodes[i] << "\t" << Stats->TSPTime[i];
	//	
	//}
}
//void MFGVRP_Solver::GetSolution() {
//	
//	for (IloInt i = 0; i < n+1; i++)
//	{
//		Cplex.getValues(ySol1[i], y[i]);
//		for (IloInt j = 0; j < n+1; j++)
//		{
//			Cplex.getValues(xSol1[i][j], x[i][j]);
//		}
//	}
//}

void MFGVRP_Solver::PrintSolution() {
	int K = 1;
	double ICEV_DIST = 0;
	if (C > 0) K = 2;
	std::cout << "\n\n";
	for (IloInt k = 0; k < K; k++)
	{
		std::cout << "\n\nVehicle " << k << " variables:";
		for (IloInt i = 0; i < n + 1; i++)
		{
			std::cout << "\n";
			for (IloInt j = 0; j < n + 1; j++)
			{
				std::cout << Cplex.getValue(x[i][j][k]) << " ";
				if (k > 0) ICEV_DIST += c[i][j] * Cplex.getValue(x[i][j][k]);
			}
		}
	}
	std::cout << "\n\nY variables:";
	int chargerused;
	for (IloInt i = 0; i < n + 1; i++)
	{
		std::cout << "\n";
		for (IloInt j = 0; j < n + 1; j++)
		{
			chargerused = 0;
			for (IloInt r = 0; r < R[i][j].size(); r++)
			{
				if (Cplex.getValue(y[i][j][r]) > 0.9) chargerused = R[i][j][r];
			}
			std::cout << chargerused << " ";
		}
	}

	///*std::cout << "\n\nTime variables:";

	//std::cout << "\n";
	//for (IloInt j = 1; j < n + 1; j++)
	//{
	//	std::cout << Cplex.getValue(time[j]) << " ";
	//}*/

	printf("\n\nICEV DISTANCE = %f", ICEV_DIST);
}


int MFGVRP_Solver::NonDominatedChargingPaths() {
	//std::vector<std::vector<std::vector<int>>> R(n+1, std::vector<std::vector<int>>(n+1, std::vector<int>(1,NULL))); // Non-dominated charging paths.
	bool NonDominated;
	double MinDist;
	int MinChar;
	int numY = 0;
	int Depot = 0;
	for (IloInt i = 0; i < n + 1; i++)
	{
		Rj.clear();
		for (IloInt j = 0; j < n + 1; j++)
		{
			Rij.clear();
			if (i == 0 || j == 0) Depot = 1; else Depot = 0;
			for (IloInt r = n + 1 + Depot; r < n + f + 1; r++)
			{
				NonDominated = true;

				for (IloInt rr = n + 1+Depot; rr < n + f + 1; rr++)
				{

					//if ((i == 0 && r == n + 1) || (j == 0 && r == n + 1)) { NonDominated = false; continue; }

					if ((c[r][i] >= c[rr][i] && c[r][j] >= c[rr][j] && r != rr)) { NonDominated = false; break; }
				}
				if (NonDominated) { Rij.push_back(r); numY++; }
			}

			Rj.push_back(Rij);
		}
		MinDist = 9999;
		for (IloInt r = n + 1; r < n + f + 1; r++)
		{
			if (MinDist > c[r][i]) {
				MinDist = c[r][i];
				MinChar = r;
			}
		}
		if (i == 0)MinCharDist.push_back(0);
		else MinCharDist.push_back(MinDist);
		R.push_back(Rj);
	}
	//for (int i = 0; i < f; i++)
	//{
	//	MinCharDist.push_back(0);
	//}
	//std::cout << "\n\n";
	//for (IloInt i = 0; i < n + 1; i++)
	//{
	//	std::cout << "\n";
	//	for (IloInt j = 0; j < n + 1; j++)
	//	{

	//		std::cout << "{-";
	//		for (IloInt r = 0; r < R[i][j].size(); r++)
	//		{
	//			std::cout << R[i][j][r] << "-";
	//		}
	//		std::cout << "} ";
	//	}
	//}
	////std::sort(MinCharDist.begin(), MinCharDist.end());
	//std::cout << "\n\nIndex\tDistance:";
	//for (int i = 0; i < n + 1; i++)
	//{
	//	std::cout << "\n" << i << "\t" << MinCharDist[i];
	//}
	return numY;
}



//void MFGVRP_Solver::TournamentInequalities(int NoOfCustomers, int Battery, int NoOfEdges, int* EdgeTail, int* EdgeHead, double* EdgeX, int* ChargingStation, CnstrMgrPointer MyCuts) {
//	std::vector<bool> temp(NoOfEdges, false);
//	std::vector<int> P;
//	std::vector<std::vector<bool>> res;
//	std::vector<std::vector<int>> resP;
//	//std::vector<int> CutTail;
//	//std::vector<int> CutHead;
//	double EnergyConsumption;
//	double MinCharConsumption; double SecondMinCharConsumption; //Minimum and second minimum distance to charging stations from P
//	int CloseChargersFound = 0;
//	std::cout << "\n\nOld algo path:\n";
//	FindCutSets(NoOfEdges, EdgeTail, EdgeHead, EdgeX, ChargingStation, 0, 0, 0, res, temp, resP, P, -1, -1, false, false, MyCuts);
//
//	//std::cout << "\n\nEdge sets: ";
//	//for (int i = 0; i < res.size(); i++)
//	//{
//	//	std::cout << "\n";
//	//	for (int j = 0; j < res[i].size(); j++)
//	//	{
//	//		std::cout << res[i][j] << " ";
//	//	}
//	//}
//
//	//std::cout << "\n\nP sets: ";
//	//for (int i = 0; i < resP.size(); i++)
//	//{
//	//	std::cout << "\n";
//	//	for (int j = 0; j < resP[i].size(); j++)
//	//	{
//	//		std::cout << resP[i][j] << " ";
//	//	}
//	//}
//
//	////Check if path is infeasible:
//	//for (int i = 0; i < res.size(); i++)
//	//{
//	//	CloseChargersFound = 0;
//	//	CutTail.clear();
//	//	CutHead.clear();
//	//	for (int j = 0; j < n+1; j++)
//	//	{
//	//		if (std::count(resP[i].begin(), resP[i].end(), MinCharDist[j].second))
//	//		{
//	//			CloseChargersFound++;
//	//			if (CloseChargersFound == 1) MinCharConsumption = g*MinCharDist[j].first;
//	//			else if (CloseChargersFound == 2) SecondMinCharConsumption = g*MinCharDist[j].first;
//	//			else break;
//	//		}
//	//	}
//	//	EnergyConsumption = 0;
//	//	
//	//	for (int j = 0; j < res[i].size(); j++) {
//	//		if (res[i][j]) {
//	//			EnergyConsumption += g*c[EdgeTail[j]][EdgeHead[j]];
//	//		}
//	//	}
//	//	if ((EnergyConsumption + MinCharConsumption + SecondMinCharConsumption) > Battery) {
//	//		//std::cout << "\n\nHey we have a violation! Edge set " << i << " are using " << EnergyConsumption << " min char=" << MinCharConsumption << " and second min char=" << SecondMinCharConsumption;
//	//		//EnergyConsumption=MinCharConsumption+SecondMinCharConsumption+RunTSPProblem(resP[i].data(), resP[i].size());
//	//		//std::cout << "\n\nMinimum consumption based on TSP problem is: " << EnergyConsumption<<" Battery capacity is: "<<Battery;
//	//		//if (EnergyConsumption>Battery)
//	//		//{
//	//			//resP[i].insert(resP[i].begin(), 0); //Just insert an extra element, so that we have the right format
//	//			//CMGR_AddCnstr(MyCuts, CMGR_ENERGY_SET, 0, resP[i].size()-1, resP[i].data(), resP.size() - 2);
//	//		//}
//	//		//else
//	//		//{
//	//			CutTail.push_back(0); CutHead.push_back(0); //Inserting an extra element to align with format
//	//			for (int j = 0; j < res[i].size(); j++)
//	//				if (res[i][j]) {
//	//					CutTail.push_back(EdgeTail[j]);
//	//					CutHead.push_back(EdgeHead[j]);
//	//				}
//	//		//	//for (IloInt j = 1; j < CutTail.size(); j++)
//	//		//	//{
//	//		//	//	std::cout << CutTail[j] << " ";
//	//		//	//}
//	//			CMGR_AddPathCnstr(MyCuts, CMGR_ENERGY_PATH, 0, CutTail.size()-1, CutTail.data(), CutHead.data(), CutTail.size() - 2);
//	//		//	
//	//		//}
//	//	}
//	//}
//
//	//CutTail.end();
//	//CutHead.end();
//
//};

//Set LatestInclusion to -1 when calling the procedure
//void MFGVRP_Solver::FindCutSets(int NoOfEdges, int* EdgeTail, int* EdgeHead, double* EdgeX, int* ChargingStation, double slack, double esum, double tsum, std::vector<std::vector<bool>>& res, std::vector<bool> temp, std::vector<std::vector<int>>& resP, std::vector<int> P, int index, int LatestInclusion, bool AddEdge, bool AddToP, CnstrMgrPointer MyCuts) {

//	if (AddEdge) //Add edge if new edge can be added to set
//	{
//		temp[index] = true;
//		if (AddToP)
//		{
//			if (P.empty()) { P.push_back(EdgeTail[LatestInclusion]); P.push_back(EdgeHead[LatestInclusion]); }
//			else if (std::count(P.begin(), P.end(), EdgeTail[LatestInclusion])) {
//				if (EdgeTail[LatestInclusion] != P[0]) P.push_back(EdgeHead[LatestInclusion]);
//				else P.insert(P.begin(), EdgeHead[LatestInclusion]);
//			}
//			else {
//				if (EdgeHead[LatestInclusion] != P[0]) P.push_back(EdgeTail[LatestInclusion]);
//				else P.insert(P.begin(), EdgeTail[LatestInclusion]);
//
//			}
//		}
//	}
//
//	index = (index + 1) % NoOfEdges; //Update index
//	if (P.empty()) //If we have not inlcuded any edges, then start by creating a new set based on each edge
//	{
//		for (IloInt i = 0; i < NoOfEdges; i++)
//		{
//			if (EdgeTail[i] == 0) continue;
//			LatestInclusion = index = i;
//			slack = (1 - EdgeX[LatestInclusion]);
//			tsum = s[EdgeTail[i]] + s[EdgeHead[i]] + c[EdgeTail[i]][EdgeHead[i]];
//			esum = g * c[EdgeTail[i]][EdgeHead[i]];
//			FindCutSets(NoOfEdges, EdgeTail, EdgeHead, EdgeX, ChargingStation, slack, esum, tsum, res, temp, resP, P, index, i, true, true, MyCuts);
//		}
//	}
//	else if (tsum + c[0][P[0]] + c[P[P.size() - 1]][0] > T || esum + g * MinCharDist[P[0]] + g * MinCharDist[P[P.size() - 1]] > B) //If we have found a violated set (on time or energy), we include a cut
//	{
//		//for (int i = 0; i < res.size(); i++) if (res[i] == temp) return; //If edgeset already included then don't add set
//		int t = 0;
//		res.push_back(temp);
//		P.insert(P.begin(), 0);
//
//		//std::cout << "\n\nOld algo path:\n";
//		std::cout << "\n";
//		for (int i = 0; i < P.size(); i++)
//		{
//			std::cout << P[i] << " ";
//		}
//
//
//		IloIntArray rhs(env, n * 2 + 2);
//
//		for (int i = 0; i < P.size(); i++) {
//			rhs[P[i]] = 1;
//			rhs[P[i] + n + 1] = 1;
//		}
//
//		//if (RunTSPProblem(rhs, P.data(), P.size()) > B) CMGR_AddCnstr(MyCuts, CMGR_ENERGY_SET, 0, P.size() - 1, P.data(), P.size() - 3);
//		//else CMGR_AddCnstr(MyCuts, CMGR_ENERGY_PATH, 0, P.size() - 1, P.data(), P.size() - 3);
//
//		//Initiate RHS for asymmetric TSP
//		rhs.end();
//
//		resP.push_back(P);
//		return;
//	}
//	else if (LatestInclusion == index) return; //If we were not able to extend the path with any of non-added , then close the path.
//	else if (EdgeTail[index] == 0) //We leave out depot edges
//	{
//		FindCutSets(NoOfEdges, EdgeTail, EdgeHead, EdgeX, ChargingStation, slack, esum, tsum, res, temp, resP, P, index, LatestInclusion, false, false, MyCuts);
//	}
//	else if (((std::count(P.begin(), P.end(), EdgeTail[index]) && !std::count(P.begin(), P.end(), EdgeHead[index])) //If EdgeTail is  in p and EdgeHead is not
//		|| (!std::count(P.begin(), P.end(), EdgeTail[index]) && std::count(P.begin(), P.end(), EdgeHead[index]))) //Or EdgeTail is not in p and EdgeHead is
//		&& (slack + 1 - EdgeX[index]) < 0.999) //And we do not exceed the slack limit
//	{
//		LatestInclusion = index;
//		slack += (1 - EdgeX[LatestInclusion]);
//		esum += g * c[EdgeTail[LatestInclusion]][EdgeHead[LatestInclusion]];
//		if (P[0] == EdgeTail[LatestInclusion] || P[P.size() - 1] == EdgeTail[LatestInclusion]) tsum += s[EdgeHead[LatestInclusion]] + c[EdgeTail[LatestInclusion]][EdgeHead[LatestInclusion]];
//		else tsum += s[EdgeTail[LatestInclusion]] + c[EdgeTail[LatestInclusion]][EdgeHead[LatestInclusion]];
//		FindCutSets(NoOfEdges, EdgeTail, EdgeHead, EdgeX, ChargingStation, slack, esum, tsum, res, temp, resP, P, index, LatestInclusion, true, true, MyCuts);
//	}
//	else //Else we do not add the current edge and go to the next edge 
//	{
//		FindCutSets(NoOfEdges, EdgeTail, EdgeHead, EdgeX, ChargingStation, slack, esum, tsum, res, temp, resP, P, index, LatestInclusion, false, false, MyCuts);
//	}
//}


double MFGVRP_Solver::RunTSPProblem(IloIntArray rhs, int* P, int NoOfCustomers, const std::string &CutType, int* FixedEdge, bool AlterObjective, bool OneWay, bool time) {

	double ObjVal = 0;
	
	//Initializing TSP model:
	TSPModel = IloModel(env);
	TSPCplex = IloCplex(TSPModel);
	double start = TSPCplex.getTime();
	if (AlterObjective)
	{
		//TSPModel.remove(TSPObj);
		//Alter objective
		//if (time) {
		//	for (IloInt i = 0; i < n + 1; i++)
		//		for (IloInt j = 0; j < n + 1; j++)
		//			TSPObj.setLinearCoef(xTSP[i][j], c[i][j]);
		//}
		if (OneWay) for (IloInt i = 1; i < NoOfCustomers; i++) TSPObj.setLinearCoef(xTSP[0][P[i]], g * c[P[0]][P[i]]);
		else
		{
			for (IloInt i = 1; i < NoOfCustomers; i++) {
				TSPObj.setLinearCoef(xTSP[0][P[i]], g * c[P[0]][P[i]]);
				TSPObj.setLinearCoef(xTSP[P[i]][0], g * c[P[i]][P[NoOfCustomers]]);
			}
		}
	}
	if (time) { TSPModel.add(TSPObjTime); }
	else { TSPModel.add(TSPObj); TSPCutOff = B; }
	//for (IloExpr::LinearIterator it = IloExpr(TSPObj.getExpr()).getLinearIterator(); it.ok(); ++it)
	//	std::cout << "variable: " << it.getVar() << " coef: " << it.getCoef() << std:: endl;
	//Update RHS:
	TSPCons.setBounds(rhs, rhs);

	if (FixedEdge != nullptr) { //Fixed edges
		if (OneWay) xTSP[0][FixedEdge[0]].setLB(1);
		else {
			xTSP[0][FixedEdge[0]].setLB(1);
			xTSP[FixedEdge[1]][0].setLB(1);
		}
	}

	//Save relevant customers
	setTSPCustomers(P, NoOfCustomers);
	TSPCustomers.insert(TSPCustomers.begin(), 0);
	TSPModel.add(TSPCons);


	//TSPCplex.extract(TSPModel);
	//TSPCplex.exportModel("MY_TSP.lp");

	TSPCplex.use(LazyTSPCuts(env, *this));
	TSPCplex.setOut(env.getNullStream());
	TSPCplex.setWarning(env.getNullStream());

	TSPCplex.solve();
	ObjVal = TSPCplex.getObjValue();
	Stats->TSPNodes.push_back(TSPCustomers.size());
	Stats->InequalityType.push_back(CutType);
	//printf("\n%s", CutType);
	//printf("\nTSP: %f\n0 ",ObjVal);
	/*printf("\nTSP SOLUTION WITH OBJ VAL %f:\n", ObjVal);
	int CurNode = 0;
	for (int i = 0; i <= n+f; i++)
	{
		if (TSPCplex.getValue(xTSP[CurNode][i]) > 0.9) {
			CurNode = i;
			printf("%d ", CurNode);
			while (CurNode != 0) {
				for (int j = 0; j <= n + f; j++)
				{
					if (TSPCplex.getValue(xTSP[CurNode][j]) > 0.9) {
						CurNode = j;
						printf("%d ", CurNode);
						break;
					}
				}
			}
			break;
		}
	}*/

	if (FixedEdge != nullptr) {
		if (OneWay) xTSP[0][FixedEdge[0]].setLB(0);
		else {
			xTSP[0][FixedEdge[0]].setLB(0);
			xTSP[FixedEdge[1]][0].setLB(0);
		}
	}


	//if (FixedEdge != nullptr) printf("\nThis was a Fixed edge TSP");
	//else if (AlterObjective) printf("\nThis was a Fixed Charger TSP");

	//printf("\nCustomers were: ");
	//for (int i = 0; i < TSPCustomers.size(); i++) printf("%d ", TSPCustomers[i]);
	if (AlterObjective)
	{
		////Reset objective:
		TSPModel.remove(TSPObj);
		//if (time) {
		//	for (IloInt i = 0; i < n + 1; i++)
		//		for (IloInt j = 0; j < n + 1; j++)
		//			if (i == 0 && j == 0) TSPObj.setLinearCoef(xTSP[i][j], 999);
		//			else if (i == 0)TSPObj.setLinearCoef(xTSP[i][j], g * MinCharDist[j]);
		//			else if (j == 0) TSPObj.setLinearCoef(xTSP[i][j], g * MinCharDist[i]);
		//			else TSPObj.setLinearCoef(xTSP[i][j], g * c[i][j]);
		//}
		if (OneWay) for (IloInt i = 1; i < NoOfCustomers; i++) TSPObj.setLinearCoef(xTSP[0][P[i]], g * MinCharDist[P[i]]);
		else
		{
			for (IloInt i = 1; i < NoOfCustomers; i++) {
				TSPObj.setLinearCoef(xTSP[0][P[i]], g * MinCharDist[P[i]]);
				TSPObj.setLinearCoef(xTSP[P[i]][0], g * MinCharDist[P[i]]);
			}
		}
		TSPModel.add(TSPObj);
	}

	Stats->TSPTime.push_back(TSPCplex.getTime() - start);
	TSPModel.end();
	TSPCplex.end();
	return ObjVal;
}

void MFGVRP_Solver::SetupTSP() {



	//TSPCplex.setParam(IloCplex::SolnPoolCapacity, 1);
	//TSPCplex.setParam(IloCplex::MemoryEmphasis, true);
	//TSPCplex.setParam(IloCplex::WorkMem, 500.0);
	//TSPCplex.setParam(IloCplex::SolnPoolReplace, 1);
	TSPObj = IloObjective(env);
	TSPObjTime = IloObjective(env);

	int count = 0;
	TSPCons = IloRangeArray(env);
	IloExpr expr(env);
	IloExpr expr1(env);
	//Variables:
	xTSP = IloVar2DMatrix(env, n + f + 1);

	for (IloInt i = 0; i < n + f + 1; i++) {
		xTSP[i] = IloNumVarArray(env, n + f + 1, 0, 1, ILOINT);
		for (IloInt j = 0; j < n + 1; j++)
		{
			std::string xName = "x_" + std::to_string(i) + "_" + std::to_string(j);
			xTSP[i][j].setName(xName.c_str());
		}
	}

	//Set objective:
	for (IloInt i = 0; i < n + f + 1; i++)
		for (IloInt j = 0; j < n + f + 1; j++)
			if (i == j) { expr += 999*xTSP[i][j]; expr1 += 999 * xTSP[i][j]; xTSP[i][j].setUB(0); }
			else if (i == 0) { if(j) expr += g * MinCharDist[j] * xTSP[i][j]; expr1 += c[i][j] * xTSP[i][j];}
			else if (j == 0) { expr += g * MinCharDist[i] * xTSP[i][j]; expr1 += c[i][j] * xTSP[i][j];}
			else { expr += g * c[i][j] * xTSP[i][j]; expr1 += c[i][j] * xTSP[i][j]; }

	TSPObj = IloMinimize(env,expr);
	TSPObjTime = IloMinimize(env, expr1);
	//Add out cons
	for (IloInt i = 0; i < n + f + 1; i++)
	{
		expr.clear();
		for (IloInt j = 0; j < n + f + 1; j++)
		{
			expr += xTSP[i][j];
		}
		TSPCons.add(expr == 1);
		++count;
	}
	//Add in cons
	for (IloInt i = 0; i < n + f + 1; i++)
	{
		expr.clear();
		for (IloInt j = 0; j < n + f + 1; j++)
		{
			expr += xTSP[j][i];
		}
		TSPCons.add(expr == 1);
		++count;
	}

	//TSPModel.add(TSPCons);
	//TSPCplex.use(LazyTSPCuts(env, *this));
	//TSPCplex.setOut(env.getNullStream());
	//TSPCplex.setWarning(env.getNullStream());
}

//Assumption: Edgetail[0], EdgeHead[0] and EdgeX[0] are all set to 0.
void MFGVRP_Solver::FindViolatedPaths(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts,bool isRoot) {
	int t, i, j, root, L, idx;
	
	int CutIdx = MyCuts->Size;
	bool PathAlreadyAdded;
	double ExtraSlack = 0;
	double ExtraSetSlack = 0;
	std::vector<int> LevelNodeVector;
	for (i = 0; i < n + f + 1; i++) AddedToPath[i] = false;



	//std::cout << "\n\nNew algo paths:\n";
	for (root = 1; root < n + 1; root++) {

		LevelNodeVector.clear(); //Clear vector for new node
		LevelNodeVector.push_back(0);
		LevelNodeVector.push_back(root);
		LevelNodeVector.push_back(0);
		double GreedyObj = 0;
		AddedToPath[root - 1] = false; //Remove previous customer to 
		AddedToPath[root] = true; //Set added to path = 1;
		LevelNode[1] = root;  //First customer on path set 
		LevelNode[0] = 0; //First node is Depot
		Esum[1] = g * MinCharDist[root]; //Set to min dist to a charging station initially (Should be changed if we do not consider energy consumption direct proportional with distance)
		Tsum[1] = c[0][root] + s[root]; //Set initial time consumption to time consumption from depot (Should be changed if we do not consider time direct proportional with distance)
		Slack[1] = 0;
		SlackSet[1] = 0;
		L = 2; //Set level to search = to 2
		LevelIdx[L] = 1; //Start the search in the begining 
		LevelNode[L] = 0; //Level 2 node to 0 until we find another node to add to path
		do
		{
			while (LevelIdx[L] <= NoOfEdges) {
				/*printf("\nEdge ID %d, Edge tail %d, Edge head %d", LevelIdx[L], EdgeTail[LevelIdx[L]], EdgeHead[LevelIdx[L]]);*/
				//Check if a non explored customer is connected to current customer
				idx = root;
				/*				if (LevelNode[L - 1] == EdgeHead[LevelIdx[L]]) {

									idx = EdgeTail[LevelIdx[L]];
								}
								else*/ if (LevelNode[L - 1] == EdgeTail[LevelIdx[L]]) {
									idx = EdgeHead[LevelIdx[L]];
								}

								if (!AddedToPath[idx])
								{
									LevelNode[L] = idx;
									LevelNodeVector[L] = idx;
									//SlackSet[L] = SlackSet[L - 1];
									ExtraSetSlack = 1;
									ExtraSlack = 1 - EdgeX[LevelIdx[L]];

									if (LevelNode[L] == EdgeTail[LevelIdx[L] - 1] && LevelNode[L - 1] == EdgeHead[LevelIdx[L] - 1])
										ExtraSlack -= EdgeX[LevelIdx[L] - 1];
									else if (LevelIdx[L] + 1<=NoOfEdges) 
										if(LevelNode[L] == EdgeTail[LevelIdx[L] + 1] && LevelNode[L - 1] == EdgeHead[LevelIdx[L] + 1])
											ExtraSlack -= EdgeX[LevelIdx[L] + 1];
									if (isRoot)
									{
										//Add additional contribution to LHS for set expansion
										for (i = 1; i <= NoOfEdges; i++) {
											//if (i == LevelIdx[L]) continue;
											if (EdgeTail[i] == idx) {
												if (std::find(LevelNodeVector.begin() + 1, LevelNodeVector.end(), EdgeHead[i]) != LevelNodeVector.end())  ExtraSetSlack -= EdgeX[i];
											}
											if (EdgeHead[i] == idx) {
												if (std::find(LevelNodeVector.begin() + 1, LevelNodeVector.end(), EdgeTail[i]) != LevelNodeVector.end()) ExtraSetSlack -= EdgeX[i];
											}
										}
									}
									else
									{
										ExtraSetSlack = ExtraSlack;
									}
									

									//If slack is strictly less than 1:
									if (Slack[L - 1] + ExtraSlack < 0.99 || SlackSet[L - 1] + ExtraSetSlack < 0.99) {
										/*Violation still possible*/


										//We update the energy and time sum (Should be updated, if we consider time and energy not direct proportional with distance)
										Esum[L] = Esum[L - 1] + g * c[LevelNode[L - 1]][LevelNode[L]];
										Tsum[L] = Tsum[L - 1] + c[LevelNode[L - 1]][LevelNode[L]] + s[LevelNode[L]];
										Slack[L] = Slack[L - 1] + ExtraSlack;
										SlackSet[L] = SlackSet[L - 1] + ExtraSetSlack;

										// If path is violated
										if (Esum[L] + g * MinCharDist[LevelNode[L]] > B)// || Tsum[L] + c[LevelNode[L]][0] > T)
										{ /*Violation found*/

											//////std::cout << "\n\nNew algo path:\n";
											//std::cout << "\n";
											//for (j = 0; j <= L; j++)
											//{
											//	std::cout << LevelNode[j] << " ";
											//}

											//Check if we have already added path
											PathAlreadyAdded = false;

											if (CutIdx < MyCuts->Size)
												for (j = CutIdx; j <= MyCuts->Size - 1; j++)
													if (MyCuts->CPL[j]->IntListSize == L)
														for (t = 1; t <= L; t++)
															if (MyCuts->CPL[j]->IntList[t] != LevelNode[L - t + 1]) break;
															else if (t == L) PathAlreadyAdded = true;


											//CMGR_AddCnstr(MyCuts, CMGR_ENERGY_PATH, 0, L, LevelNode, L - 2);

											//If Path is not already added
											if (!PathAlreadyAdded) {

												////std::cout << "\n\nNew algo path:\n";
												//std::cout << "\n";
												//for (j = 0; j <= L; j++)
												//{
												//	std::cout << LevelNode[j] << " ";
												//}

												//CMGR_AddCnstr(MyCuts, CMGR_ENERGY_PATH, 0, L, LevelNode, L - 2);
												//printf("\nSlack for set %f --- Slack for path %f\n Set:\t", SlackSet[L], Slack[L]);
												//for (i = 0; i < LevelNodeVector.size()-1; i++) printf("%d ", LevelNodeVector[i]);
												//printf("%d ", idx);
												//Set RHS for ATSP

												if (isRoot)
												{
													if (L<=300)
													{
														GreedyObj = Greedy(LevelNode, L, true, false);
														if (GreedyObj > B)
														{
															//printf("\n");
															//for (int i = 1; i <= L; i++)
															//{
															//	printf("%d ", LevelNode[i]);
															//}
															for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; }

															if (RunTSPProblem(rhs, LevelNode, L + 1, "No_charge") > B) {
																CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_SET, 0, L, LevelNode, L - 2);
																for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
																CutCounter++;
																(LevelIdx[L])++;
																continue;
																//return;

															}
															for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
														}
													}
													
													if (Slack[L] < 0.999) {
														CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_PATH, 0, L, LevelNode, L - 2);
														CutCounter++;
														//return;
													}
													//printf("\n We got through it");
												}
												else
												{
													if (L<=300 & isRoot)
													{
														GreedyObj = Greedy(LevelNode, L, true, false);
														if (GreedyObj > B)
														{
															for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; }

															if (RunTSPProblem(rhs, LevelNode, L + 1, "No_charge") > B) {
																CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_SET, 0, L, LevelNode, L - 2);
																for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
																CutCounter++;
																return;
															}
															for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
														}
													}
													
													if (Slack[L] < 0.999) {
														CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_PATH, 0, L, LevelNode, L - 2);
														CutCounter++;
														return;
													}
												}

												//Add set or path violation based on ATSP solution 

												
												//if (!isRoot && CutCounter>=10)
												//{
												//	return;
												//}
												//If ATSP violated, then add set violation
												//Else add path violation


												//Reset RHS for ATSP

											}
										}
										else
										{ /*Path is still feasible, but could be violated when extended*/
											LevelNodeVector[L] = idx;
											AddedToPath[LevelNode[L]] = true; //Add customer to path
											L++; //Go to next level

											LevelIdx[L] = 0; //Set index at next level to 0
											LevelNode[L] = 0; //Set node at next level to 0
											LevelNodeVector.push_back(0);
										}
									}
								}
								(LevelIdx[L])++;
			}

			L--;
			LevelNodeVector.pop_back();
			AddedToPath[LevelNode[L]] = false;
			LevelNode[L] = 0;
			(LevelIdx[L])++;
		} while (L > 1);

	}

	//rhs.end();
}

//Assumption: Edgetail[0], EdgeHead[0] and EdgeX[0] are all set to 0.
//void MFGVRP_Solver::FindViolatedDepotPaths(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<int>& Charger, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts) {
//	int t, i, j, L, idx, P;
//	IloIntArray rhs(env, n * 2 + 2);
//	double ReturntoDepot;
//	bool PathAlreadyAdded, OnArc;
//	bool* EdgeAdded = new bool[NoOfEdges + 1];
//	int* FixedEdge = new int[2];
//	std::vector<int> TSPCustomers;
//
//
//	for (int i = 0; i <= NoOfEdges; i++) EdgeAdded[i] = false;
//	//std::cout << "\n\nNew algo paths:\n";
//
//	LevelNode[1] = 0;  //First customer on path set 
//	LevelNode[0] = 0; //First node is Depot
//	Esum[1] = 0; //Set to min dist to a charging station initially (Should be changed if we do not consider energy consumption direct proportional with distance)
//	Tsum[1] = 0; //Set initial time consumption to time consumption from depot (Should be changed if we do not consider time direct proportional with distance)
//	Slack[1] = 0; //Set intial slack to 0
//	L = 2; //Set level to search to 2
//	LevelIdx[L] = 1; //Start the search in the begining 
//	LevelNode[L] = 0; //Level 2 node to 0 until we find another node to add to path
//
//	do
//	{
//		while (LevelIdx[L] <= NoOfEdges) {
//			/*printf("\nEdge ID %d, Edge tail %d, Edge head %d", LevelIdx[L], EdgeTail[LevelIdx[L]], EdgeHead[LevelIdx[L]]);*/
//			//Check if a non explored customer is connected to current customer
//			OnArc = false;
//
//			if (LevelNode[L - 1] == EdgeHead[LevelIdx[L]] && !EdgeAdded[LevelIdx[L]] && Charger[LevelIdx[L]] == -1) {
//				idx = EdgeTail[LevelIdx[L]];
//				OnArc = true;
//			}
//			else if (LevelNode[L - 1] == EdgeTail[LevelIdx[L]] && !EdgeAdded[LevelIdx[L]] && Charger[LevelIdx[L]] == -1) {
//				idx = EdgeHead[LevelIdx[L]];
//				//std::cout << "\n" << EdgeHead[LevelIdx[L]];
//				OnArc = true;
//			}
//
//			if (!EdgeAdded[LevelIdx[L]] && OnArc)
//			{
//				LevelNode[L] = idx;
//				EdgeAdded[LevelIdx[L]] = true;
//
//
//				//If slack is strictly less than 1:
//				if (Slack[L - 1] + (1 - EdgeX[LevelIdx[L]]) < 0.999) {
//					/*Violation still possible*/
//
//
//
//					//We update the energy and time sum (Should be updated, if we consider time and energy not direct proportional with distance)
//					Esum[L] = Esum[L - 1] + g * c[LevelNode[L - 1]][LevelNode[L]];
//
//					if (idx != 0) Tsum[L] = Tsum[L - 1] + c[LevelNode[L - 1]][LevelNode[L]] + s[LevelNode[L]];
//					else Tsum[L] = Tsum[L - 1] + c[LevelNode[L - 1]][LevelNode[L]];
//
//					if (LevelNode[L] == 0) ReturntoDepot = 0; else ReturntoDepot = c[LevelNode[L]][0];
//					Slack[L] = Slack[L - 1] + 1 - EdgeX[LevelIdx[L]];
//					// If path is violated
//					if (Esum[L] + g * MinCharDist[LevelNode[L]] > B || Tsum[L] + ReturntoDepot > T)
//					{ /*Violation found*/
//
//						//If we are back at the depot 
//						if (idx == 0) {
//							TSPCustomers.clear();
//							//Set RHS and customers for for ATSP
//							for (j = 2; j <= L - 1; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + 1] = 1; TSPCustomers.push_back(LevelNode[j]); }
//
//							//Set fixed Edge for ATSP
//							FixedEdge[0] = LevelNode[2]; FixedEdge[1] = LevelNode[L - 1];
//
//							//Add set or path violation based on ATSP solution 
//							if (RunTSPProblem(rhs, TSPCustomers.data(), TSPCustomers.size(), "Depot_Path", FixedEdge) + g * c[0][FixedEdge[0]] + g * c[FixedEdge[1]][0] - g * c[FixedEdge[0]][FixedEdge[1]] > B) {
//								TSPCustomers.insert(TSPCustomers.begin(), 0);
//								CMGR_AddCnstr(MyCuts, CMGR_ENERGY_DEPOT, 0, TSPCustomers.size() - 1, TSPCustomers.data(), 0);
//							}
//							else CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_PATH, 0, L, LevelNode, L - 2);
//
//							//Reset RHS for ATSP
//							for (j = 2; j <= L - 1; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + 1] = 0; }
//							break;
//						}
//						else CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_PATH, 0, L, LevelNode, L - 2);
//					}
//					//If we are back at the depot 
//					else if (idx == 0) break;
//					else
//					{ /*Path is still feasible, but could be violated when extended*/
//						AddedToPath[LevelNode[L]] = true; //Add customer to path
//
//						L++; //Go to next level
//
//						LevelIdx[L] = 0; //Set index at next level to 0
//						LevelNode[L] = 0; //Set node at next level to 0
//					}
//				}
//			}
//			(LevelIdx[L])++;
//		}
//
//		L--;
//		(LevelIdx[L])++;
//
//	} while (L > 1);
//	//TSPCustomers.end();
//	//rhs.end();
//	delete[] EdgeAdded;
//	delete[] FixedEdge;
//}

//Assumptions: Edgetail[0], EdgeHead[0] and EdgeX[0] are all set to 0.
//			   Inequalities have been added on infeasible charging paths
//void MFGVRP_Solver::CheckRouteFeasibility(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<int>& Charger, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts) {
//	int cnt, i, j, t, L, idx, P;
//	bool* EdgeAdded = new bool[NoOfEdges + 1];
//	double* NoChargePathConsumption = new double[n + 1];
//	int* LevelCharger = new int[n + 1];
//	double tsum, SOC;
//	std::vector<int> Violated_Path;
//	int ViolatedPathIdx;
//	for (i = 0; i <= NoOfEdges; i++) EdgeAdded[i] = false;
//	//std::cout << "\n\nNew algo paths:\n";
//	tsum = 0;
//	LevelNode[0] = 0; LevelCharger[0] = -1;
//	LevelNode[1] = 0;  //First node is depot
//	LevelCharger[1] = -1; //Set charger at depot =-1
//	L = 2; //Set level to search to 2
//	P = 1; //Set path level on route to 0
//	NoChargePathConsumption[P] = 0;
//	LevelIdx[L] = 1; //Start the search at index 1 
//	LevelCharger[L] = -1; //Initlize with Charger =-1
//	LevelNode[L] = 0; //Level 2 node to 0 until we find another node to add to path
//
//	do
//	{
//		while (LevelIdx[L] <= NoOfEdges) {
//			/*printf("\nEdge ID %d, Edge tail %d, Edge head %d", LevelIdx[L], EdgeTail[LevelIdx[L]], EdgeHead[LevelIdx[L]]);*/
//			//Check if a non explored customer is connected to current customer
//			idx = -1;
//
//			if (LevelNode[L - 1] == EdgeHead[LevelIdx[L]] && !EdgeAdded[LevelIdx[L]]) {
//				idx = EdgeTail[LevelIdx[L]];
//			}
//			else if (LevelNode[L - 1] == EdgeTail[LevelIdx[L]] && !EdgeAdded[LevelIdx[L]]) {
//				idx = EdgeHead[LevelIdx[L]];
//			}
//
//			if (!EdgeAdded[LevelIdx[L]] && idx != -1)
//			{
//				LevelNode[L] = idx;
//				if (Charger[LevelIdx[L]] != -1) LevelCharger[L] = Charger[LevelIdx[L]]; else LevelCharger[L] = -1;
//				EdgeAdded[LevelIdx[L]] = true;
//
//				//Update path requirements:
//				if (LevelCharger[L] != -1) {
//					NoChargePathConsumption[P] += g * c[LevelCharger[L]][LevelNode[L - 1]];
//
//					//printf("\nEnd costumer %d charging station %d with consumption %f", LevelNode[L-1],LevelCharger[L], NoChargePathConsumption[P]);
//					P++;
//					NoChargePathConsumption[P] = g * c[LevelCharger[L]][LevelNode[L]];
//					tsum += c[LevelCharger[L]][LevelNode[L]] + c[LevelCharger[L]][LevelNode[L - 1]];
//					if (LevelNode[L] != 0) tsum += s[LevelNode[L]];
//				}
//				else {
//					NoChargePathConsumption[P] += g * c[LevelNode[L - 1]][LevelNode[L]];
//					tsum += c[LevelNode[L - 1]][LevelNode[L]];
//					if (LevelNode[L] != 0) tsum += s[LevelNode[L]];
//				}
//
//				if (idx == 0) {
//					// If we are back at the depot, then check route feasibility:
//
//
//					// First check each no-charge path violates battery capacity:
//					for (i = 1; i <= P; i++)
//					{
//						if (NoChargePathConsumption[i] > B)
//						{
//							ViolatedPathIdx = 1;
//							cnt = 0;
//							for (j = 2; j <= L; j++)
//							{
//								if (LevelCharger[j] != -1 || LevelNode[j] == 0) { cnt++; if (cnt < i) { ViolatedPathIdx = j; } }
//								if (cnt == i) {
//
//									Violated_Path.clear();
//									Violated_Path.push_back(0);
//									if (LevelCharger[ViolatedPathIdx] != -1) { Violated_Path.push_back(LevelNode[ViolatedPathIdx]); Violated_Path.push_back(LevelCharger[ViolatedPathIdx]); }
//									for (t = ViolatedPathIdx; t < j; t++)
//									{
//										Violated_Path.push_back(LevelNode[t]);
//									}
//									if (LevelCharger[j] != -1) { Violated_Path.push_back(LevelCharger[j]); Violated_Path.push_back(LevelNode[j]); }
//									else Violated_Path.push_back(LevelNode[j]);
//
//									CMGR_AddCnstr(MyCuts, CMGR_FINAL_FEASIBILITY, 0, Violated_Path.size() - 1, Violated_Path.data(), Violated_Path.size() - 3);
//									break;
//								}
//							}
//						}
//					}
//
//
//					// Then check time feasibility if we are using charging stations (else these would have been taken care of in the other constraints):
//
//					if (P > 1) {
//
//						SOC = B - NoChargePathConsumption[1];
//						for (i = 2; i <= P; i++)
//						{
//							tsum += r * std::max(NoChargePathConsumption[i] - SOC, 0.0);
//							SOC = std::max(SOC - NoChargePathConsumption[i], 0.0);
//						}
//						if (ObjValTemporary >= 837.6715 && ObjValTemporary <= 837.6717) printf("\n\nT sum is: %f\n\n", tsum);
//						if (tsum > T) {
//							//Violation found
//							Violated_Path.clear();
//							Violated_Path.push_back(0);
//
//							for (i = 1; i <= L; i++)
//							{
//								if (LevelCharger[i] != -1) Violated_Path.push_back(LevelCharger[i]);
//								Violated_Path.push_back(LevelNode[i]);
//							}
//
//							//Add path cut:
//							CMGR_AddCnstr(MyCuts, CMGR_FINAL_FEASIBILITY, 0, Violated_Path.size() - 1, Violated_Path.data(), L - 2);
//						}
//					}
//
//					P = 1; L = 3;
//					tsum = 0;
//					NoChargePathConsumption[P] = 0;
//					LevelCharger[L] = -1;
//					break;
//					//And we go back to L=3
//				}
//				else
//				{ /*Path is still feasible, but could be violated when extended*/
//					AddedToPath[LevelNode[L]] = true; //Add customer to path
//
//					L++; //Go to next level
//
//					LevelIdx[L] = 0; //Set index at next level to 0
//					LevelNode[L] = 0; //Set node at next level to 0
//					LevelCharger[L] = -1;
//				}
//			}
//			(LevelIdx[L])++;
//		}
//
//		L--;
//		(LevelIdx[L])++;
//
//	} while (L > 1);
//
//	//Violated_Path.end();
//	delete[] EdgeAdded;
//	delete[] LevelCharger;
//	delete[] NoChargePathConsumption;
//}
//Seperates inequalities that relies on finding a set, that we cannot serve with energy feasibility if we drive certain paths to customers:
//void MFGVRP_Solver::TwoNodeInequalities(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<int>& Charger, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts) {
//	int t, i, j, root, L, idx, rootidx, endidx;
//	double rhssum, lhssum;
//	IloIntArray rhs(env, n * 2 + 2);
//	int CutIdx = MyCuts->Size;
//	bool PathAlreadyAdded;
//	int* LevelCharger = new int[n + 1];
//	std::vector<int> RootNodes;
//	std::vector<int> RootChargers;
//	std::vector<int> S;
//	int* FixedEdge = new int[2];
//
//	for (i = 0; i < n + f + 1; i++)
//	{
//		AddedToPath[i] = false;
//	}
//
//	for (i = 0; i <= NoOfEdges; i++)
//	{
//		if (Charger[i] > -1)
//		{
//			if (EdgeTail[i] > 0)
//			{
//				if (std::find(RootNodes.begin(), RootNodes.end(), EdgeTail[i]) == RootNodes.end()) {
//					RootNodes.push_back(EdgeTail[i]);
//					RootChargers.push_back(Charger[i]);
//				}
//			}
//
//			if (std::find(RootNodes.begin(), RootNodes.end(), EdgeHead[i]) == RootNodes.end()) {
//				RootNodes.push_back(EdgeHead[i]);
//				RootChargers.push_back(Charger[i]);
//			}
//		}
//		else if (EdgeTail[i] == 0)
//		{
//			//if (std::find(RootNodes.begin(), RootNodes.end(), EdgeHead[i]) == RootNodes.end()) RootNodes.push_back(std::make_pair(EdgeTail[i], 0));
//			if (std::find(RootNodes.begin(), RootNodes.end(), EdgeHead[i]) == RootNodes.end()) {
//				RootNodes.push_back(EdgeHead[i]);
//				RootChargers.push_back(0);
//			}
//		}
//	}
//
//	//printf("\n\n");
//	//for (int i = 0; i < RootNodes.size(); i++)
//	//{
//	//	printf("\n%d\t%d", RootNodes[i], RootChargers[i]);
//	//}
//	//std::cout << "\n\nNew algo paths:\n";
//	for (root = 1; root < RootNodes.size(); root++) {
//
//		for (i = 0; i <= NoOfEdges; i++)
//		{
//			if (EdgeTail[i] == RootNodes[root] || EdgeHead[i] == RootNodes[root])
//			{
//				if (Charger[i] == RootChargers[root] || EdgeTail[i] == RootChargers[root])
//				{
//					rootidx = i;
//					rhssum = EdgeX[i];
//
//					break;
//
//				}
//			}
//		}
//		//if (Charger[root]<0 && EdgeTail[root]>0) continue; //Continue, if root node is nor the depot or a charger
//		AddedToPath[RootNodes[root - 1]] = false; //Remove previous customer to 
//		AddedToPath[RootNodes[root]] = true; //Set added to path = 1;
//		LevelNode[1] = RootNodes[root];  //First customer on path set
//		LevelNode[0] = 0; //First node is opposite edge
//		Esum[1] = g * c[RootChargers[root]][RootNodes[root]]; //Set dist to initial charging option
//		Slack[1] = 0;
//		L = 2; //Set level to search = to 2
//		LevelIdx[L] = 1; //Start the search in the begining 
//		LevelNode[L] = 0; //Level 2 node to 0 until we find another node to add to path
//		do
//		{
//			while (LevelIdx[L] <= NoOfEdges) {
//				/*printf("\nEdge ID %d, Edge tail %d, Edge head %d", LevelIdx[L], EdgeTail[LevelIdx[L]], EdgeHead[LevelIdx[L]]);*/
//				//Do not include root edge:
//				if (LevelIdx[L] != rootidx)
//				{
//
//
//					//Check if a unexplored customer is connected to the current customer
//					idx = RootNodes[root];
//					if (LevelNode[L - 1] == EdgeHead[LevelIdx[L]]) {
//
//						idx = EdgeTail[LevelIdx[L]];
//					}
//					else if (LevelNode[L - 1] == EdgeTail[LevelIdx[L]] && EdgeTail[LevelIdx[L]] != 0) {
//						idx = EdgeHead[LevelIdx[L]];
//					}
//
//					if (!AddedToPath[idx])
//					{
//						LevelNode[L] = idx;
//
//
//
//						//We update the energy and time sum (Should be updated, if we consider time and energy not direct proportional with distance)
//						if (Charger[LevelIdx[L]] == -1) Esum[L] = Esum[L - 1] + g * c[LevelNode[L - 1]][LevelNode[L]];
//						else Esum[L] = Esum[L - 1] + g * c[Charger[LevelIdx[L]]][LevelNode[L - 1]];
//
//						//Tsum[L] = Tsum[L - 1] + c[LevelNode[L - 1]][LevelNode[L]] + s[LevelNode[L]];
//						if (Charger[LevelIdx[L]] == -1 || !EdgeTail[idx] == 0) Slack[L] = Slack[L - 1] + 1 - EdgeX[LevelIdx[L]];
//						// If path is violated
//						if (Charger[LevelIdx[L]] != -1 || EdgeTail[LevelIdx[L]] == 0)
//						{ /*Violation found*/
//							if (Esum[L] > B)
//							{
//								//If Path is not already added
//								S.clear();
//								endidx = LevelIdx[L];
//								rhssum += EdgeX[LevelIdx[L]];
//								lhssum = 0;
//								S.push_back(RootChargers[root]);
//
//								if (Charger[LevelIdx[L]] != -1) {
//									for (i = 1; i < L; i++) S.push_back(LevelNode[i]);
//									S.push_back(Charger[LevelIdx[L]]);
//								}
//								else for (i = 1; i <= L; i++) S.push_back(LevelNode[i]);
//								//printf("\n\n");
//								//for (i = 0; i < S.size(); i++) printf("%d\n", S[i]);
//
//								for (i = 0; i <= NoOfEdges; i++)
//								{
//									if (i == rootidx || i == endidx)
//										continue;
//									if (std::find(S.begin() + 1, S.end() - 1, EdgeTail[i]) != S.end() - 1 && std::find(S.begin() + 1, S.end() - 1, EdgeHead[i]) == S.end() - 1) lhssum += EdgeX[i];
//									if (std::find(S.begin() + 1, S.end() - 1, EdgeHead[i]) != S.end() - 1 && std::find(S.begin() + 1, S.end() - 1, EdgeTail[i]) == S.end() - 1) lhssum += EdgeX[i];
//								}
//								//printf("\n\n%f", rhssum);
//								//printf("\n\n%f", lhssum);
//
//								//Check if we have already added path
//								PathAlreadyAdded = false;
//
//								if (CutIdx < MyCuts->Size) {
//									for (j = CutIdx; j <= MyCuts->Size - 1; j++) {
//										if (MyCuts->CPL[j]->IntListSize == S.size()) {
//											for (t = 1; t <= S.size(); t++) {
//												if (MyCuts->CPL[j]->IntList[t] != S[S.size() - t]) {
//													break;
//												}
//												else if (t == S.size()) {
//													PathAlreadyAdded = true;
//												}
//											}
//										}
//									}
//								}
//
//								if (lhssum < rhssum && !PathAlreadyAdded)
//								{
//									TSPCustomers.clear();
//									//Set RHS and customers for for ATSP
//									for (j = 1; j < S.size() - 1; j++) { rhs[S[j]] = 1; rhs[S[j] + n + 1] = 1; TSPCustomers.push_back(S[j]); }
//
//									//Set fixed Edge for ATSP
//									FixedEdge[0] = S[1]; FixedEdge[1] = S[S.size() - 2];
//									//printf("\n\n%d, %d, %d, %d", S[0], S[1], S[S.size() - 2], S[S.size() - 1]);
//									//Add set or path violation based on ATSP solution 
//									if (RunTSPProblem(rhs, TSPCustomers.data(), TSPCustomers.size(), "Two_Node", FixedEdge) + g * c[S[0]][FixedEdge[0]] + g * c[S[S.size() - 1]][FixedEdge[1]] - g * c[FixedEdge[0]][FixedEdge[1]] > B) {
//										S.insert(S.begin(), 0);
//										//printf("\n\n%d", MyCuts->Size);
//										CMGR_AddCnstr(MyCuts, CMGR_FIXED_EDGES, 0, S.size() - 1, S.data(), 0);
//										//printf("\n\n%d", MyCuts->Size);
//									}
//
//									//Reset RHS for ATSP
//									for (j = 2; j < S.size() - 1; j++) { rhs[S[j]] = 0; rhs[S[j] + n + 1] = 0; }
//								}
//							}
//
//						}
//						else
//						{ /*Path is still feasible, but could be violated when extended*/
//							AddedToPath[LevelNode[L]] = true; //Add customer to path
//							L++; //Go to next level
//
//							LevelIdx[L] = 0; //Set index at next level to 0
//							LevelNode[L] = 0; //Set node at next level to 0
//						}
//
//					}
//				}
//				(LevelIdx[L])++;
//			}
//
//			L--;
//			AddedToPath[LevelNode[L]] = false;
//			(LevelIdx[L])++;
//		} while (L > 1);
//
//	}
//	delete[] FixedEdge;
//}

//Seperates inequalities based on paths connected to charging options:
void MFGVRP_Solver::FixedChargers(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts, bool isRoot) {
	int t, i, j, root, L, idx;
	int CutIdx = MyCuts->Size;
	bool PathAlreadyAdded;
	int prev_root = 0;
	int* FixedEdge = new int[2];
	double violation = 0;
	double violationFC = 0;
	std::vector<int> CutSet;
	std::vector<int> ExtSet;
	double AdditionalEnergyRequirement = 0;
	double ExtraSlack = 0;

	for (i = 0; i < n + f + 1; i++)
	{
		AddedToPath[i] = false;

	}


	//We loop over all edges.
	for (root = 1; root < NoOfEdges; root++) {

		//If an edge contains a charger or the depot, then we start building a path
		if (EdgeTail[root] == 0 || EdgeTail[root] > n)
		{
			if (EdgeHead[root] == 0 || EdgeHead[root] > n) continue; //If the edge connects the depot and a charger then continue (Assumes reprocessing has taken care of all depot-charger edges that violate the battery capacity)
			ExtraSlack = 0;
			AddedToPath[EdgeTail[prev_root]] = false; //Remove previous nodes
			AddedToPath[EdgeHead[prev_root]] = false; //Remove previous nodes 
			AddedToPath[EdgeHead[root]] = true; //Set added to path = 1 for first customer; 
			LevelNode[1] = EdgeHead[root];  //First customer on path set 
			LevelNode[0] = EdgeTail[root]; //First node is the depot/charger
			for (int i = 1; i < NoOfEdges; i++) if (EdgeHead[root] == EdgeTail[i] && EdgeTail[root] == EdgeHead[i]){	ExtraSlack = EdgeX[i];break;}
			Esum[1] = g * c[EdgeTail[root]][EdgeHead[root]]; //Set to energy consumption to required energy for travelling from initial charger/depot to first customer.
			Tsum[1] = c[EdgeTail[root]][EdgeHead[root]] + s[EdgeHead[root]]; //Set initial time consumption to time consumption from depot (Should be changed if we do not consider time direct proportional with distance)
			Slack[1] = 1 - EdgeX[root] - ExtraSlack;
			
			

			L = 2; //Set level to search = to 2
			LevelIdx[L] = 1; //Start the search in the begining 
			LevelNode[L] = 0; //Level 2 node to 0 until we find another node to add to path
			prev_root = root;
			do
			{
				while (LevelIdx[L] < NoOfEdges) {
					/*printf("\nEdge ID %d, Edge tail %d, Edge head %d", LevelIdx[L], EdgeTail[LevelIdx[L]], EdgeHead[LevelIdx[L]]);*/
					//Check if a non explored customer is connected to current customer
					idx = EdgeHead[root];
					//if (CutCounter >= 10) return;
					if (LevelNode[L - 1] == EdgeTail[LevelIdx[L]] && LevelIdx[L] != prev_root) {
						idx = EdgeHead[LevelIdx[L]];
					}

					if (!AddedToPath[idx])
					{
						LevelNode[L] = idx;
						ExtraSlack = 0;
						if (LevelNode[L] == EdgeTail[LevelIdx[L] - 1] && LevelNode[L - 1] == EdgeHead[LevelIdx[L] - 1])
							ExtraSlack = EdgeX[LevelIdx[L] - 1];
						else if (LevelIdx[L] + 1<NoOfEdges) 
							if(LevelNode[L] == EdgeTail[LevelIdx[L] + 1] && LevelNode[L - 1] == EdgeHead[LevelIdx[L] + 1])
								ExtraSlack = EdgeX[LevelIdx[L] + 1];

						//If slack is strictly less than 1:
						//if (Slack[L - 1] + (1 - EdgeX[LevelIdx[L]]) < 0.999) {
							/*Violation still possible*/


							//We update the energy and time sum (Should be updated, if we consider time and energy not direct proportional with distance)
						Esum[L] = Esum[L - 1] + g * c[LevelNode[L - 1]][LevelNode[L]];
						if (LevelNode[L] <= n) Tsum[L] = Tsum[L - 1] + c[LevelNode[L - 1]][LevelNode[L]] + s[LevelNode[L]]; else Tsum[L] = Tsum[L - 1] + c[LevelNode[L - 1]][LevelNode[L]];
						Slack[L] = Slack[L - 1] + 1 - EdgeX[LevelIdx[L]]-ExtraSlack;
						if (LevelNode[L] <= n) AdditionalEnergyRequirement = MinCharDist[LevelNode[L]]; else AdditionalEnergyRequirement = 0;

						if (Esum[L] + AdditionalEnergyRequirement > B)
						{ /*Potential violation found*/

						  //Initiate vector:
							CutSet.clear(); CutSet.push_back(0);
							CutSet.push_back(LevelNode[0]);
							for (j = 1; j < L; j++) { CutSet.push_back(LevelNode[j]); }
							CutSet.push_back(LevelNode[L]);
							////std::cout << "\n\nNew algo path:\n";
							  //std::cout << "\n";
							  //for (j = 0; j <= L; j++)
							  //{
							//	std::cout << LevelNode[j] << " ";
							  //}

							 //Check if we have already added path
							PathAlreadyAdded = false;

							if (CutIdx < MyCuts->Size)
								for (j = CutIdx; j <= MyCuts->Size - 1; j++)
									if (MyCuts->CPL[j]->IntListSize == L+1)
										for (t = 0; t <= L+1; t++)
											if (MyCuts->CPL[j]->IntList[t] != LevelNode[L - t + 1]) break;
											else if (t == L) PathAlreadyAdded = true;

							//If Path is not already added
							

							
							if (PathAlreadyAdded) { 
								(LevelIdx[L])++;
								continue;
							}
							//If one way path:

							if (L<=300 & isRoot)
							{
								//Add path contraint if we can to avoid running TSP:							
								
								if (LevelNode[L] > 0 && LevelNode[L] <= n) { //If path is only connected to one charger:
									

									//Try adding fixed edge:
									violation = 0;
									violationFC = 0;

									for (i = 1; i < NoOfEdges; i++)
									{
										if (std::find(CutSet.begin() + 2, CutSet.end(), EdgeTail[i]) != CutSet.end() && std::find(CutSet.begin() + 2, CutSet.end(), EdgeHead[i]) != CutSet.end()) {
											violation += EdgeX[i];
											violationFC += 2 * EdgeX[i];
										}
										else if (LevelNode[0] == EdgeTail[i] && std::find(CutSet.begin() + 2, CutSet.end(), EdgeHead[i]) != CutSet.end()){
											if(EdgeHead[i] == LevelNode[1]) violation += EdgeX[i];
											violationFC += EdgeX[i];
										}
									}


									if (violationFC > L * 2 - 1.99) {
										ExtSet.clear();
										rhs[0] = 1;
										rhs[n + f + 1] = 1;
										//printf("\n%d ", LevelNode[0]);
										for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; }// printf("%d ", LevelNode[j]);}
										ExtSet.insert(ExtSet.begin(), 0);
										//printf("\nStarted single charger:");
										if (RunTSPProblem(rhs, LevelNode, L + 1, "Single_charger", nullptr, true, true) > B && L > 2) {
											for (i = 1; i < n + 1; i++) if (std::find(CutSet.begin(), CutSet.end(), i) == CutSet.end()) ExtSet.push_back(i);
											/*printf("\nWE FOUND IT");*/
											CMGR_AddCnstr(MyCuts, CMGR_FIXED_CHARGER_SINGLE, 0, L + 1, CutSet.data(), L * 2 - 2);
											rhs[0] = 0;
											rhs[n + f + 1] = 0;
											//printf("\Finished single charger:");
											for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
											if (isRoot)	return;
											(LevelIdx[L])++;
											continue;
											//return;

										}
										//printf("\Finished single charger:");
										rhs[0] = 0;
										rhs[n +  f + 1] = 0;
										for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
									}


									if (violation >= L - 0.99)
									{
										//Update edges that should be fixed:
										FixedEdge[0] = LevelNode[1];
										//printf("\n%d ", LevelNode[0]);
										rhs[0] = 1;
										rhs[n + f + 1] = 1;
										for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; }// printf("%d ", LevelNode[j]);}

										/*printf("\nViolation is %f", violation);*/
										//Add fixed charger:
										//printf("\nStarted single edge:");
										if (RunTSPProblem(rhs, LevelNode, L + 1, "Single_edge", FixedEdge, true, true) > B && L > 2) {
											CMGR_AddCnstr(MyCuts, CMGR_FIXED_SINGLEEDGE, 0, L + 1, CutSet.data(), L - 1);
											rhs[0] = 0;
											rhs[n + f + 1] = 0;
											for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }

											if (isRoot)	return;
											//printf("\nAdded like fixed edge");
											(LevelIdx[L])++;
											continue;
										}
										//printf("\nFinished single edge:");
										rhs[0] = 0;
										rhs[n + f + 1] = 0;
										for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
									}


									if (Slack[L] < 0.999) {
										CMGR_AddCnstr(MyCuts, CMGR_CHARGER_PATH, 0, L + 1, CutSet.data(), L - 1);
										//printf("\nAdded like Path");
										//return;
										if (isRoot)	return;
										(LevelIdx[L])++;
										continue;
									}
								}
								else { //If path connects two chargers:

									//Check if cut is violated
									violation = 0;
									violationFC = 0;

									for (i = 1; i < NoOfEdges; i++)
									{
										if (std::find(CutSet.begin() + 2, CutSet.end() - 1, EdgeTail[i]) != CutSet.end() - 1 && std::find(CutSet.begin() + 2, CutSet.end() - 1, EdgeHead[i]) != CutSet.end() - 1) {
											violation += EdgeX[i];
											violationFC += 3 * EdgeX[i];
										}
										else if (EdgeTail[i] == LevelNode[0] && std::find(CutSet.begin() + 2, CutSet.end() - 1, EdgeHead[i]) != CutSet.end() - 1) {
											if (EdgeHead[i] == LevelNode[1]) violation += EdgeX[i];
											violationFC += EdgeX[i];
										}
										else if (std::find(CutSet.begin() + 2, CutSet.end() - 1, EdgeTail[i]) != CutSet.end() - 1 && EdgeHead[i] == LevelNode[L]) {
											if (EdgeTail[i] == LevelNode[L - 1]) violation += EdgeX[i];
											violationFC += EdgeX[i];
										}
									}

									if (violationFC > (L-1)*3-1.99) { //Try to add fixed chargers if violation is found

										rhs[0] = 1;
										rhs[n + f + 1] = 1;
										//printf("\n%d ", LevelNode[0]);
										for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; }// printf("\n%d ", LevelNode[j]);}
										//printf("\n%d ", LevelNode[L]);

										if (RunTSPProblem(rhs, LevelNode, L, "Fixed_chargers", nullptr, true) > B && L > 2) {
											ExtSet.clear();
											ExtSet.insert(ExtSet.begin(), 0);
											for (i = 1; i < n + 1; i++) if (std::find(CutSet.begin(), CutSet.end(), i) == CutSet.end()) ExtSet.push_back(i);
											/*printf("\nWE FOUND IT");*/
											CMGR_AddCnstr(MyCuts, CMGR_FIXED_CHARGERS, 0, L + 1, CutSet.data(),  (L-1)*3 - 2);
											rhs[0] = 0;
											rhs[n + f + 1] = 0;
											for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
											if (isRoot)	return;
											(LevelIdx[L])++;
											continue;
											//return;
										}
										rhs[0] = 0;
										rhs[n + f + 1] = 0;
										for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
										//	//printf("\nAdded like fixed charger");
									}
									
									if (violation > L - 0.99)//Add fixed chargers if violated:
									{
										//Update edges that should be fixed:
										FixedEdge[0] = LevelNode[1];
										//printf("\n%d ", LevelNode[0]);
										FixedEdge[1] = LevelNode[L - 1];
										rhs[0] = 1;
										rhs[n + f + 1] = 1;
										for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; }// printf("%d ", LevelNode[j]);}
										//printf("%d ", LevelNode[L]);


										if (RunTSPProblem(rhs, LevelNode, L, "Fixed_edges", FixedEdge, true) > B && L > 2) {
											CMGR_AddCnstr(MyCuts, CMGR_FIXED_EDGES, 0, L + 1, CutSet.data(), L - 1);
											rhs[0] = 0;
											rhs[n + f + 1] = 0;
											for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
											//std::cout << "\nViolation: " << violation<<" LHS: "<<L-1 << "\n";
											//for (j = 0; j <= L; j++)
											//{
											//	std::cout << LevelNode[j] << " ";
											//}
											//return;
											if (isRoot)	return;
											(LevelIdx[L])++;
											continue;
										}
										rhs[0] = 0;
										rhs[n + f + 1] = 0;
										for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
									}


									if (Slack[L] < 0.999) {
										CMGR_AddCnstr(MyCuts, CMGR_CHARGER_PATH, 0, L + 1, CutSet.data(), L - 1);
										//printf("\nAdded like Path");
										//return;
										if (isRoot)	return;
										(LevelIdx[L])++;
										continue;
									}
								}							// If we have a customer, then continure extending the route, until we reach a charger or the depot 

							}
							else
							{
								//Add path contraint if we can to avoid running TSP:							
								if (Slack[L] < 0.999) {
									CMGR_AddCnstr(MyCuts, CMGR_CHARGER_PATH, 0, L + 1, CutSet.data(), L - 1);
									//printf("\nAdded like Path");
									if (isRoot)	return;
									(LevelIdx[L])++;
									continue;
								}
								//if (LevelNode[L] > 0 && LevelNode[L] <= n) { //If path is only connected to one charger:
								//	//Try adding fixed edge:
								//	//Try adding fixed edge:
								//	violation = 0;
								//	violationFC = 0;

								//	for (i = 1; i < NoOfEdges; i++)
								//	{
								//		if (std::find(CutSet.begin() + 2, CutSet.end(), EdgeTail[i]) != CutSet.end() && std::find(CutSet.begin() + 2, CutSet.end(), EdgeHead[i]) != CutSet.end()) {
								//			violation += EdgeX[i];
								//			violationFC += 2 * EdgeX[i];
								//		}
								//		else if (LevelNode[0] == EdgeTail[i] && std::find(CutSet.begin() + 2, CutSet.end(), EdgeHead[i]) != CutSet.end()) {
								//			if (EdgeHead[i] == LevelNode[1]) violation += EdgeX[i];
								//			violationFC += EdgeX[i];
								//		}
								//	}


								//	if (violationFC > L * 2 - 1.99) {
								//		ExtSet.clear();
								//		rhs[0] = 1;
								//		rhs[n + f + 1] = 1;
								//		//printf("\n%d ", LevelNode[0]);
								//		for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; }// printf("%d ", LevelNode[j]);}
								//		ExtSet.insert(ExtSet.begin(), 0);
								//		//printf("\nStarted single charger:");
								//		if (RunTSPProblem(rhs, LevelNode, L + 1, "Single_charger", nullptr, true, true) > B && L > 2) {
								//			for (i = 1; i < n + 1; i++) if (std::find(CutSet.begin(), CutSet.end(), i) == CutSet.end()) ExtSet.push_back(i);
								//			/*printf("\nWE FOUND IT");*/
								//			CMGR_AddCnstr(MyCuts, CMGR_FIXED_CHARGER_SINGLE, 0, L + 1, CutSet.data(), L * 2 - 2);
								//			rhs[0] = 0;
								//			rhs[n + f + 1] = 0;
								//			//printf("\Finished single charger:");
								//			for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
								//			(LevelIdx[L])++;
								//			return;
								//			//return;
								//			CutCounter++;

								//		}
								//		//printf("\Finished single charger:");
								//		rhs[0] = 0;
								//		rhs[n + f + 1] = 0;
								//		for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
								//	}


								//	if (violation >= L - 0.99)
								//	{
								//		//Update edges that should be fixed:
								//		FixedEdge[0] = LevelNode[1];
								//		//printf("\n%d ", LevelNode[0]);
								//		rhs[0] = 1;
								//		rhs[n + f + 1] = 1;
								//		for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; }// printf("%d ", LevelNode[j]);}

								//		/*printf("\nViolation is %f", violation);*/
								//		//Add fixed charger:
								//		//printf("\nStarted single edge:");
								//		if (RunTSPProblem(rhs, LevelNode, L + 1, "Single_edge", FixedEdge, true, true) > B && L > 2) {
								//			CMGR_AddCnstr(MyCuts, CMGR_FIXED_SINGLEEDGE, 0, L + 1, CutSet.data(), L - 1);
								//			rhs[0] = 0;
								//			rhs[n + f + 1] = 0;
								//			for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }

								//			//return;
								//			//printf("\nAdded like fixed edge");
								//			CutCounter++;
								//			(LevelIdx[L])++;
								//			return;
								//		}
								//		//printf("\nFinished single edge:");
								//		rhs[0] = 0;
								//		rhs[n + f + 1] = 0;
								//		for (j = 1; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
								//	}


							//	}
							//	else { //If path connects two chargers:
							//		//try adding fixed edges
							//		//Else we start checking if fixed charger inequality is violated:
							//		//Check if cut is violated
							//		//Check if cut is violated
							//		violation = 0;
							//		violationFC = 0;

							//		for (i = 1; i < NoOfEdges; i++)
							//		{
							//			if (std::find(CutSet.begin() + 2, CutSet.end() - 1, EdgeTail[i]) != CutSet.end() - 1 && std::find(CutSet.begin() + 2, CutSet.end() - 1, EdgeHead[i]) != CutSet.end() - 1) {
							//				violation += EdgeX[i];
							//				violationFC += 3 * EdgeX[i];
							//			}
							//			else if (EdgeTail[i] == LevelNode[0] && std::find(CutSet.begin() + 2, CutSet.end() - 1, EdgeHead[i]) != CutSet.end() - 1) {
							//				if (EdgeHead[i] == LevelNode[1]) violation += EdgeX[i];
							//				violationFC += EdgeX[i];
							//			}
							//			else if (std::find(CutSet.begin() + 2, CutSet.end() - 1, EdgeTail[i]) != CutSet.end() - 1 && EdgeHead[i] == LevelNode[L]) {
							//				if (EdgeTail[i] == LevelNode[L - 1]) violation += EdgeX[i];
							//				violationFC += EdgeX[i];
							//			}
							//		}

							//		if (violationFC > (L - 1) * 3 - 1.99) { //Try to add fixed chargers if violation is found

							//			rhs[0] = 1;
							//			rhs[n + f + 1] = 1;
							//			//printf("\n%d ", LevelNode[0]);
							//			for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; }// printf("\n%d ", LevelNode[j]);}
							//			//printf("\n%d ", LevelNode[L]);

							//			if (RunTSPProblem(rhs, LevelNode, L, "Fixed_chargers", nullptr, true) > B && L > 2) {
							//				ExtSet.clear();
							//				ExtSet.insert(ExtSet.begin(), 0);
							//				for (i = 1; i < n + 1; i++) if (std::find(CutSet.begin(), CutSet.end(), i) == CutSet.end()) ExtSet.push_back(i);
							//				/*printf("\nWE FOUND IT");*/
							//				CMGR_AddCnstr(MyCuts, CMGR_FIXED_CHARGERS, 0, L + 1, CutSet.data(), (L - 1) * 3 - 2);
							//				rhs[0] = 0;
							//				rhs[n + f + 1] = 0;
							//				for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
							//				CutCounter++;
							//				(LevelIdx[L])++;
							//				return;
							//				//return;
							//			}
							//			rhs[0] = 0;
							//			rhs[n + f + 1] = 0;
							//			for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
							//			//	//printf("\nAdded like fixed charger");
							//		}

							//		if (violation > L - 0.99)//Add fixed chargers if violated:
							//		{
							//			//Update edges that should be fixed:
							//			FixedEdge[0] = LevelNode[1];
							//			//printf("\n%d ", LevelNode[0]);
							//			FixedEdge[1] = LevelNode[L - 1];
							//			rhs[0] = 1;
							//			rhs[n + f + 1] = 1;
							//			for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; }// printf("%d ", LevelNode[j]);}
							//			//printf("%d ", LevelNode[L]);


							//			if (RunTSPProblem(rhs, LevelNode, L, "Fixed_edges", FixedEdge, true) > B && L > 2) {
							//				CMGR_AddCnstr(MyCuts, CMGR_FIXED_EDGES, 0, L + 1, CutSet.data(), L - 1);
							//				rhs[0] = 0;
							//				rhs[n + f + 1] = 0;
							//				for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
							//				//std::cout << "\nViolation: " << violation<<" LHS: "<<L-1 << "\n";
							//				//for (j = 0; j <= L; j++)
							//				//{
							//				//	std::cout << LevelNode[j] << " ";
							//				//}
							//				//return;
							//				CutCounter++;
							//				(LevelIdx[L])++;
							//				return;
							//			}
							//			rhs[0] = 0;
							//			rhs[n + f + 1] = 0;
							//			for (j = 1; j < L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
							//		}

							//	}							// If we have a customer, then continure extending the route, until we reach a charger or the depot 
							}
							

						}
						else if (LevelNode[L] > 0 && LevelNode[L] <= n) {
							AddedToPath[LevelNode[L]] = true; //Add customer to path
							L++; //Go to next level

							LevelIdx[L] = 0; //Set index at next level to 0
							LevelNode[L] = 0; //Set node at next level to 0
						}
					}
					
					(LevelIdx[L])++;
				}

				L--;
				AddedToPath[LevelNode[L]] = false;
				(LevelIdx[L])++;
			} while (L > 1);
		}
	}

	//rhs.end();
}
//Assumes that green zone coordinates is sorted in terms of connecting points and that there should also be placed a line from the first and last coordinate of the green zone:
void MFGVRP_Solver::SetupGreenZone() {
	bool Intersects;
	std::vector < std::pair < Point, Point>> GZLines;
	double xdiff = 0;
	double ydiff = 0;
	Point CPnt1 = { 0,0 }, CPnt2 = { 0,0 }, EntryPoint = { 0,0 }, GZPnt1 = { 0,0 }, GZPnt2 = { 0, 0 };
	//Loop through all green zones
	for (int gz = 0; gz < G; gz++)
	{
		SetupSPP(gz);
		GZLines.clear();

		//Save the lines, that defines the green zone
		for (int j = 0; j < GZPoints[gz].size(); j++)
		{
			if (j == GZPoints[gz].size() - 1) {
				GZPnt1 = { GZPoints[gz][0].first, GZPoints[gz][0].second }, GZPnt2 = { GZPoints[gz][j].first, GZPoints[gz][j].second };
				GZLines.push_back(std::make_pair(GZPnt1, GZPnt2));
				cSPP[j + 1][1] = cSPP[1][j + 1] = sqrt(pow(GZPnt1.x - GZPnt2.x, 2) + pow(GZPnt1.y - GZPnt2.y, 2));
				SPPUBs[j + 1][1] = SPPUBs[1][j + 1] = 1;
			}
			else {
				GZPnt1 = { GZPoints[gz][j].first, GZPoints[gz][j].second }, GZPnt2 = { GZPoints[gz][j + 1].first, GZPoints[gz][j + 1].second };
				GZLines.push_back(std::make_pair(GZPnt1, GZPnt2));
				cSPP[j + 1][j + 2] = cSPP[j + 2][j + 1] = sqrt(pow(GZPnt1.x - GZPnt2.x, 2) + pow(GZPnt1.y - GZPnt2.y, 2));
				SPPUBs[j + 1][j + 2] = SPPUBs[j + 2][j + 1] = 1;
			}
		}

		//Find all customers that are located inside a green zone (Assumes that a customers cannot be part of multiple green zones)
		for (int i = 0; i < n; i++)
		{
			if (pnpoly(GZPoints[gz], CustomersPoints[i])) GreenZoneCustomers.push_back(i + 1);
		}

		//Loop throug all customers
		for (int i = 0; i < n; i++)
		{
			//If the tail customers is not part of a green zone
			if (!std::count(GreenZoneCustomers.begin(), GreenZoneCustomers.end(), i)) {
				//Loop through all customers with index larger than j
				for (int j = i + 1; j < n + 1; j++)
				{
					//If the head customer is not part of the green zone
					if (!std::count(GreenZoneCustomers.begin(), GreenZoneCustomers.end(), j)) {
						//Check if the direct path between j and t enters a green zone:
						CPnt1 = { xCoord[i], yCoord[i] }, CPnt2 = { xCoord[j], yCoord[j] };
						for (int l = 0; l < GZLines.size(); l++)
						{
							//If the direct line between j,t intersects with the green zone:
							if (doIntersect(GZLines[l].first, GZLines[l].second, CPnt1, CPnt2)) {
								//printf("\n\nCustomer %d and %d -> GZ points are %d,%d and %d,%d. Customer points are %d,%d and %d,%d", j, t, GZLines[l].first.x, GZLines[l].first.y, GZLines[l].second.x, GZLines[l].second.y, CPnt1.x, CPnt1.y, CPnt2.x, CPnt2.y);

								//Initially set UB to 1 for all links from the customers to the green zone:
								for (int ii = 1; ii < GZPoints[gz].size() + 1; ii++)
								{
									SPPUBs[0][ii] = SPPUBs[ii][0] = 1;
									SPPUBs[ii][GZPoints[gz].size() + 1] = SPPUBs[GZPoints[gz].size() + 1][ii] = 1;
								}
								//Find valid entry points for SPP:
								for (int p = 0; p < GZPoints[gz].size(); p++)
								{
									EntryPoint = { GZPoints[gz][p].first, GZPoints[gz][p].second };
									//Check if entry link instersects with any green zone border:
									for (int ll = 0; ll < GZLines.size(); ll++)
									{
										/*printf("\nEntry = %d,%d and Line is (%d,%d),(%d,%d)", EntryPoint.x, EntryPoint.y, GZLines[ll].first.x, GZLines[ll].first.y, GZLines[ll].second.x, GZLines[ll].second.y);*/
										//If entry point is not part of the line defining the greenzone
										if ((EntryPoint.x != GZLines[ll].first.x || EntryPoint.y != GZLines[ll].first.y) &&
											(EntryPoint.x != GZLines[ll].second.x || EntryPoint.y != GZLines[ll].second.y)) {
											if (doIntersect(CPnt1, EntryPoint, GZLines[ll].first, GZLines[ll].second)) {
												//printf("\nCustomer %d to (%d,%d) intersects with line (%d,%d),(%d,%d). Therefore, (%d,%d) is not part of the entry points for %d ", j,EntryPoint.x,EntryPoint.y,GZLines[ll].first.x,GZLines[ll].first.y, GZLines[ll].second.x, GZLines[ll].second.y, EntryPoint.x, EntryPoint.y, j);
												SPPUBs[0][p + 1] = SPPUBs[p + 1][0] = 0;
											}
											if (doIntersect(CPnt2, EntryPoint, GZLines[ll].first, GZLines[ll].second)) {
												//printf("\nCustomer %d to (%d,%d) intersects with line (%d,%d),(%d,%d). Therefore, (%d,%d) is not part of the entry points for %d ", t, EntryPoint.x, EntryPoint.y, GZLines[ll].first.x, GZLines[ll].first.y, GZLines[ll].second.x, GZLines[ll].second.y, EntryPoint.x, EntryPoint.y, t);
												SPPUBs[GZPoints[gz].size() + 1][p + 1] = SPPUBs[p + 1][GZPoints[gz].size() + 1] = 0;
											}

										}
										//Set costs for entry points in SPP:
										cSPP[0][p + 1] = cSPP[p + 1][0] = sqrt(pow(CPnt1.x - EntryPoint.x, 2) + pow(ydiff = CPnt1.y - EntryPoint.y, 2));
										cSPP[GZPoints[gz].size() + 1][p + 1] = cSPP[p + 1][GZPoints[gz].size() + 1] = sqrt(pow(CPnt2.x - EntryPoint.x, 2) + pow(ydiff = CPnt2.y - EntryPoint.y, 2));
									}

								}

								RunSPP();
								//HERE WE MUST UPDATE OUR COSTS FOR THE CUSTOMERS DRIVING AROUND THE GREEN ZONE. WE NEED TO ADD ANOTHER COST MATRIX, WHICH ONLY APPLIES TO ICEV's
								printf("\nOriginal cost: %f\tWith Greenzone: %f", c[i][j], SPPCplex.getObjValue());
								break;
							}

						}
					}
				}
			}

		}

		//Reset SPP model:
		ResetSPP();

	}

	//GZLines.end();
}


void MFGVRP_Solver::SetupSPP(int idx) {

	//Setup variables and costs for SPP:
	xSPP = IloVar2DMatrix(env, GZPoints[idx].size() + 2);
	cSPP = IloNum2DMatrix(env, GZPoints[idx].size() + 2);
	SPPUBs = IloNum2DMatrix(env, GZPoints[idx].size() + 2);
	for (int i = 0; i < GZPoints[idx].size() + 2; i++) {
		xSPP[i] = IloNumVarArray(env, GZPoints[idx].size() + 2, 0, 1, ILOINT); //Set up variables
		cSPP[i] = IloNumArray(env, GZPoints[idx].size() + 2);
		SPPUBs[i] = IloNumArray(env, GZPoints[idx].size() + 2);
	}


	//Initializing SPP model:
	SPPModel = IloModel(env);
	SPPCplex = IloCplex(SPPModel);
	SPPObj = IloObjective(env);
	SPPObj = IloMinimize(env);
	SPPModel.add(SPPObj);
	int n = xSPP.getSize();

	SPPCons = IloRangeArray(env);
	IloExpr expr(env);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			expr += cSPP[i][j] * xSPP[i][j];

	expr.clear();

	//Add origin constraints
	for (int i = 0; i < n; i++)
		expr += xSPP[0][i] - xSPP[i][0];

	SPPCons.add(expr == 1);

	expr.clear();
	//Add distination constraints
	for (int i = 0; i < n; i++)
		expr += xSPP[n - 1][i] - xSPP[i][n - 1];

	SPPCons.add(expr == -1);

	//Add flow constraints
	for (int i = 1; i < n - 1; i++)
	{
		expr.clear();
		for (int j = 0; j < n; j++) {
			expr += xSPP[i][j] - xSPP[j][i];
		}

		SPPCons.add(expr == 0);
	}

	SPPModel.add(SPPCons);
	SPPCplex.setOut(env.getNullStream());
	SPPCplex.setWarning(env.getNullStream());
}

void MFGVRP_Solver::RunSPP() {


	//Set obj coefficients
	for (int i = 0; i < xSPP.getSize(); i++)
		SPPObj.setLinearCoefs(xSPP[i], cSPP[i]);

	//Set possible flows:
	for (int i = 0; i < xSPP.getSize(); i++)
		for (int j = 0; j < xSPP[i].getSize(); j++)
			xSPP[i][j].setUB(SPPUBs[i][j]);

	SPPModel.add(SPPObj);
	SPPCplex.solve();
}

void MFGVRP_Solver::ResetSPP() {
	//Remove all data from SPP 
	SPPModel.remove(SPPCons);
	SPPModel.remove(SPPObj);
	SPPCons.clear();
	SPPObj.end();
	SPPCons.end();
	SPPCplex.end();
	SPPModel.end();

	for (int i = 0; i < xSPP.getSize(); i++) { xSPP[i].end(); cSPP[i].end(); SPPUBs.end(); }

	xSPP.end(), cSPP.end(), SPPUBs.end();
}
//THE FOLLOWING FUNCTION IS TAKEN DIRECTLY FROM https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
// Given three collinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool MFGVRP_Solver::onSegment(Point p, Point q, Point r)
{
	if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
		q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y))
		return true;

	return false;
}


//THE FOLLOWING FUNCTION IS TAKEN DIRECTLY FROM https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int MFGVRP_Solver::orientation(Point p, Point q, Point r)
{
	// See https://www.geeksforgeeks.org/orientation-3-ordered-points/
	// for details of below formula.
	int val = (q.y - p.y) * (r.x - q.x) -
		(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0;  // collinear

	return (val > 0) ? 1 : 2; // clock or counterclock wise
}

//THE FOLLOWING FUNCTION IS TAKEN DIRECTLY FROM https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool MFGVRP_Solver::doIntersect(Point p1, Point q1, Point p2, Point q2)
{
	// Find the four orientations needed for general and
	// special cases
	int o1 = orientation(p1, q1, p2);
	int o2 = orientation(p1, q1, q2);
	int o3 = orientation(p2, q2, p1);
	int o4 = orientation(p2, q2, q1);

	// General case
	if (o1 != o2 && o3 != o4)
		return true;

	// Special Cases
	// p1, q1 and p2 are collinear and p2 lies on segment p1q1
	if (o1 == 0 && onSegment(p1, p2, q1)) return true;

	// p1, q1 and q2 are collinear and q2 lies on segment p1q1
	if (o2 == 0 && onSegment(p1, q2, q1)) return true;

	// p2, q2 and p1 are collinear and p1 lies on segment p2q2
	if (o3 == 0 && onSegment(p2, p1, q2)) return true;

	// p2, q2 and q1 are collinear and q1 lies on segment p2q2
	if (o4 == 0 && onSegment(p2, q1, q2)) return true;

	return false; // Doesn't fall in any of the above cases
}


//Export Cplex Model:
void MFGVRP_Solver::ExportCplexModel() {
	Cplex.extract(Model);
	Cplex.exportModel("MY_MODEL_WITH_CUTS.lp");
}

void MFGVRP_Solver::FindViolatedChargerPaths(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<double>& EdgeX, CnstrMgrPointer MyCuts) {
	int t, i, j, root, L, idx;
	IloIntArray rhs(env, n * 2 + 2);
	int CutIdx = MyCuts->Size;
	bool PathAlreadyAdded;
	int prev_root = 0;
	int* FixedEdge = new int[2];
	std::vector<int> CutSet;
	std::vector<int> ExtSet;

	for (i = 0; i < n + f + 1; i++)
	{
		AddedToPath[i] = false;

	}
	//printf("\nEdge tail\tEdge head\tEdge x\n");
	//for (i = 0; i < NoOfEdges; i++)
	//{
	//	printf("\n%d\t%d\t%f", EdgeTail[i], EdgeHead[i], EdgeX[i]);
	//}

	//We loop over all edges.
	for (root = 1; root < NoOfEdges; root++) {

		//If an edge contains a charger or the depot, then we start building a path
		if (EdgeTail[root] == 0 || EdgeTail[root] > n)
		{
			if (EdgeHead[root] == 0 || EdgeHead[root] > n) continue; //If the edge connects the depot and a charger then continue (Assumes reprocessing has taken care of all depot-charger edges that violate the battery capacity)

			AddedToPath[EdgeTail[prev_root]] = false; //Remove previous nodes
			AddedToPath[EdgeHead[prev_root]] = false; //Remove previous nodes 
			AddedToPath[EdgeHead[root]] = true; //Set added to path = 1 for first customer; 
			LevelNode[1] = EdgeHead[root];  //First customer on path set 
			LevelNode[0] = EdgeTail[root]; //First node is the depot/charger
			Esum[1] = g * c[EdgeTail[root]][EdgeHead[root]]; //Set to energy consumption to required energy for travelling from initial charger/depot to first customer.
			Tsum[1] = c[EdgeTail[root]][EdgeHead[root]] + s[EdgeHead[root]]; //Set initial time consumption to time consumption from depot (Should be changed if we do not consider time direct proportional with distance)
			Slack[1] = 1 - EdgeX[root];
			L = 2; //Set level to search = to 2
			double violation = 0;
			double ExtraSlack = 0;
			LevelIdx[L] = 1; //Start the search in the begining 
			LevelNode[L] = 0; //Level 2 node to 0 until we find another node to add to path
			prev_root = root;
			do
			{
				while (LevelIdx[L] < NoOfEdges) {
					/*printf("\nEdge ID %d, Edge tail %d, Edge head %d", LevelIdx[L], EdgeTail[LevelIdx[L]], EdgeHead[LevelIdx[L]]);*/
					//Check if a non explored customer is connected to current customer
					idx = EdgeHead[root];

					if (LevelNode[L - 1] == EdgeTail[LevelIdx[L]] && LevelIdx[L] != prev_root) {
						idx = EdgeHead[LevelIdx[L]];
					}

					if (!AddedToPath[idx])
					{
						LevelNode[L] = idx;

						//If slack is strictly less than 1:
						//if (Slack[L - 1] + (1 - EdgeX[LevelIdx[L]]) < 0.999) {
							/*Violation still possible*/
						ExtraSlack = 1 - EdgeX[LevelIdx[L]];
						for (i = 1; i < NoOfEdges; i++) if (LevelNode[L - 1] == EdgeHead[i] && idx == EdgeTail[i]) { ExtraSlack -= EdgeX[i]; break; }


						if (Slack[L - 1] + ExtraSlack < 0.999) {
							//We update the energy and time sum (Should be updated, if we consider time and energy not direct proportional with distance)
							Esum[L] = Esum[L - 1] + g * c[LevelNode[L - 1]][LevelNode[L]];
							if (LevelNode[L] <= n) Tsum[L] = Tsum[L - 1] + c[LevelNode[L - 1]][LevelNode[L]] + s[LevelNode[L]]; else Tsum[L] = Tsum[L - 1] + c[LevelNode[L - 1]][LevelNode[L]];
							Slack[L] = Slack[L - 1] + ExtraSlack;


							if (Esum[L] > B)
							{ /*Potential violation found*/

								////std::cout << "\n\nNew algo path:\n";

								//std::cout << "\n";
								//printf("Slack: %f nodes: %d\n", Slack[L], L);
								//for (j = 0; j <= L; j++)
								//{
								//	std::cout << LevelNode[j] << " ";
								//}

								////Check if we have already added path
								PathAlreadyAdded = false;

								//if (CutIdx < MyCuts->Size)
								//	for (j = CutIdx; j <= MyCuts->Size - 1; j++)
								//		if (MyCuts->CPL[j]->IntListSize == L+1)
								//			for (t = 1; t <= L+1; t++)
								//				if (MyCuts->CPL[j]->IntList[t] != LevelNode[L - t + 1]) break;
								//				else if (t == L) PathAlreadyAdded = true;


								//CMGR_AddCnstr(MyCuts, CMGR_ENERGY_PATH, 0, L, LevelNode, L - 2);

								//If Path is not already added
								if (!PathAlreadyAdded) {
									//Set RHS for ATSP and initiate set for cut
									CutSet.clear(); CutSet.push_back(0);
									CutSet.push_back(LevelNode[0]);
									//rhs[0] = 1;
									//rhs[n+1] = 1;

									//printf("\nRoute 1: %f  %f\t", Esum[L], Slack[L]);
									//for (i = 0; i <= L; i++)
									//{
									//	printf("%d ", LevelNode[i]);
									//}

									for (j = 1; j < L; j++) {
										CutSet.push_back(LevelNode[j]);
										//rhs[LevelNode[j]] = 1; 
										//rhs[LevelNode[j] + n + 1] = 1; 
									}
									CutSet.push_back(LevelNode[L]);

									//Update edges that should be fixed:
									FixedEdge[0] = LevelNode[1];
									FixedEdge[1] = LevelNode[L - 1];
									//Add set or path violation based on ATSP solution 


									CMGR_AddCnstr(MyCuts, CMGR_CHARGER_PATH, 0, L + 1, CutSet.data(), L - 1);
									//printf("\nAdded like Path");
									delete[] FixedEdge;
									return;


								}
							}
							else if (LevelNode[L] > 0 && LevelNode[L] <= n) // If we have a customer, then continure extending the route, until we reach a charger or the depot 

							{
								AddedToPath[LevelNode[L]] = true; //Add customer to path
								L++; //Go to next level

								LevelIdx[L] = 0; //Set index at next level to 0
								LevelNode[L] = 0; //Set node at next level to 0
							}
						}
					}
					(LevelIdx[L])++;
				}

				L--;
				AddedToPath[LevelNode[L]] = false;
				(LevelIdx[L])++;
			} while (L > 1);
		}
	}

	//rhs.end();
}

void MFGVRP_Solver::FindTimeInfeasiblePaths(int NoOfEdges, std::vector<int>& EdgeTail, std::vector<int>& EdgeHead, std::vector<double>& EdgeX, std::vector<int>& Origin, bool ICEV, CnstrMgrPointer MyCuts,bool isRoot) {
	
	try
	{
		int t, i, j, root, L, idx;
		//IloIntArray rhs(env, n * 2 + 2);
		int CutIdx = MyCuts->Size;
		bool PathAlreadyAdded;
		double ExtraSlack = 0;
		double ExtraSetSlack;
		double ServiceTime = 0;
		double TSPObjVal = 0;
		std::vector<int> LevelNodeVector;
		std::vector<int> Chargers;
		std::vector<int> Customers;
		double RechargingTime;
		double GreedyObj = 0;
		int nRootNodes = 0;
		std::set<int> UniqueChargers;
		bool ConnectionFound = false;
		if (ICEV) nRootNodes = n + 1; else nRootNodes = n + f + 1;
		for (i = 0; i < n + f + 1; i++) AddedToPath[i] = false;



		//std::cout << "\n\nNew algo paths:\n";
		for (root = 1; root < nRootNodes; root++) {
			
			LevelNodeVector.clear(); //Clear vector for new node
			LevelNodeVector.push_back(0);
			LevelNodeVector.push_back(root);
			LevelNodeVector.push_back(0);
			//printf("\n");
			//for (i = 0; i < n + f + 1; i++) printf(" %d",AddedToPath[i]);
			AddedToPath[root - 1] = false; //Remove previous customer to 
			if (root <= n)AddedToPath[root] = true; //Set added to path = 1;
			LevelNode[1] = root;  //First customer on path set 
			LevelNode[0] = 0; //First node is Depot
			Tsum[1] = c[0][root]; //Set initial time consumption to time consumption from depot (Should be changed if we do not consider time direct proportional with distance)
			if (root <= n) Tsum[1] += +s[root];
			Esum[1] = g * c[0][root]; //Set dist to initial charging option
			Slack[1] = 0;
			SlackSet[1] = 0;
			L = 2; //Set level to search = to 2
			LevelIdx[L] = 1; //Start the search in the begining 
			LevelNode[L] = 0; //Level 2 node to 0 until we find another node to add to path
			do
			{
				while (LevelIdx[L] <= NoOfEdges) {
					/*printf("\nEdge ID %d, Edge tail %d, Edge head %d", LevelIdx[L], EdgeTail[LevelIdx[L]], EdgeHead[LevelIdx[L]]);*/
					//Check if a non explored customer is connected to current customer
					idx = root;
					ConnectionFound = false;
					//Check if a connection is found:
					if (LevelNode[L - 1] == EdgeTail[LevelIdx[L]] && EdgeHead[LevelIdx[L]] != 0) { 
						
						if (LevelNode[L - 1] > n && L > 2) { 
							if (Origin[LevelIdx[L]] == Origin[LevelIdx[L - 1]]) { 
								idx = EdgeHead[LevelIdx[L]]; 
								ConnectionFound = true; 
							}
						} else if (EdgeHead[LevelIdx[L]] > n) {
							if (Origin[LevelIdx[L]] == Origin[LevelIdx[L]-1])
							{
								if (!AddedToPath[EdgeHead[LevelIdx[L] - 1]])
								{
									idx = EdgeHead[LevelIdx[L]];
									ConnectionFound = true;
								}
							}
							else if (Origin.size() > LevelIdx[L] + 1)
							{
								if (Origin[LevelIdx[L]] == Origin[LevelIdx[L] + 1]) {
									if (!AddedToPath[EdgeHead[LevelIdx[L] + 1]])
									{
										idx = EdgeHead[LevelIdx[L]];
										ConnectionFound = true;
									}
								}

							}
						}
						else { idx = EdgeHead[LevelIdx[L]]; ConnectionFound = true; }

					}

					if (!AddedToPath[idx] && ConnectionFound)
					{
						LevelNodeVector[L] = idx;
						LevelNode[L] = idx;
						ExtraSetSlack = 1;
						ExtraSlack = 1 - EdgeX[LevelIdx[L]];
						if (LevelNode[L] == EdgeTail[LevelIdx[L] - 1] && LevelNode[L - 1] == EdgeHead[LevelIdx[L] - 1])
							ExtraSlack -= EdgeX[LevelIdx[L] - 1];
						else if (LevelIdx[L] + 1<=NoOfEdges)
							if(LevelNode[L] == EdgeTail[LevelIdx[L] + 1] && LevelNode[L - 1] == EdgeHead[LevelIdx[L] + 1])
								ExtraSlack -= EdgeX[LevelIdx[L] + 1];
						if (isRoot)
						{
							//Add additional contribution to LHS for set expansion
							for (i = 1; i <= NoOfEdges; i++) {
								//if (i == LevelIdx[L]) continue;
								if (EdgeTail[i] == idx) {
									if (std::find(LevelNodeVector.begin() + 1, LevelNodeVector.end(), EdgeHead[i]) != LevelNodeVector.end()) ExtraSetSlack -= EdgeX[i];
								}
								if (EdgeHead[i] == idx) {
									if (std::find(LevelNodeVector.begin() + 1, LevelNodeVector.end(), EdgeTail[i]) != LevelNodeVector.end()) ExtraSetSlack -= EdgeX[i];
								}
							}
						}
						else
						{
							ExtraSetSlack = ExtraSlack;
						}

						//If slack is strictly less than 1:
						if (Slack[L - 1] + ExtraSlack < 0.999 || SlackSet[L - 1] + ExtraSetSlack < 0.999) {
							/*Violation still possible*/


							//We update the energy and time sum (Should be updated, if we consider time and energy not direct proportional with distance)
							Tsum[L] = Tsum[L - 1] + c[LevelNode[L - 1]][LevelNode[L]];
							if (LevelNode[L] <= n)Tsum[L] += +s[LevelNode[L]];
							Esum[L] = Esum[L - 1] + g * c[LevelNode[L - 1]][LevelNode[L]];
							RechargingTime = (Esum[L] + g * c[LevelNode[L]][0] - B) * r;
							Slack[L] = Slack[L - 1] + ExtraSlack;
							SlackSet[L] = SlackSet[L - 1] + ExtraSetSlack;
							if (RechargingTime < 0 || ICEV) RechargingTime = 0;

							// If path is violated
							if (Tsum[L] + c[LevelNode[L]][0] + RechargingTime > T)
							{ /*Violation found*/

								//////std::cout << "\n\nNew algo path:\n";
								//std::cout << "\n";
								//for (j = 0; j <= L; j++)
								//{
								//	std::cout << LevelNode[j] << " ";
								//}

								//Check if we have already added path
								PathAlreadyAdded = false;

								if (CutIdx < MyCuts->Size)
									for (j = CutIdx; j <= MyCuts->Size - 1; j++)
										if (MyCuts->CPL[j]->IntListSize == L)
											for (t = 1; t <= L; t++)
												if (MyCuts->CPL[j]->IntList[t] != LevelNode[L - t + 1]) break;
												else if (t == L) PathAlreadyAdded = true;


								//CMGR_AddCnstr(MyCuts, CMGR_ENERGY_PATH, 0, L, LevelNode, L - 2);

								//If Path is not already added
								if (!PathAlreadyAdded) {

									////std::cout << "\n\nNew algo path:\n";
									//std::cout << "\n";
									//for (j = 0; j <= L; j++)
									//{
									//	std::cout << LevelNode[j] << " ";
									//}

									//CMGR_AddCnstr(MyCuts, CMGR_ENERGY_PATH, 0, L, LevelNode, L - 2);
									//printf("\nSlack for set %f --- Slack for path %f\n Set:\t", SlackSet[L], Slack[L]);
									//for (i = 0; i < LevelNodeVector.size()-1; i++) printf("%d ", LevelNodeVector[i]);
									//printf("%d ", idx);
									//Set RHS for ATSP


									//Add set or path violation based on ATSP solution 

									if (Slack[L] < 0.999) {
										if (ICEV) {

											if (L <= 300 & isRoot)
											{
												GreedyObj = Greedy(LevelNode, L, false, ICEV);
												if (GreedyObj > T)
												{
													ServiceTime = 0;
													for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; ServiceTime += s[LevelNode[j]]; }

													TSPCutOff = T - ServiceTime;
													if (RunTSPProblem(rhs, LevelNode, L + 1, "Time_set", nullptr, false, false, true) + ServiceTime > T)
													{

														CMGR_AddCnstr(MyCuts, CMGR_TIME_INFEASIBLE_SET, 0, L, LevelNode, L - 2);
														for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }

														if (!isRoot) return; else { (LevelIdx[L])++; continue; }
													}
													for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
												}
											}
											CMGR_AddCnstr(MyCuts, CMGR_INFEASIBLE_ICEV_PATH, 0, L, LevelNode, L - 2);
											if (!isRoot) return;
										}
										else {
											if (L <= 300 & isRoot)
											{
												Chargers.clear();
												Customers.clear();
												Chargers.push_back(0);
												Customers.push_back(0);
												UniqueChargers.clear();

												for (j = 1; j <= L; j++) { if (LevelNode[j] <= n) Customers.push_back(LevelNode[j]); else Chargers.push_back(LevelNode[j]); }
												
												UniqueChargers.insert(Chargers.begin(), Chargers.end());
												/*printf("\nSize chargers %d Unique Chargers %d", Chargers.size(), UniqueChargers.size());*/
												GreedyObj = Greedy(LevelNode, L, false, ICEV);
												
												if (GreedyObj > T && UniqueChargers.size()==Chargers.size())
												{
													ServiceTime = 0;

													for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; if (LevelNode[j] <= n) ServiceTime += s[LevelNode[j]]; }

													TSPCutOff = T - ServiceTime;
													TSPObjVal = RunTSPProblem(rhs, LevelNode, L + 1, "Time_EV_set", nullptr, false, false, true);

													RechargingTime = (TSPObjVal / g - B) * r;
													if (RechargingTime < 0 || ICEV) RechargingTime = 0;

													
													
													if (TSPObjVal + ServiceTime > T && Chargers.size() == 1)
													{
														CMGR_AddCnstr(MyCuts, CMGR_TIME_INFEASIBLE_SET, 0, L, LevelNode, L - 2);
														for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }

														if (!isRoot) return; else { (LevelIdx[L])++; continue; }
													}


													if (TSPObjVal + RechargingTime + ServiceTime > T)
													{


														CMGR_AddExtCnstr(MyCuts, CMGR_TIME_INFEASIBLE_EV_SET, 0, Customers.size() - 1, Customers.data(), Chargers.size() - 1, Chargers.data(), L - 2);
														for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }

														if (!isRoot) return; else { (LevelIdx[L])++; continue; }
													}
													for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
												}
											}

											CMGR_AddCnstr(MyCuts, CMGR_INFEASIBLE_EV_PATH, 0, L, LevelNode, L - 2);
											if (!isRoot) return;
										}
									}
									else if (ICEV && L <= 300 && isRoot)
									{
										GreedyObj = Greedy(LevelNode, L, false, ICEV);
										if (GreedyObj > T)
										{
											ServiceTime = 0;
											for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; ServiceTime += s[LevelNode[j]]; }

											TSPCutOff = T - ServiceTime;
											if (RunTSPProblem(rhs, LevelNode, L + 1, "Time_set", nullptr, false, false, true) + ServiceTime > T)
											{

												CMGR_AddCnstr(MyCuts, CMGR_TIME_INFEASIBLE_SET, 0, L, LevelNode, L - 2);
												for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }

												if (!isRoot) return;
											}
											for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
										}

									}
									else if (L <= 300 && isRoot)
									{
										GreedyObj = Greedy(LevelNode, L, false, ICEV);
										Chargers.clear();
										Customers.clear();
										Chargers.push_back(0);
										Customers.push_back(0);
										UniqueChargers.clear();

										for (j = 1; j <= L; j++) { if (LevelNode[j] <= n) Customers.push_back(LevelNode[j]); else  Chargers.push_back(LevelNode[j]); }

										UniqueChargers.insert(Chargers.begin(), Chargers.end());
										
										if (GreedyObj > T && UniqueChargers.size() == Chargers.size())
										{
											ServiceTime = 0;
											for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + f + 1] = 1; if (LevelNode[j] <= n) ServiceTime += s[LevelNode[j]]; }

											TSPCutOff = T - ServiceTime;
											TSPObjVal = RunTSPProblem(rhs, LevelNode, L + 1, "Time_EV_set", nullptr, false, false, true);

											RechargingTime = TSPObjVal / g - B;
											if (RechargingTime < 0 || ICEV) RechargingTime = 0;

											Chargers.clear();
											Customers.clear();
											Chargers.push_back(0);
											Customers.push_back(0);

											for (j = 1; j <= L; j++) { if (LevelNode[j] <= n) Customers.push_back(LevelNode[j]); else Chargers.push_back(LevelNode[j]); }


											if (TSPObjVal + ServiceTime > T && Chargers.size() == 1)
											{
												CMGR_AddCnstr(MyCuts, CMGR_TIME_INFEASIBLE_SET, 0, L, LevelNode, L - 2);
												for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }

												(LevelIdx[L])++; continue;
											}


											if (TSPObjVal + RechargingTime + ServiceTime > T)
											{
												CMGR_AddExtCnstr(MyCuts, CMGR_TIME_INFEASIBLE_EV_SET, 0, Customers.size() - 1, Customers.data(), Chargers.size() - 1, Chargers.data(), L - 2);
												for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }

											}
											for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + f + 1] = 0; }
										}
									}
									/*else {
										for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 1; rhs[LevelNode[j] + n + 1] = 1; }

										if (RunTSPProblem(rhs, LevelNode, L + 1) > B) {
											CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_SET, 0, L, LevelNode, L - 2);
											for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + 1] = 0; }

											return;
										}
										for (j = 0; j <= L; j++) { rhs[LevelNode[j]] = 0; rhs[LevelNode[j] + n + 1] = 0; }

									}*/

									//If ATSP violated, then add set violation
									//Else add path violation


									//Reset RHS for ATSP

								}
							}
							else
							{ /*Path is still feasible, but could be violated when extended*/
								LevelNodeVector[L] = idx;
								if (LevelNode[L] <= n)AddedToPath[LevelNode[L]] = true; //Add customer to path
								L++; //Go to next level

								LevelIdx[L] = 0; //Set index at next level to 0
								LevelNode[L] = 0; //Set node at next level to 0
								LevelNodeVector.push_back(0);
							}
						}
					}
					(LevelIdx[L])++;
				}

				L--;
				LevelNodeVector.pop_back();
				AddedToPath[LevelNode[L]] = false;
				LevelNode[L] = 0;
				(LevelIdx[L])++;
			} while (L > 1);

		}
	}
	catch (const std::exception&)
	{
		std::cerr << "Problem in time";
	}
	

}

void MFGVRP_Solver::AgressiveSetSearch(int NoOfEdges, int* EdgeTail, int* EdgeHead, double* EdgeX, CnstrMgrPointer MyCuts, bool Energy,bool ICEV) {
	int nNodes;
	if (Energy || ICEV) nNodes = n + 1;
	else nNodes = n + f + 1;

	//std::vector<std::pair<double, std::pair<int, int>>> sorted_pairs;

	//for (int i = 1; i < NoOfEdges; i++)
	//{
	//	sorted_pairs.push_back({ EdgeX[i],{EdgeTail[i],EdgeHead[i]} });
	//}
	//std::sort(sorted_pairs.begin(), sorted_pairs.end(), std::greater<std::pair<double, std::pair<int, int>>>());
	//
	//for (int i = 0; i < NoOfEdges-1; i++)
	//{
	//	EdgeX[i+1] = sorted_pairs[i].first;
	//	EdgeTail[i+1] = sorted_pairs[i].second.first;
	//	EdgeHead[i+1] = sorted_pairs[i].second.second;
	//}

	Graph *g1 = new Graph(nNodes,f);
	g1->addEdge(EdgeTail, EdgeHead, EdgeX, NoOfEdges,1);
	g1->ICEV = ICEV;
	g1->Energy = Energy;
	g1->MFGVRP = this;
	g1->MyCuts = MyCuts;
	g1->connectedComponents();
	delete g1;
}


void MFGVRP_Solver::AgressiveSetSearch_FC(int NoOfEdges,int NoOfCharEdges, int* EdgeTail, int* EdgeHead, double* EdgeX,  CnstrMgrPointer MyCuts) {
	Graph* g1 = new Graph(n + 1,f);
	printf("\nChargers");
	for (int i = 0; i < NoOfCharEdges; i++)
	{
		printf("\n%d\t%d\t%f",EdgeTail[i],EdgeHead[i],EdgeX[i]);
	}
	printf("\nDirects");
	for (int i = NoOfCharEdges; i < NoOfEdges; i++)
	{
		printf("\n%d\t%d\t%f", EdgeTail[i], EdgeHead[i], EdgeX[i]);
	}

	g1->addEdge(EdgeTail, EdgeHead, EdgeX, NoOfEdges,NoOfCharEdges);
	g1->MFGVRP = this;
	g1->MyCuts = MyCuts;

	std::vector<int> ChargerTails;
	std::vector<int> ChargerHeads;
	std::vector<double> ChargerVals;
	for (int i = 1; i < NoOfCharEdges; i++)
	{
		ChargerTails.push_back(EdgeTail[i]);
		ChargerHeads.push_back(EdgeHead[i]);
		ChargerVals.push_back(EdgeX[i]);
	}
	g1->ChargerTails = &ChargerTails;
	g1->ChargerHeads = &ChargerHeads;
	g1->ChargerVals = &ChargerVals;

	g1->connectedComponentsFC();
	delete g1;
}


double MFGVRP_Solver::Graph::sumEdgeValues(const std::vector<int>& components) {
	double sum = 0;
	std::list<int>::iterator t;
	std::list<double>::iterator e;
	for (int i = 0; i < components.size(); i++) {
		for (int j = i + 1; j < components.size(); j++) {
			for (t = adj[components[i]].begin(), e = edgeValues[components[i]].begin(); t != adj[components[i]].end(); ++t, ++e)
				if (*t == components[j]) {
					sum += *e;
				}
		}
	}

	return sum;
}

double MFGVRP_Solver::Graph::sumEdgesToChargers(const std::vector<int>& components, bool tail) {
	double sum = 0;

	if (tail) edge1 = 0; else edge2 = 0;
	
	for (int i = 0; i <= MFGVRP->f; i++) {
		sums[i] = 0;
		Customers[i] = 0;
		MinDists[i] = 100000;
	}

	if (tail)
	{
		for (int i = 0; i < ChargerHeads->size(); i++)
		{
			for (int t = 0; t < components.size(); t++)
			{
				if (ChargerHeads->at(i)==components[t])
				{
					if (ChargerTails->at(i)==0)
					{
						sums[ChargerTails->at(i)] += ChargerVals->at(i);
						if (MFGVRP->c[ChargerTails->at(i)][ChargerHeads->at(i)] < MinDists[ChargerTails->at(i)])
						{
							Customers[ChargerTails->at(i)] = ChargerHeads->at(i);
							MinDists[ChargerTails->at(i)] = MFGVRP->c[ChargerTails->at(i)][ChargerHeads->at(i)];
						}
					}
					else
					{
						sums[ChargerTails->at(i) - MFGVRP->n] += ChargerVals->at(i);
						if (MFGVRP->c[ChargerTails->at(i)][ChargerHeads->at(i)] < MinDists[ChargerTails->at(i) - MFGVRP->n])
						{
							Customers[ChargerTails->at(i) - MFGVRP->n] = ChargerHeads->at(i);
							MinDists[ChargerTails->at(i) - MFGVRP->n] = MFGVRP->c[ChargerTails->at(i)][ChargerHeads->at(i)];
						}
					}

					if (ChargerVals->at(i) > sum)
					{
						Charger1_Edge = ChargerTails->at(i);
						Customer1_Edge = ChargerHeads->at(i);
						edge1 = ChargerVals->at(i);
						sum = ChargerVals->at(i);
					}
				}
			}

		}
	}
	else {
		for (int i = 0; i < ChargerTails->size(); i++)
		{
			for (int t = 0; t < components.size(); t++)
			{
				if (ChargerTails->at(i) == components[t])
				{
					if (ChargerHeads->at(i) == 0)
					{
						sums[ChargerHeads->at(i)] += ChargerVals->at(i);
						if (MFGVRP->c[ChargerTails->at(i)][ChargerHeads->at(i)] < MinDists[ChargerHeads->at(i)])
						{
							Customers[ChargerHeads->at(i)] = ChargerTails->at(i);
							MinDists[ChargerHeads->at(i)] = MFGVRP->c[ChargerTails->at(i)][ChargerHeads->at(i)];
						}
					}
					else
					{
						sums[ChargerHeads->at(i) - MFGVRP->n] += ChargerVals->at(i);
						if (MFGVRP->c[ChargerTails->at(i)][ChargerHeads->at(i)] < MinDists[ChargerHeads->at(i) - MFGVRP->n])
						{
							Customers[ChargerHeads->at(i) - MFGVRP->n] = ChargerTails->at(i);
							MinDists[ChargerHeads->at(i) - MFGVRP->n] = MFGVRP->c[ChargerTails->at(i)][ChargerHeads->at(i)];
						}
					}

					if (ChargerVals->at(i) > sum)
					{
						Customer2_Edge = ChargerTails->at(i);
						Charger2_Edge = ChargerHeads->at(i);
						edge2 = ChargerVals->at(i);
						sum = ChargerVals->at(i);
					}
				}
			}
		}
	}
	sum = 0;
	for (int i = 0; i <= MFGVRP->f; i++)
	{
		if (sums[i] > sum) {
			sum = sums[i];
			if (i==0)
				if (tail) Charger1 = 0; else Charger2 = 0;
			else
				if (tail) Charger1 = i+MFGVRP->n; else Charger2 = i + MFGVRP->n;

			if (tail) Customer1 = Customers[i]; else Customer2 = Customers[i];
		}
	}
	return sum;
}


void MFGVRP_Solver::Graph::connectedComponentsFC() {
	LargestViolatedSet.clear();
	MaxSum = 0;
	double Cap = MFGVRP->B;
	MaxGreedy = 0;
	bool* visited = new bool[V];
	std::vector<int> Chargers;
	std::vector<int> Customers;
	double TSPObjVal = 0;
	int* FixedEdge = new int[2];


	for (int v = 1; v < V; v++) {
		for (int v = 0; v < V; v++)
			visited[v] = false;

		ConnectedComponents.clear();
		SetUpdated = false;
		MaxGreedy = 0;
		MaxSum = 0;
		DFSUtilFC(v, visited,0,0);
		if (!SetUpdated) continue;


		//Check violation:
		CharSum1 = sumEdgesToChargers(LargestViolatedSet, true);
		CharSum2 = sumEdgesToChargers(LargestViolatedSet, false);
		if (CharSum1 + MaxSum * 2 > 2 * LargestViolatedSet.size() -2 && MaxGreedy + MFGVRP->c[Charger1][Customer1] + MFGVRP->MinCharDist[LargestViolatedSet.back()] > Cap)
		{
			LargestViolatedSet.insert(LargestViolatedSet.begin(), Charger1);
			//LargestViolatedSet.insert(LargestViolatedSet.begin(), 0);
			//MaxGreedy = MFGVRP->c[Charger1][LargestViolatedSet[0]] + MFGVRP->MinCharDist[LargestViolatedSet.back()];
			OneWay = true; FixedArcs = false;
		}
		else if (MaxSum + edge1 > LargestViolatedSet.size() - 0.99 && MaxGreedy + MFGVRP->c[Charger1_Edge][Customer1_Edge] + MFGVRP->MinCharDist[LargestViolatedSet.back()] > Cap)
		{
			LargestViolatedSet.insert(LargestViolatedSet.begin(), Charger1_Edge);
			//LargestViolatedSet.insert(LargestViolatedSet.begin(), 0);
			//MaxGreedy = MFGVRP->c[Charger1][LargestViolatedSet[0]] + MFGVRP->MinCharDist[LargestViolatedSet.back()];
			OneWay = true; FixedArcs = true;
		}
		else if (CharSum1 + CharSum2 + MaxSum * 3 > 3 * LargestViolatedSet.size()-1.99 && MaxGreedy + MFGVRP->c[Charger1][Customer1] + MFGVRP->c[Customer2][Charger2] > Cap)
		{
			LargestViolatedSet.insert(LargestViolatedSet.begin(), Charger1);
			//LargestViolatedSet.insert(LargestViolatedSet.begin(), 0);
			LargestViolatedSet.push_back(Charger2);
			//MaxGreedy = MFGVRP->c[Charger1][LargestViolatedSet[0]] + MFGVRP->MinCharDist[LargestViolatedSet.back()];
			OneWay = false; FixedArcs = false;
		}
		else if (MaxSum + edge1 + edge2 > LargestViolatedSet.size() + 0.01 && MaxGreedy + MFGVRP->c[Charger1_Edge][Customer1_Edge] + MFGVRP->c[Customer2_Edge][Charger2_Edge] > Cap)
		{
			LargestViolatedSet.insert(LargestViolatedSet.begin(), Charger1_Edge);
			//LargestViolatedSet.insert(LargestViolatedSet.begin(), 0);
			LargestViolatedSet.push_back(Charger2_Edge);
			//MaxGreedy = MFGVRP->c[Charger1][LargestViolatedSet[0]] + MFGVRP->MinCharDist[LargestViolatedSet.back()];
			OneWay = false; FixedArcs = true;
		}
		else continue;


		//Validate cut:
		if (OneWay && FixedArcs)
		{
			//Update edges that should be fixed:
			FixedEdge[0] = Customer1_Edge;
			MFGVRP->rhs[0] = 1;
			MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 1;
			for (int j = 1; j < LargestViolatedSet.size() ; j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 1; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 1; }
			//printf("%d ", LevelNode[L]);

			if (MFGVRP->RunTSPProblem(MFGVRP->rhs, LargestViolatedSet.data(), LargestViolatedSet.size(), "Single_edge", FixedEdge, true,true) > MFGVRP->B) {
				LargestViolatedSet.insert(LargestViolatedSet.begin(), 0);
				CMGR_AddCnstr(MyCuts, CMGR_FIXED_SINGLEEDGE, 0, LargestViolatedSet.size()-1 , LargestViolatedSet.data(), LargestViolatedSet.size() - 3);
				LargestViolatedSet.erase(LargestViolatedSet.begin());
				MFGVRP->rhs[0] = 0;
				MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 0;
				for (int j = 1; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
				delete[] visited;
				delete[] FixedEdge;
				return;

			}
			MFGVRP->rhs[0] = 0;
			MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 0;
			for (int j = 1; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
		}
		else if (!OneWay && FixedArcs)
		{
			//Update edges that should be fixed:
			FixedEdge[0] = Customer1_Edge;
			//printf("\n%d ", LevelNode[0]);
			FixedEdge[1] = Customer2_Edge;
			MFGVRP->rhs[0] = 1; 
			MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 1;
			for (int j = 1; j < LargestViolatedSet.size()-1; j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 1; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 1; }
			//printf("%d ", LevelNode[L]);

			if (MFGVRP->RunTSPProblem(MFGVRP->rhs, LargestViolatedSet.data(), LargestViolatedSet.size()-1, "Fixed_edges", FixedEdge, true) > MFGVRP->B) {
				LargestViolatedSet.insert(LargestViolatedSet.begin(), 0);
				CMGR_AddCnstr(MyCuts, CMGR_FIXED_EDGES, 0, LargestViolatedSet.size()-1, LargestViolatedSet.data(), LargestViolatedSet.size()-3);
				LargestViolatedSet.erase(LargestViolatedSet.begin());
				MFGVRP->rhs[0] = 0;
				MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 0;
				for (int j = 1; j < LargestViolatedSet.size() - 1; j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
				delete[] visited;
				delete[] FixedEdge;
				return;

			}
			MFGVRP->rhs[0] = 0;
			MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 0;
			for (int j = 1; j < LargestViolatedSet.size() - 1; j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
		}
		else if (OneWay && !FixedArcs)
		{
			MFGVRP->rhs[0] = 1;
			MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 1;
			for (int j = 1; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 1; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 1; }

			if (MFGVRP->RunTSPProblem(MFGVRP->rhs, LargestViolatedSet.data(), LargestViolatedSet.size(), "Single_charger", nullptr, true, true) > MFGVRP->B) {
				LargestViolatedSet.insert(LargestViolatedSet.begin(), 0);
				CMGR_AddCnstr(MyCuts, CMGR_FIXED_CHARGER_SINGLE, 0, LargestViolatedSet.size() - 1, LargestViolatedSet.data(), 2*(LargestViolatedSet.size() - 2)-2);
				LargestViolatedSet.erase(LargestViolatedSet.begin());
				MFGVRP->rhs[0] = 0;
				MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 0;
				for (int j = 1; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
				delete[] visited;
				return;

			}
			MFGVRP->rhs[0] = 0;
			MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 0;
			for (int j = 1; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
		}
		else
		{
			MFGVRP->rhs[0] = 1;
			MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 1;
			for (int j = 1; j < LargestViolatedSet.size() - 1; j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 1; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 1; }

			if (MFGVRP->RunTSPProblem(MFGVRP->rhs, LargestViolatedSet.data(), LargestViolatedSet.size() - 1, "Fixed_chargers", nullptr, true) > MFGVRP->B) {
				LargestViolatedSet.insert(LargestViolatedSet.begin(), 0);
				CMGR_AddCnstr(MyCuts, CMGR_FIXED_CHARGERS, 0, LargestViolatedSet.size() - 1, LargestViolatedSet.data(), 3*(LargestViolatedSet.size()-3 ) - 2);
				LargestViolatedSet.erase(LargestViolatedSet.begin());
				MFGVRP->rhs[0] = 0;
				MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 0;
				for (int j = 1; j < LargestViolatedSet.size() - 1; j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
				delete[] visited;
				delete[] FixedEdge;
				return;

			}
			MFGVRP->rhs[0] = 0;
			MFGVRP->rhs[MFGVRP->n + MFGVRP->f + 1] = 0;
			for (int j = 1; j < LargestViolatedSet.size() - 1; j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
		}
	}



	delete[] visited;
	delete[] FixedEdge;
}



// Method to print connected components in an
// undirected graph
void MFGVRP_Solver::Graph::connectedComponents()
{
	LargestViolatedSet.clear();
	MaxSum = 0;
	
	if (Energy) Cap=MFGVRP->B; else Cap=MFGVRP->T;
	
	ConnectedComponents.reserve(MFGVRP->n + MFGVRP->f);
	bool* visited = new bool[V];
	std::vector<int> Chargers;
	std::vector<int> Customers;
	double RechargingTime;
	double TSPObjVal = 0;
	double GreedyVal;
	
	for (int v = 1; v < V; v++) {
		if (visited[v]) continue;
		for (int i = 0; i < V; i++)
			visited[i] = false;
		ConnectedComponents.clear();
		MaxGreedy = 0;
		SetUpdated = false;
		DFSUtil(v, visited,0,0,0);
		if (!SetUpdated) continue;
		//if (LargestViolatedSet.size()>25)
		//{
		//	continue;
		//}
		//if (LargestViolatedSet.size() > 20)
		//{
		//	break;
		//}
		LargestViolatedSet.insert(LargestViolatedSet.begin(), 0);
		GreedyVal=MFGVRP->Greedy(LargestViolatedSet.data(), LargestViolatedSet.size(), Energy, ICEV);
		if (ICEV)
		{
			//IF ICEV
			if (GreedyVal>Cap)
			{
				ServiceTime = 0;
				for (int j = 0; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 1; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 1; ServiceTime += MFGVRP->s[LargestViolatedSet[j]]; }

				MFGVRP->TSPCutOff = Cap - ServiceTime;
				if (MFGVRP->RunTSPProblem(MFGVRP->rhs, LargestViolatedSet.data(), LargestViolatedSet.size(), "Time_set", nullptr, false, false, true) + ServiceTime > Cap)
				{
					CMGR_AddCnstr(MyCuts, CMGR_TIME_INFEASIBLE_SET, 0, LargestViolatedSet.size()-1, LargestViolatedSet.data(), LargestViolatedSet.size() - 3);
					for (int j = 0; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
					delete[] visited;
					return;
				}
				for (int j = 0; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
			}
		}
		else
		{
			if (Energy)
			{
				//IF energy comsumption
				if (GreedyVal > Cap)
				{
					for (int j = 0; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 1; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 1; }
					/*MFGVRP->Greedy(LargestViolatedSet.data(), LargestViolatedSet.size(), Energy);*/
					MFGVRP->TSPCutOff = Cap;
					if (MFGVRP->RunTSPProblem(MFGVRP->rhs, LargestViolatedSet.data(), LargestViolatedSet.size(), "No_charge") > Cap)
					{

						CMGR_AddCnstr(MyCuts, CMGR_NO_CHARGE_SET, 0, LargestViolatedSet.size() - 1, LargestViolatedSet.data(), LargestViolatedSet.size() - 3);
						for (int j = 0; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
						delete[] visited;
						return;
					}
					for (int j = 0; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
				}
			}
			else
			{
				ServiceTime = 0;
				//RechargingTime = (MaxGreedy/MFGVRP->g - Cap) * MFGVRP->r;
				//if (RechargingTime < 0) RechargingTime = 0;

				if (GreedyVal > Cap){
					for (int j = 0; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 1; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 1; if (LargestViolatedSet[j]<=MFGVRP->n) ServiceTime += MFGVRP->s[LargestViolatedSet[j]]; }

					MFGVRP->TSPCutOff = Cap - ServiceTime;
					TSPObjVal = MFGVRP->RunTSPProblem(MFGVRP->rhs, LargestViolatedSet.data(), LargestViolatedSet.size(), "Time_set", nullptr, false, false, true);

					RechargingTime = TSPObjVal / MFGVRP->g - MFGVRP->B;
					if (RechargingTime < 0 ) RechargingTime = 0;

					Chargers.clear();
					Customers.clear();
					Chargers.push_back(0);
					Customers.push_back(0);

					for (int j = 1; j < LargestViolatedSet.size(); j++) { if (LargestViolatedSet[j] <= MFGVRP->n) Customers.push_back(LargestViolatedSet[j]); else Chargers.push_back(LargestViolatedSet[j]); }

					if (TSPObjVal + ServiceTime > Cap && Chargers.size() == 1)
					{
						CMGR_AddCnstr(MyCuts, CMGR_TIME_INFEASIBLE_SET, 0, LargestViolatedSet.size() - 1, LargestViolatedSet.data(), LargestViolatedSet.size() - 3);
						for (int j = 0; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
						delete[] visited;
						return;
					}


					else if (TSPObjVal + RechargingTime + ServiceTime > Cap)
					{
						CMGR_AddExtCnstr(MyCuts, CMGR_TIME_INFEASIBLE_EV_SET, 0, Customers.size() - 1, Customers.data(), Chargers.size() - 1, Chargers.data(), LargestViolatedSet.size() - 3);
						for (int j = 0; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
						delete[] visited;
						return;

					}
					for (int j = 0; j < LargestViolatedSet.size(); j++) { MFGVRP->rhs[LargestViolatedSet[j]] = 0; MFGVRP->rhs[LargestViolatedSet[j] + MFGVRP->n + MFGVRP->f + 1] = 0; }
					}
			}
		}
	}
	delete[] visited;
}


void MFGVRP_Solver::Graph::DFSUtil(int v, bool visited[],double dist,double service,double CurSum)
{
	// Mark the current node as visited and print it
	visited[v] = true;
	//printf(" %d ", v);
	ConnectedComponents.push_back(v);
	//if (ConnectedComponents.size()>25)
	//{
	//	return;
	//}
	//if (ConnectedComponents.size() > 20) return;
	//if (ConnectedComponents.size() > 30) return;
	if (!Energy && v <= MFGVRP->n) service += MFGVRP->s[v];
	
	if (Energy) ExtraDist = MFGVRP->MinCharDist[ConnectedComponents[0]] + MFGVRP->MinCharDist[ConnectedComponents.back()];
	else ExtraDist = MFGVRP->c[0][ConnectedComponents[0]] + MFGVRP->c[ConnectedComponents.back()][0];

	if (service+dist+ExtraDist>Cap)
	{
	//	CurSum = sumEdgeValues(ConnectedComponents);
		if (CurSum > ConnectedComponents.size() - 1.99 && service + dist + ExtraDist > MaxGreedy)
		{                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
			//ConnectedComponents.insert(ConnectedComponents.begin(), 0);
			//GreedyObj = MFGVRP->Greedy(ConnectedComponents.data(), ConnectedComponents.size(), Energy, ICEV);
			//if (service + dist + ExtraDist > MaxGreedy) {
				LargestViolatedSet = ConnectedComponents;
				MaxSum = CurSum;
				MaxGreedy = service + dist + ExtraDist;
				SetUpdated = true;
			//}
			//ConnectedComponents.erase(ConnectedComponents.begin());
		}
	}
	
	
	// Recur for all the vertices
	// adjacent to this vertex
	std::list<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); ++i)
		if (!visited[*i]) {
			ExtraSum = 0;
			for (int j = 0; j < adj[*i].size(); j++)
				if (visited[*std::next(adj[*i].begin(), j)])
					ExtraSum += *std::next(edgeValues[*i].begin(), j);
					
			DFSUtil(*i, visited,dist+MFGVRP->c[v][*i],service,CurSum+ExtraSum);
			// add the corresponding edge value to the sum variable
		}
}
//
void MFGVRP_Solver::Graph::DFSUtilFC(int v, bool visited[],double dist,double CurSum) {
	// Mark the current node as visited and print it
	visited[v] = true;
	FCFound = false;
	//printf(" %d ", v);
	ConnectedComponents.push_back(v);
	//CurSum = sumEdgeValues(ConnectedComponents);
	if (CurSum > ConnectedComponents.size() - 1.99 && dist>MaxGreedy)
	{	
		LargestViolatedSet = ConnectedComponents;
		SetUpdated = true;
		MaxGreedy = dist;
		MaxSum = CurSum;
		//Check connection to one charger:
		//CharSum=sumEdgesToChargers(ConnectedComponents, true);
		////Check if fix charger is violated:
		//if (CurSum*2+CharSum>ConnectedComponents.size()-1.99)
		//{
		//	ConnectedComponents.insert(ConnectedComponents.begin(), Charger1);
		//	GreedyObj = MFGVRP->GreedyFixedChargers(ConnectedComponents.data(),ConnectedComponents.size(),true,false);

		//	if (GreedyObj>MaxGreedy)
		//	{
		//		LargestViolatedSet = ConnectedComponents;
		//		MaxSum = CurSum;
		//		MaxGreedy = GreedyObj;
		//		FCFound=true;
		//		SetUpdated = true;
		//		OneWay = true;
		//		FixedArcs = false;
		//	}

		//	ConnectedComponents.erase(ConnectedComponents.begin());
		//}
		////Check if fixed charging arc is violated:
		//if (!FCFound && CurSum+CharSum>ConnectedComponents.size()-0.99)
		//{
		//	ConnectedComponents.insert(ConnectedComponents.begin(), Charger1);
		//	GreedyObj = MFGVRP->GreedyFixedChargers(ConnectedComponents.data(), ConnectedComponents.size(), true, true);

		//	if (GreedyObj > MaxGreedy)
		//	{
		//		LargestViolatedSet = ConnectedComponents;
		//		MaxSum = CurSum;
		//		MaxGreedy = GreedyObj;
		//		FCFound = true;
		//		SetUpdated = true;
		//		OneWay = true;
		//		FixedArcs = true;
		//	}

		//	ConnectedComponents.erase(ConnectedComponents.begin());
		//}
		//
		//CharSum += sumEdgesToChargers(ConnectedComponents, false);
		////Check connection to two chargers:
		////Check if fix charger is violated:
 
		//{
		//	ConnectedComponents.insert(ConnectedComponents.begin(), Charger1);
		//	ConnectedComponents.push_back(Charger2);
		//	GreedyObj = MFGVRP->GreedyFixedChargers(ConnectedComponents.data(), ConnectedComponents.size()-1, false, false);

		//	if (GreedyObj > MaxGreedy)
		//	{
		//		LargestViolatedSet = ConnectedComponents;
		//		MaxSum = CurSum;
		//		MaxGreedy = GreedyObj;
		//		FCFound = true;
		//		SetUpdated = true;
		//		OneWay = false;
		//		FixedArcs = false;
		//	}

		//	ConnectedComponents.erase(ConnectedComponents.begin());
		//	ConnectedComponents.pop_back();
		//}
		//
		////Check if fix charging arcs is violated:
		//if (!FCFound && CurSum + CharSum > ConnectedComponents.size()-0.99)
		//{
		//	ConnectedComponents.insert(ConnectedComponents.begin(), Charger1);
		//	ConnectedComponents.push_back(Charger2);
		//	GreedyObj = MFGVRP->GreedyFixedChargers(ConnectedComponents.data(), ConnectedComponents.size()-1, false, true);

		//	if (GreedyObj > MaxGreedy)
		//	{
		//		LargestViolatedSet = ConnectedComponents;
		//		MaxSum = CurSum;
		//		MaxGreedy = GreedyObj;
		//		FCFound = true;
		//		SetUpdated = true;
		//		OneWay = false;
		//		FixedArcs = true;
		//	}

		//	ConnectedComponents.erase(ConnectedComponents.begin());
		//	ConnectedComponents.pop_back();
		//}
	}

	// Recur for all the vertices
	// adjacent to this vertex
	std::list<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); ++i)
		if (!visited[*i]) {
			ExtraSum = 0;
			for (int j = 0; j < adj[*i].size(); j++)
				if (visited[*std::next(adj[*i].begin(), j)])
					ExtraSum += *std::next(edgeValues[*i].begin(), j);

			DFSUtilFC(*i, visited, dist + MFGVRP->c[v][*i], CurSum + ExtraSum);
			// add the corresponding edge value to the sum variable
		}
}

MFGVRP_Solver::Graph::Graph(int V,int f)
{
	this->V = V;
	adj = new std::list<int>[V];
	edgeValues = new std::list<double>[V]; // added array to store edge values
	sums = new double[f + 1];
	MinDists = new double[f + 1];
	Customers = new int[f + 1];
}

MFGVRP_Solver::Graph::~Graph() {
	delete[] adj;
	delete[] edgeValues;
	delete[] sums;
	delete[] Customers;
	delete[] MinDists;
}

// method to add an undirected edge
void MFGVRP_Solver::Graph::addEdge(int* v, int* w, double* edgevalue, int m, int start)
{
	for (int i = start; i < m; i++)
	{
		if (v[i]*w[i] == 0) continue;
		adj[v[i]].push_back(w[i]);
		adj[w[i]].push_back(v[i]);
		edgeValues[v[i]].push_back(edgevalue[i]); // added edgevalue to array
		edgeValues[w[i]].push_back(edgevalue[i]);
		
	}

}
double MFGVRP_Solver::Greedy(int* Components,int nCities, bool Energy, bool ICEV) {
	double GreedyVal = 0;
	// Vector to store the solution path
	std::vector<int> path;

	// Set the starting city
	int current_city = 0;
	path.push_back(current_city);
	double dist = 0;
	// Set all cities as unvisited
	std::vector<bool> visited(nCities, false);
	visited[current_city] = true;
	double tot_distance = 0;
	
	// Repeat until all cities have been visited
	while (path.size() < nCities) {
		// Set the next city as the first unvisited city
		int next_city = -1;
		double min_distance = 1e9;

		// Find the closest unvisited city
		for (int i = 0; i < nCities; i++) {
			if (!visited[i]) {
				if (Energy) { if (current_city == 0) dist = g * MinCharDist[Components[i]]; else dist = g * c[Components[current_city]][Components[i]]; }
				else dist = c[current_city][Components[i]];
				
				if (dist < min_distance) {
					next_city = i;
					min_distance = dist;
				}
			}
		}

		// Mark the next city as visited
		visited[next_city] = true;
		if (Components[next_city] <= n) GreedyVal += min_distance + s[Components[next_city]]; else GreedyVal += min_distance;
		tot_distance += min_distance;
		// Add the next city to the solution path
		path.push_back(Components[next_city]);

		// Set the next city as the current city
		current_city = next_city;
	}
	if (Energy) GreedyVal += g * MinCharDist[path.back()];
	else GreedyVal += c[path.back()][0];

	
	if (!ICEV & !Energy & tot_distance / g > B) GreedyVal += r * (tot_distance / g - B);
	return GreedyVal;
}

double MFGVRP_Solver::GreedyFixedChargers(int* Components, int nCities, bool OneWay, bool FixedArc) {
	double GreedyVal = 0;
	// Vector to store the solution path
	std::vector<int> path;

	// Set the starting city
	int current_city = 0;
	path.push_back(Components[current_city]);
	double dist = 0;
	// Set all cities as unvisited
	std::vector<bool> visited(nCities, false);
	visited[current_city] = true;
	int start_city = 1;
	if (FixedArc)
	{	
		start_city = 2;
		GreedyVal = g * c[Components[0]][Components[1]];
		path.push_back(Components[1]);
		current_city = 1;
		visited[current_city] = true;

	}
	// Repeat until all cities have been visited
	while (path.size() < nCities) {
		// Set the next city as the first unvisited city
		int next_city = -1;
		double min_distance = 1e9;

		// Find the closest unvisited city
		for (int i = start_city; i < nCities; i++) {
			if (!visited[i]) {
				dist = g * c[Components[current_city]][Components[i]]; 


				if (dist < min_distance) {
					next_city = i;
					min_distance = dist;
				}
			}
		}

		// Mark the next city as visited
		visited[next_city] = true;
		GreedyVal += min_distance;
		// Add the next city to the solution path
		path.push_back(Components[next_city]);

		// Set the next city as the current city
		current_city = next_city;
	}
	if (OneWay) GreedyVal += g * MinCharDist[path.back()];
	else GreedyVal += c[path.back()][Components[nCities]];


	return GreedyVal;
}

void MFGVRP_Solver::AddResultsToFile() {
	std::ofstream outfile("file.txt", std::ios_base::app);

	std::string name = Stats->instanceName;
	name.erase(name.find("."), 4);
	std::size_t PathID = name.find("\\");
	if (PathID != std::string::npos) name.erase(0, PathID + 1);
	std::ifstream out("file.txt", std::ios::app);
	//outfile << "Name" << "\t" << "n" << "\t" << "r" << "\t" << "EVs" <<
	//	"\t" << "ICEVs" << "\t" << "Upper bound" << "\t" << "Lower bound" << "\t" << "Time" <<
	//	"\t" << "Optimal" << "\t" << "Infeasible" << "\t" << "LB root" <<
	//	"\t" << "Time root" << "\t" << "Cap root" << "\t" <<
	//	"FCI root" << "\t" << "MS root" << "\t" << "Comb root" <<
	//	"\t" << "NCP root" << "\t" << "NCS root" << "\t" << "CP root"
	//	<< "\t" << "FC root" << "\t" << "FCS root"
	//	<< "\t" << "FE root" << "\t" << "FSE root"
	//	<< "\t" << "ICEV_T" << "\t" << "EV_T" << "\t" << "time_set" << "\t" << "time_EV_set"
	//	<< "\t" << "Tree size" << "\t" << "Max depth" << "\t" << "UB depth" << "\t" << "Feas depth"
	//	<< "\t" << "CAP" << "\t" <<
	//	"FCI" << "\t" << "MS" << "\t" << "Comb" <<
	//	"\t" << "NCP" << "\t" << "NCS" << "\t" << "CP"
	//	<< "\t" << "FC" << "\t" << "FCS"
	//	<< "\t" << "FE" << "\t" << "FSE"  << "\t" << "ICEV_T" <<"\t" << "EV_T" <<"\t" << "time_set" << "\t" << "time_EV_set" << "\t" << "Num TSP"
	//	<< "\t" << "Avg TSP" << "\t" << "Max TSP" << "\t" << "Min TSP" "\t" << "Total TSP time"<< "\t" << "Total Cap time" << "\t" <<"Total enumeration time" <<"\t" << "Total Sep time" << "\t" << "Fetch data time" << "\t" << "Add cuts time" << "\t" << "Convert vector time" << "\t" << "Seperation Vector time"
	//	<< std::endl;

	outfile << name << "\t" << Stats->n << "\t" << Stats->r << "\t" << Stats->EVs <<
		"\t" << Stats->ICEVs << "\t" << Stats->bestUpperBound << "\t" << Stats->bestLowerBound << "\t" << Stats->time <<
		"\t" << Stats->isOptimal << "\t" << Stats->isInfeasible << "\t" << Stats->LBAtRootNode <<
		"\t" << Stats->timeRootNode << "\t" << Stats->numberOfCutsAtRootNode->CapCuts << "\t" <<
		Stats->numberOfCutsAtRootNode->FCI << "\t" << Stats->numberOfCutsAtRootNode->MultiStar << "\t" << Stats->numberOfCutsAtRootNode->Comb <<
		"\t" << Stats->numberOfCutsAtRootNode->NoChargePath << "\t" << Stats->numberOfCutsAtRootNode->NoChargeSet << "\t" << Stats->numberOfCutsAtRootNode->ChargerPath
		<< "\t" << Stats->numberOfCutsAtRootNode->FixedChargers << "\t" << Stats->numberOfCutsAtRootNode->FixedChargerSingle
		<< "\t" << Stats->numberOfCutsAtRootNode->FixedEdges << "\t" << Stats->numberOfCutsAtRootNode->FixedSingleEdge
		<< "\t" << Stats->numberOfCutsAtRootNode->TimeInfeasibleICEVPath << "\t" << Stats->numberOfCutsAtRootNode->TimeInfeasibleEVPath << "\t" << Stats->numberOfCutsAtRootNode->TimeInfeasibleSet << "\t" << Stats->numberOfCutsAtRootNode->TimeInfeasibleEVSet
		<< "\t" << Stats->treeSize << "\t" << Stats->treeDepth << "\t" << Stats->treeUBDepth << "\t" << Stats->treeFeasDepth
		<< "\t" << Stats->totalNumberOfCuts->CapCuts << "\t" <<
		Stats->totalNumberOfCuts->FCI << "\t" << Stats->totalNumberOfCuts->MultiStar << "\t" << Stats->totalNumberOfCuts->Comb <<
		"\t" << Stats->totalNumberOfCuts->NoChargePath << "\t" << Stats->totalNumberOfCuts->NoChargeSet << "\t" << Stats->totalNumberOfCuts->ChargerPath
		<< "\t" << Stats->totalNumberOfCuts->FixedChargers << "\t" << Stats->totalNumberOfCuts->FixedChargerSingle
		<< "\t" << Stats->totalNumberOfCuts->FixedEdges << "\t" << Stats->totalNumberOfCuts->FixedSingleEdge << "\t" << Stats->totalNumberOfCuts->TimeInfeasibleICEVPath <<"\t" << Stats->totalNumberOfCuts->TimeInfeasibleEVPath << "\t" << Stats->totalNumberOfCuts->TimeInfeasibleSet << "\t" << Stats->totalNumberOfCuts->TimeInfeasibleEVSet << "\t" << Stats->TSPSolved
		<< "\t" << Stats->avgTSPTime << "\t" << Stats->maxTSPTime << "\t" << Stats->minTSPtime << "\t" << Stats->TSPSolved * Stats->avgTSPTime << "\t" << Stats->CapSepTime <<"\t" <<Stats->EnumerationTime- Stats->TSPSolved * Stats->avgTSPTime << "\t" << Stats->TotalSepTime << "\t" << Stats->TotalLoadSepData << "\t" << Stats->AddCutsTime << "\t" << Stats->ConvertVals << "\t" << Stats->SetupSepVectors
		<< std::endl;

	outfile.close();

	//Solution file:

	if (Stats->TSPSolved > 0) {
		std::ofstream TSPfile("TSP/" + name + "_TSP.json");

		TSPfile << "{" << std::endl;
		TSPfile << "\t\"Instance\":\t\"" << name << "\"," << std::endl;
		TSPfile << "\t\"noOfTSP\":\t" << Stats->TSPSolved << "," << std::endl;

		TSPfile << "\t\"inqType\":\t[\"";
		for (int i = 0; i < Stats->InequalityType.size() - 1; i++) { TSPfile << Stats->InequalityType[i] << "\", \""; }
		TSPfile << Stats->InequalityType[Stats->InequalityType.size() - 1] << "\"]," << std::endl;

		TSPfile << "\t\"numOfNodes\":\t[";
		for (int i = 0; i < Stats->TSPNodes.size() - 1; i++) { TSPfile << Stats->TSPNodes[i] << ", "; }
		TSPfile << Stats->TSPNodes[Stats->TSPNodes.size() - 1] << "]," << std::endl;

		TSPfile << "\t\"compTime\":\t[";
		for (int i = 0; i < Stats->TSPTime.size() - 1; i++) { TSPfile << Stats->TSPTime[i] << ", "; }
		TSPfile << Stats->TSPTime[Stats->TSPTime.size() - 1] << "]" << std::endl;

		TSPfile << "}";
		TSPfile.close();
	}

	if (Stats->isOptimal)
	{
		IloNum3DMatrix xTemp(env, n+1);
		IloNum3DMatrix yTemp(env, n+1);


		IloNumArray xTemp1d(env, xDummy.getSize());
		Cplex.getValues(xTemp1d, xDummy);

		IloNumArray yTemp1d(env, yDummy.getSize());
		Cplex.getValues(yTemp1d, yDummy);


		int cnt = 0;

		for (int i = 0; i < n+1; i++)
		{
			xTemp[i] = IloNum2DMatrix(env, n+1);
			for (int j = 0; j < n+1; j++)
			{
				xTemp[i][j] = IloNumArray(env, 2);
				for (int k = 0; k < 2; k++)
				{
					xTemp[i][j][k] = xTemp1d[cnt];
					cnt++;
				}
			}
		}

		xTemp1d.end();
		cnt = 0;

		for (int i = 0; i < n+1; i++)
		{
			yTemp[i] = IloNum2DMatrix(env, n+1);
			for (int j = 0; j < n+1; j++)
			{
				yTemp[i][j] = IloNumArray(env, R[i][j].size());
				for (int r = 0; r < R[i][j].size(); r++)
				{
					yTemp[i][j][r] = yTemp1d[cnt];
					cnt++;
				}
			}
		}
		yTemp1d.end();


		std::ofstream Solfile("Solutions/" + name + "_sol.json");
		int CurNode = 0;
		int NoOfRoutes = 0;
		Solfile << "{" << std::endl;
		Solfile << "\t\"Instance\":\t\"" << name << "\"," << std::endl;
		Solfile << "\t\"ObjVal\":\t" << Stats->bestUpperBound << "," << std::endl;
		Solfile << "\t\"noOfRoutesEVRoutes\":\t" << Stats->EVs << "," << std::endl;
		Solfile << "\t\"noOfRoutesICEVRoutes\":\t" << Stats->ICEVs << "," << std::endl;
		Solfile << "\t\"EVRoutes\":\t[";
		for (int i = 0; i < n + 1; i++)
		{
			if (xTemp[CurNode][i][0] > 0.9) {
				Solfile << "[" << CurNode << ", " << i;
				CurNode = i;
				for (int j = 0; j < n + 1; j++)
				{
					if (xTemp[CurNode][j][0] > 0.9) {
						Solfile << ", " << j;
						CurNode = j;
						j = -1;
					}
					if (j != -1)
					{
						for (int r = 0; r < R[CurNode][j].size(); r++)
						{
							if (yTemp[CurNode][j][r] > 0.9)
							{
								Solfile << ", " << R[CurNode][j][r] << ", " << j;
								CurNode = j;
								j = -1;
								break;
							}
						}
					}

					if (CurNode == 0)
					{
						NoOfRoutes++;
						if (NoOfRoutes < Stats->EVs) {
							Solfile << "], ";
							break;
						}
						else
						{
							Solfile << "]";
							break;
						}


					}
				}

			}
			for (int r1 = 0; r1 < R[CurNode][i].size(); r1++)
			{
				if (yTemp[CurNode][i][r1] > 0.9) {
					Solfile << "[" << CurNode << ", " << R[CurNode][i][r1] << ", " << i;
					CurNode = i;
					for (int j = 0; j < n + 1; j++)
					{
						if (xTemp[CurNode][j][0] > 0.9) {
							Solfile << ", " << j;
							CurNode = j;
							j = -1;
						}
						if (j != -1) {
							for (int r = 0; r < R[CurNode][j].size(); r++)
							{
								if (yTemp[CurNode][j][r] > 0.9)
								{
									Solfile << ", " << R[CurNode][j][r] << ", " << j;
									CurNode = j;
									j = -1;
									break;
								}
							}
						}

						if (CurNode == 0)
						{
							NoOfRoutes++;
							if (NoOfRoutes < Stats->EVs) {
								Solfile << "], ";
								break;
							}
							else
							{
								Solfile << "]";
								break;
							}
						}
					}
				}
			}



		}
		Solfile << "]," << std::endl;
		Solfile << "\t\"ICEVRoutes\":\t[";
		NoOfRoutes = 0;
		for (int i = 0; i < n + 1; i++)
		{
			if (xTemp[CurNode][i][1] > 0.9) {
				Solfile << "[" << CurNode << ", " << i;
				CurNode = i;
				for (int j = 0; j < n + 1; j++)
				{
					if (xTemp[CurNode][j][1] > 0.9) {
						Solfile << ", " << j;
						CurNode = j;
						j = -1;
					}

					if (CurNode == 0)
					{
						NoOfRoutes++;
						if (NoOfRoutes < Stats->ICEVs) {
							Solfile << "], ";
							break;
						}
						else
						{
							Solfile << "]";
							break;
						}


					}
				}

			}
		}
		Solfile << "]" << std::endl;
		Solfile << "}";
		Solfile.close();

		for (int i = 0; i < n+1; i++)
		{
			for (int j = 0; j < n+1; j++)
			{
				xTemp[i][j].end();
			}
			xTemp[i].end();
		}
		xTemp.end();

		for (int i = 0; i < n+1; i++)
		{
			for (int j = 0; j < n+1; j++)
			{
				yTemp[i][j].end();
			}
			yTemp[i].end();
		}
		yTemp.end();
	}
	//outfile << Stats->instanceName << "\t" << Stats->n << "\t" << Stats->r << "\t" << Stats->EVs <<
	//	"\t" << Stats->ICEVs << "\t" << Stats->bestUpperBound << "\t" << Stats->bestLowerBound << "\t" << Stats->time <<
	//	"\t" << Stats->isOptimal << "\t" << Stats->isInfeasible << "\t" << Stats->LBAtRootNode <<
	//	"\t" << Stats->timeRootNode << "\t" << Stats->numberOfCutsAtRootNode->CapCuts << "\t" <<
	//	Stats->numberOfCutsAtRootNode->FCI << "\t" << Stats->numberOfCutsAtRootNode->MultiStar <<
	//	"\t" << Stats->numberOfCutsAtRootNode->NoChargePath << "\t" << Stats->numberOfCutsAtRootNode->NoChargeSet << "\t" << Stats->numberOfCutsAtRootNode->ChargerPath
	//	<< "\t" << Stats->numberOfCutsAtRootNode->FixedChargers << "\t" << Stats->numberOfCutsAtRootNode->FixedChargerSingle
	//	<< "\t" << Stats->numberOfCutsAtRootNode->FixedEdges << "\t" << Stats->numberOfCutsAtRootNode->FixedSingleEdge
	//	<< "\t" << Stats->treeSize << "\t" << Stats->treeDepth << "\t" << Stats->treeUBDepth << "\t" << Stats->treeFeasDepth
	//	<< "\t" << Stats->totalNumberOfCuts->CapCuts << "\t" <<
	//	Stats->totalNumberOfCuts->FCI << "\t" << Stats->totalNumberOfCuts->MultiStar <<
	//	"\t" << Stats->totalNumberOfCuts->NoChargePath << "\t" << Stats->totalNumberOfCuts->NoChargeSet << "\t" << Stats->totalNumberOfCuts->ChargerPath
	//	<< "\t" << Stats->totalNumberOfCuts->FixedChargers << "\t" << Stats->totalNumberOfCuts->FixedChargerSingle
	//	<< "\t" << Stats->totalNumberOfCuts->FixedEdges << "\t" << Stats->totalNumberOfCuts->FixedSingleEdge << "\t" << Stats->TSPSolved
	//	<< "\t" << Stats->avgTSPTime << "\t" << Stats->maxTSPTime << Stats->minTSPtime
	//	<< std::endl;

}