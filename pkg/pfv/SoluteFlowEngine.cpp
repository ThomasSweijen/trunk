/*************************************************************************
*  Copyright (C) 2013 by T. Sweijen (T.sweijen@uu.nl)                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef YADE_CGAL
#ifdef FLOW_ENGINE

#define SOLUTE_FLOW
#ifdef SOLUTE_FLOW

#define TEMPLATE_FLOW_NAME SoluteFlowEngineT
#include <yade/pkg/pfv/FlowEngine.hpp>
#include<Eigen/SparseLU>


#undef TEMPLATE_FLOW_NAME

class SoluteCellInfo : public FlowCellInfo
{	
	public:
	Real solute_concentration;
	SoluteCellInfo (void) : FlowCellInfo() {solute_concentration=0;}
	inline Real& solute (void) {return solute_concentration;}
	inline const Real& solute (void) const {return solute_concentration;}
	inline void getInfo (const SoluteCellInfo& otherCellInfo) {FlowCellInfo::getInfo(otherCellInfo); solute()=otherCellInfo.solute();}
	
};

typedef TemplateFlowEngine<SoluteCellInfo,FlowVertexInfo> SoluteFlowEngineT;
REGISTER_SERIALIZABLE(SoluteFlowEngineT);
YADE_PLUGIN((SoluteFlowEngineT));

class SoluteFlowEngine : public SoluteFlowEngineT
{
	public :
		
	Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor, int>, Eigen::AMDOrdering<int> > eSolver2;
		
		void SoluteAction();
		void solveConcentration();
		void InitializeSoluteTransport();
		void solveSoluteTransportMatrix ();
		double getConcentration(unsigned int id){return solver->T[solver->currentTes].cellHandles[id]->info().solute();	}
		double insertConcentration(unsigned int id,double conc){
			solver->T[solver->currentTes].cellHandles[id]->info().solute() = conc;
			return conc;}
		void soluteBC();
		double getConcentrationPlane (double Yobs,double Yr, int xyz);
		double checkMassBalance();
		double massIntegrationStorage();
		void massIntegrationBC();
		///Elaborate the description as you wish
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(SoluteFlowEngine,SoluteFlowEngineT,"A variant of :yref:`FlowEngine` with solute transport).",
		///No additional variable yet, else input here
// 		((Vector3r, gradP, Vector3r::Zero(),,"Macroscopic pressure gradient"))
       		((double,D,0.0,,"Diffusion coefficient, if 0 no diffusion is considered"))
		((unsigned int,bcid1,3,,"Boundary condition 1, enter the O.bodies identifyer of wall, at the upstream side"))
		((unsigned int,bcid2,2,,"Boundary condition 2, enter the O.bodies identifyer of wall (only required for diffusion)"))
		((double,bcconcentration1,1.0,,"Concentration at boundary condition 1"))
		((double,bcconcentration2,0.0,,"Concentration at boundary condition 2"))
		((double,initialConcentration,0.0,,"The initial solute concentration found over the whole domain"))
		((bool,MassStorage,false,,"How to store your data per pore, in mass (true) or volume (false)"))
		((double,totalmassBCin,0.0,,"total mass which entered via the boundary conditions."))
		((double,totalmassBCout,0.0,,"total mass which left via the boundary conditions."))
		((double,massinBC,0.0,,"Mass present in BC pores"))
		((bool,MassIntegration,true,,"Switch for keeping track of mass flux into the system"))
		,,
		
		,
		.def("checkMassBalance",&SoluteFlowEngine::checkMassBalance,"Returns the ratio of mass in system - mass went in via boundary over mass in system in %")
		.def("massIntegrationStorage",&SoluteFlowEngine::massIntegrationStorage,"Mass balance check: integration of mass in the system")
		.def("massIntegrationBC",&SoluteFlowEngine::massIntegrationBC,"Mass balance check: integration of mass in the system")
		.def("solveSoluteTransportMatrix",&SoluteFlowEngine::solveSoluteTransportMatrix,"Solute transport (advection and diffusion) engine for diffusion use a diffusion coefficient (D) other than 0.")
		.def("getConcentration",&SoluteFlowEngine::getConcentration,(python::arg("id")),"get concentration of pore with ID")
		.def("insertConcentration",&SoluteFlowEngine::insertConcentration,(python::arg("id"),python::arg("conc")),"Insert Concentration (ID, Concentration)")
		.def("solute_BC",&SoluteFlowEngine::soluteBC,"Run Boundary Conditions")
		.def("InitializeSoluteTransport",&SoluteFlowEngine::InitializeSoluteTransport,"Reset all concentrations in cells")
		.def("solveConcentration",&SoluteFlowEngine::solveConcentration,"runSolver without recomputing")
		.def("SoluteAction",&SoluteFlowEngine::SoluteAction,"runSolver without recomputing")
		.def("getConcentrationPlane",&SoluteFlowEngine::getConcentrationPlane,(python::arg("Yobs"),python::arg("Yr"),python::arg("xyz")),"Measure concentration within a plane. Concentrations get weighted to their distance to the plane")
		)
};
REGISTER_SERIALIZABLE(SoluteFlowEngine);


// PeriodicFlowEngine::~PeriodicFlowEngine(){}


void SoluteFlowEngine::SoluteAction(){
  bool localDebug = false;

  
  //If this is the first computation, some preparation has to be done
  if (firstSoluteEngine){
      InitializeSoluteTransport();
      solveSoluteTransportMatrix();
      firstSoluteEngine=false;
      if(localDebug){cerr<<"First time in solute Engine!"<<endl;}
  }
  
  
  //After triangulation, matrix has to be resolved
  if(updateSoluteEngine){
    updateVolumes ( *solver );	//NOTE This function updates info()->volume(), and info()->dv(), which is required for getting info()->invvoid()
    solveSoluteTransportMatrix();
    updateSoluteEngine=false;
    if(localDebug){cerr<<"Matrix recomputed due to triangulation"<<endl;}
  }
  
  //Normal calculation of concentration
  soluteBC();
  if(MassIntegration){massIntegrationBC();}
  solveConcentration();
  if(localDebug){cerr<<"Solved for concentration!"<<endl;}
      
      
    
}





void SoluteFlowEngine::InitializeSoluteTransport ()
{
	//Fill the concentration list with 0.0
	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{
	 cell->info().solute() = initialConcentration;
	}
}

void SoluteFlowEngine::solveConcentration(){
	//Solve for concentration using the already solved coefficient matrix
	int ncells=0;
	ncells=solver->T[solver->currentTes].cellHandles.size();
	Eigen::VectorXd eb2(ncells); Eigen::VectorXd ex2(ncells);
  	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
	eb2[cell->info().id]=cell->info().solute();}
	
	ex2 = eSolver2.solve(eb2);
    
	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){ 
	cell->info().solute()= ex2[cell->info().id]; 	}
}

void SoluteFlowEngine::solveSoluteTransportMatrix()

{	
	//In this function the coefficient matrix is solved for advection and diffusion
	//Definitions used for filling the matrix
	double coeff = 0.00;
	double coeff1 = 0.00;	//Coefficient for off-diagonal element
	double coeff2 = 0.00;  //Coefficient for diagonal element
	double qin = 0.00;	//Flux into the pore per pore throat 
	double Qout=0.0;	//Total Flux out of the pore
	int ncells;		//Number of cells
	int ID = 0;
	double invdistance = 0.0;	//Fluid facet area divided by pore throat length for each pore throat
	double invdistancelocal = 0.0; //Sum of invdistance
	ncells=solver->T[solver->currentTes].cellHandles.size();

	//Definitions for solving the matrix
	Eigen::SparseMatrix<double, Eigen::ColMajor,int> Aconc;
	typedef Eigen::Triplet<double> ETriplet2;
	std::vector<ETriplet2> tripletList2;
	Aconc.resize(ncells,ncells);

		
	// Fill coefficient matrix
	
	//Loop over all pores
	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
	cell->info().invVoidVolume() = 1.0 / ( abs(cell->info().volume()) - abs(solver->volumeSolidPore(cell) ) );
	if (cell->info().id == 10){cerr<<endl<<cell->info().invVoidVolume();}
		invdistance=0.0;
	

		//loop over all neighbors of a pore
		unsigned int i=cell->info().id;
		for (unsigned int ngb=0;ngb<4;ngb++){
	  
	  //Calculate distance, facet area between two pores
		      CGT::Point& p2 = cell ->neighbor(ngb)->info();
		      CGT::Point& p1 = cell->info();
		      CGT::CVector l = p1-p2;
		      CGT::Real fluidSurf = sqrt(cell->info().facetSurfaces[ngb].squared_length())*cell->info().facetFluidSurfacesRatio[ngb];
		      invdistancelocal = (fluidSurf/sqrt(l.squared_length()));
		      invdistance+=(fluidSurf/sqrt(l.squared_length()));
		      
		      //Fill off-diagonal coefficients
		      coeff = scene->dt*cell->info().invVoidVolume();
		      ID = cell->neighbor(ngb)->info().id;
		      qin=abs(cell->info().kNorm() [ngb])* ( cell->neighbor ( ngb )->info().p()-cell->info().p());
		      Qout=Qout+max(qin,0.0);
		      coeff1=-1*coeff*(abs(max(qin,0.0))-(D*invdistancelocal));
		      if (coeff1 != 0.0){
			  tripletList2.push_back(ETriplet2(i,ID,coeff1));
		      }
	      }
	      //Fill diagonal coefficients
	      coeff2=1.0+(cell->info().dv()*scene->dt*cell->info().invVoidVolume())+coeff*(abs(Qout)+(D*invdistance));
	      if (coeff2 != 0.0){
		tripletList2.push_back(ETriplet2(i,i,coeff2));
	      }

	      Qout=0.0;
	      qin=0.0;
	}   
	
	//Solve the coefficient Matrix
	Aconc.setFromTriplets(tripletList2.begin(), tripletList2.end());
	eSolver2.compute(Aconc);
	tripletList2.clear();
  }

  
void SoluteFlowEngine::soluteBC()
{
	//Boundary conditions according to soluteTransport.
	//It simply assigns boundary concentrations to cells with a common vertices (e.g. infinite large sphere which makes up the boundary condition in flowEngine)
	//Diffusion requires two boundary conditions whereas advection only requires one!
	//NOTE (bruno): cell cirulators can be use to get all cells having bcid2 has a vertex more efficiently (see e.g. FlowBoundingSphere.ipp:721)
	
	
	  
    	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{
		for (unsigned int ngb=0;ngb<4;ngb++){ 
			if (cell->vertex(ngb)->info().id() == bcid1){
			  cell->info().solute() = bcconcentration1;
			}
			if (D != 0){if (cell->vertex(ngb)->info().id() == bcid2){
			cell->info().solute() = bcconcentration2;			  
			}}
			  
		}
	}
	
}

double SoluteFlowEngine::getConcentrationPlane (double Yobs,double Yr, int xyz)
{
	//Get the concentration within a certain plane (Y_obs), whilst the cells are located in a small volume around this plane
	//The concentration in cells are weighed for their distance to the observation point
	//for a point on the x-axis (xyz=0), point on the y-axis (xyz=1), point on the z-axis (xyz=2)
	double sumConcentration = 0.0;
	double sumFraction=0.0;
	double concentration=0.0;
	//Find cells within designated volume

	 FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{
	CGT::Point& p1 = cell->info();
	if (abs(p1[xyz]) < abs(abs(Yobs) + abs(Yr))){
	  if(abs(p1[xyz]) > abs(abs(Yobs) - abs(Yr))){
		sumConcentration += cell->info().solute()*(1-(abs(p1[xyz])-abs(Yobs))/abs(Yr));
		sumFraction += (1-(abs(p1[xyz])-abs(Yobs))/abs(Yr));	
	}
	}
	}
      concentration = sumConcentration / sumFraction;
      return concentration;
}

double SoluteFlowEngine::checkMassBalance()
{
  double percentage = 0.0;
  double temp = 0.0;
  temp = massIntegrationStorage() - massinBC + totalmassBCout;
  percentage = 100.00*(temp - totalmassBCin)/temp;
  return percentage;
}


void SoluteFlowEngine::massIntegrationBC()
{    	
  bool alreadydone = false; 
  double masstotalin = 0.0;
  double masstotalout = 0.0;
  double flux=0.0;
  double fluxout = 0.0;
  double qin = 0.0;
  double q = 0.0;
  massinBC = 0.0;
  FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{
	for (unsigned int vertexngb=0;vertexngb<4;vertexngb++){ 
	if (cell->vertex(vertexngb)->info().id() == bcid1 && alreadydone == false){
			 massinBC += bcconcentration1*(abs(cell->info().volume()) - abs(solver->volumeSolidPore(cell)));
			 for (unsigned int ngb=0;ngb<4;ngb++){
			 q = (cell->info().kNorm() [ngb])* ( cell->neighbor ( ngb )->info().p()-cell->info().p());
			 qin=abs(min(0.0,q));
			 flux += qin*cell->info().solute()*scene->dt;
			 }
			 alreadydone = true;
			
	}
	if (cell->vertex(vertexngb)->info().id() == bcid2 && alreadydone == false){
			 for (unsigned int ngb=0;ngb<4;ngb++){
			 q = (cell->info().kNorm() [ngb])* ( cell->neighbor ( ngb )->info().p()-cell->info().p());
			 qin=abs(max(0.0,q));
			 fluxout += qin*cell->neighbor (ngb)->info().solute()*scene->dt;
			 }
			 alreadydone = true;
	
	}
	}
	alreadydone=false;
	masstotalin +=flux;
	masstotalout += fluxout;
	flux = 0.0;
	fluxout = 0.0;
	}
	totalmassBCin += masstotalin;
	totalmassBCout += masstotalout;
	//cerr << endl << "in "<< totalmassBCin <<"out "<<totalmassBCout << "total "<<massIntegrationStorage();
	
}	

double SoluteFlowEngine::massIntegrationStorage()
{
  double mass = 0.0;
  FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{
	 mass += cell->info().solute()*(abs(cell->info().volume()) - abs(solver->volumeSolidPore(cell)));
	}
  return mass;
}



YADE_PLUGIN ( ( SoluteFlowEngine ) );

#endif //SOLUTE_FLOW
#endif //FLOW_ENGINE

#endif /* YADE_CGAL */