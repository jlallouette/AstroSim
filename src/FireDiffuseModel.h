/*----------------------------------------------------------------------------
  AstroSim: Simulation of astrocyte networks Ca2+ dynamics
  Copyright (c) 2016-2017 Jules Lallouette
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/

#ifndef IFMODEL_H
#define IFMODEL_H

#include <vector>
#include "Model.h"
#include "ODEProblems.h"
#include "ParamHandler.h"
#include "CouplingFunction.h"
#include "StimulationStrat.h"
#include "Network.h"
#include "MetricComputeStrat.h"
#include "ChIModelMetrics.h"
#include "StimulationMetrics.h"

#define FIREDIFFUSEMODEL_NBVALS_PER_CELL 1

namespace AstroModel
{

/**********************************************************************/
/* Fire Diffuse Cell                                                  */
/**********************************************************************/
	class FireDiffuseCell : public SaveAndLoadFromStream
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::FireDiffuseNetFunct;
		friend class ODE::FireDiffuseCellFunct;
		friend class FireDiffuseModel;

	public:
		static std::string ClassName;
		static unsigned int NbValsPerCell;
		//===========================================================||
		// Enumeration binding indices to values names               ||
		//===========================================================||
		enum DynValNames {
			C = 0// Messenger concentration
		};

		//===========================================================||
		// Static constant equilibrium values for the variable       ||
		//===========================================================||
		static double DefaultC;

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		FireDiffuseCell(const ODENetworkDynamicsModel<CouplingFunction, 
			FireDiffuseCell> * _model, double *_dv = 0, bool _fv = false);
		// Copy constructor
		FireDiffuseCell(const FireDiffuseCell & c);
		// Destructor
		virtual ~FireDiffuseCell();

		//===========================================================||
		// Model type methods                                        ||
		//===========================================================||
		// Initializes the cell to default values
		virtual void Initialize();
		// Set the cell to equilibrium
		virtual void SetToEquilibrium();
		// Update value names in ODEProblem
		virtual void SetValNamesPostfix(std::string pf);
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the cell from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the cell to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		// Change dynamic values to given pointer
		virtual void SetDynVals(double *_dv, bool _f);
		// Return the number of dyn vals per cell
		virtual unsigned int GetNbDynVals();
		// Return the value of a dyn val
		inline double GetDynVal(int name) const
			{ return dynVals[name]; }
		// Return the number of dyn vals
		inline unsigned int GetNbDynVals() const
			{ return nbDynVals;}

	protected:
		static double DefaultThresh;
		static double DefaultActQuant;
		static double DefaultTauRefr;
		static double DefaultDegrad;
		static double DefaultActivDelay;

		//===========================================================||
		// Links to other objects                                    ||
		//===========================================================||
		ODE::FireDiffuseCellFunct * funct; // Derivative computation function
		const FireDiffuseModel * model;     // Pointer to parent model

		//===========================================================||
		// Cell biochemical parameters                               ||
		//===========================================================||
		double thresh; // Firing thresh
		double actQuant; // generated quantity on firing
		double tauRefr; // Time to spend in refractory period
		double degrad; // Degradation rate of messenger
		double activDelay; // Delay between thresh crossing and activation

		double refractTime;
		bool willActivate;
		double actTime;

		//===========================================================||
		// Dynamic values                                            ||
		//===========================================================||
		double *dynVals; // Dynamic values
		bool freeDynVals;
		unsigned int nbDynVals;
		double totFlux;  // Total fluxes of messenger from the cell
	};

/**********************************************************************/
/* Fire Diffuse Model                                                 */
/**********************************************************************/
	class FireDiffuseModel : 
		public ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>, 
		public StimulableCellNetwork
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::FireDiffuseNetFunct;

	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		FireDiffuseModel(ParamHandler & h = ParamHandler::GlobalParams);
		// Constructor from stream
		FireDiffuseModel(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams);
		// Destructor
		virtual ~FireDiffuseModel();

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		// Sets up cells and allocate data for ODEs
		virtual void SetUpCellsAndODEs(unsigned int nbCells);
		// Initializes the model
		virtual void Initialize(ResultSaver saver = ResultSaver::NullSaver);
		// Method called before solving the ODE problem
		virtual bool PreSimulationCall(ResultSaver saver);
		// Method called after solving the ODE problem
		virtual bool PostSimulationCall(ResultSaver saver);
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Load funct from stream
		virtual bool LoadFunctFromStream(std::ifstream & stream);
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model computations                                ||
		//===========================================================||
		// Computes fluxes across cells
		virtual void ComputeFluxes(double t);

		//===========================================================||
		// Setters and callbacks                                     ||
		//===========================================================||
		// Changes the flux of a cell
		void ModifFluxes(unsigned int i, double flux);
		// Notifies the model that all values have been updated for timestep t
		virtual void UpdateVals(double t);

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Add a network metric to the network
		virtual bool AddMetric(Metric *_m, bool _f = false);

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const;
		// Returns the total fluxes going out of the cell
		virtual double GetTotalFlux(unsigned int i) const;
		// Return the dynamic value that constitutes the excitable part of the system
		virtual double GetExcDynVal(unsigned int cellNb) const;

		// Returns a const ref on the network
		virtual const AbstractNetwork & GetNetwork() const 
			{return NetworkDynamicsModel<CouplingFunction>::GetNetwork();}
		// Returns tStart
		virtual const double & GetTStart() const 
			{ return ODE::ODEProblem<double, double>::GetTStart(); }
		// Return tEnd
		virtual const double & GetTEnd() const 
			{ return ODE::ODEProblem<double, double>::GetTEnd(); }
		// Returns the number of cells in the model
		inline unsigned int GetNbCells() const 
			{ return ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::GetNbCells(); }
		// Set all cells to equilibrium
		virtual void SetAllCellsToEquilibrium()
			{ return ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::SetAllCellsToEquilibrium(); }
		// Returns the neighbors of cell i
		virtual const std::vector<unsigned int> & GetNeighbors(unsigned int i) const
			{ return ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::GetNeighbors(i); }
		// Is the given cell currently stimulated ?
		virtual bool IsStimulated(unsigned int ind) const
			{ return Stimulable::IsStimulated(ind); }

	protected:

		virtual double getModelVersionNum() const { return 0.1; }

		//===========================================================||
		// Network parameters                                        ||
		//===========================================================||

		//===========================================================||
		// Associated ODE Problem                                    ||
		//===========================================================||
		ODE::ODESolver<double, double> * solver;

		//===========================================================||
		// Metrics                                                   ||
		//===========================================================||
		SortedMetrics<FireDiffuseModel> metrics;
	};
}

#endif
