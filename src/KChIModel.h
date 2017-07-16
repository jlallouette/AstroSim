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

#ifndef KCHIMODEL_H
#define KCHIMODEL_H

#include "ChIModel.h"

namespace ODE
{
	// Forward declarations
	class KChICellFunct;
	class KChINetworkFunct;
}

namespace AstroModel
{

// Do not use this in simulation, model not correctly calibrated yet
/**********************************************************************/
/* KChI Model                                                         */
/**********************************************************************/
	class KChIModel : public ChIModel
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::KChICellFunct;
		friend class ODE::KChINetworkFunct;
		friend class KChICell;

	public:
		static std::string ClassName;
		enum GJCCompModel {
			SimpleEq = 0, // Considers that PKC can cancel CamKII phosphorylation
			DoubleEq      // Considers that only phosphatase can unphosphorylate
		};
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		KChIModel(ParamHandler & h = ParamHandler::GlobalParams);
		// Constructor from stream
		KChIModel(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams);
		// Destructor
		virtual ~KChIModel();

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Sets up cells and allocate data for ODEs
		virtual void SetUpCellsAndODEs(unsigned int nbCells);
		// Change ODE vals to given pointer (for external use)
		virtual void SetVals(double *_vals = 0, unsigned int _nbVals = 0);
		// Initializes the model
		virtual void Initialize(ResultSaver saver = ResultSaver::NullSaver);
		// Compute other parameters
		virtual void ComputeOtherParameters();
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Load cells from string and sets function
		virtual bool LoadFunctAndCellsFromStream(std::ifstream & stream);
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
		// Changes the extracellular potassium flux of a cell
		void ModifKoutFluxes(unsigned int i, double flux);

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		// Gives the total number of desired dyn vals (not equivalent to GetNbVals
		// from OPEProblem<double, double>
		virtual unsigned int GetTotNbDynVals() const;

	protected:
		static double DefaultSij;
		static double DefaultT;
		static GJCCompModel DefaultGJCComp;
		static double DefaultKDiffVoltThr; 
		static double DefaultVKirH;
		static double DefaultVKirS;
		static double DefaultVLTmHalf;  
		static double DefaultVLTmSlope; 
		static double DefaultVLThHalf;  
		static double DefaultVLThSlope; 
		static double DefaultKMP;       
		static double DefaultalphaP;       
		static double DefaultalphaM;       
		static double DefaultIP3BasalPerm; 
		static double DefaultKBasalPerm;   

		virtual double getModelVersionNum() const { return 4.2; }

		//===========================================================||
		// Network parameters                                        ||
		//===========================================================||
		// Geometric parameters
		double Sij; // Default surface of contact between astrocytes;

		// Other parameters
		double T; // Temparature in Kelvin
		GJCCompModel GJCComp; // Model used for computing GJC phosphorylations
		double KDiffVoltThr; // Threshold on intercellular voltage under which K+ does not diffuse

		// [K+] intake
		double VKirH;   // Half activating potential for Kir channels (V_h)
		double VKirS;   // Activating potential slope for Kir channels (V_s)

		// [Ca2+] fluxes with extracellular space
		double VLTmHalf;  // L-type channel half activating potential (V_{m,1/2}))
		double VLTmSlope; // L-type channel activating slope (V_{m,slope}))
		double VLThHalf;  // L-type channel half inactivating potential (V_{h,1/2}))
		double VLThSlope; // L-type channel inactivating slope (V_{h,slope}))
		double KMP;       // Ca2+ affinity of PMCA pumps (K_{MP})

		// GJC Phosphorylation parameters
		double alphaP;       // Ratio of increased GJC conductance over basal one (\alpha^+)
		double alphaM;       // Ratio of decreased GJC conductance over basal one (\alpha^-)
		double IP3BasalPerm; // Basal Permeability for intercellular IP3 diffusion
		double KBasalPerm;   // Basal Permeability for intercellular K+ diffusion

		//===========================================================||
		// Computed network parameters (not saved)                   ||
		//===========================================================||
		double FoRT; // F / (R * T)

		//===========================================================||
		// Metrics                                                   ||
		//===========================================================||
	};
}

#endif
