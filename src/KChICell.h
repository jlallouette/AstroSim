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

#ifndef KCHICELL_H
#define KCHICELL_H

#include "ChICell.h"

#define KCHIMODEL_NBVALS_PER_CELL 9

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
/* KChI Cell                                                          */
/**********************************************************************/
	class KChICell : public ChICell
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::KChICellFunct;
		friend class ODE::KChINetworkFunct;
		friend class KChIModel;

	public:
		static std::string ClassName;
		//===========================================================||
		// Enumeration binding indices to values names               ||
		//===========================================================||
		enum KDynValNames {
			Ko = 3, // K+ concentration in extracellular space
			Ki,     // K+ concentration in intracellular space
			Vm,     // Membrane potential
			Cer,    // Ca2+ Concentration in Endoplasmic Reticulum
			Gp,     // GJC enhancing phosphorylation ratio
			Gm      // GJC decreasing phosphorylation ratio
		};

		//===========================================================||
		// Static constant equilibrium values for the 6 variables    ||
		//===========================================================||
		static double DefaultKo;
		static double DefaultKi;
		static double DefaultVm;
		static double DefaultCer;
		static double DefaultGp;
		static double DefaultGm;

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		KChICell(const KChIModel * _model, double *_dv = 0, bool _fv = false);
		// Copy constructor
		KChICell(const KChICell & c);
		// Destructor
		virtual ~KChICell();

		//===========================================================||
		// Model type methods                                        ||
		//===========================================================||
		// Initializes the cell to default values
		virtual void Initialize();
		// Compute other parameters
		virtual void ComputeOtherParameters();
		// Set the cell to equilibrium
		virtual void SetToEquilibrium();
		// Update value names in ODEProblem
		virtual void SetValNamesPostfix(std::string pf);

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

	protected:
		static double DefaultVa;     
		static double DefaultSa;     
		static double DefaultalphExt;
		static double DefaultKoBl; 
		static double DefaultCaOut;
		static double DefaultClOut;
		static double DefaultClIn; 
		static double DefaultNaOut;
		static double DefaultNaIn; 
		static double DefaultCap;
		static double DefaultGClLeak;
		static double DefaultGNaLeak;
		static double DefaultOmegaK; 
		static double DefaultGKLeak; 
		static double DefaultGKirMax;
		static double DefaultKNa;
		static double DefaultKK;
		static double DefaultJNaKATPaseMax;
		static double DefaultGCaLeak; 
		static double DefaultPCaLType;
		static double DefaultOMP;     
		static double Defaultkpact;
		static double Defaultkpinh;
		static double DefaultkPhos;
		static bool   DefaultEnsureCaEq;

		//===========================================================||
		// Links to other objects                                    ||
		//===========================================================||

		//===========================================================||
		// Cell biochemical parameters                               ||
		//===========================================================||
		// Geometric parameters
		double Va;      // Cell volume (V_a)
		double Sa;      // Cell surface (S_a)
		double alphExt; // Ratio of extracellular space volume over cell volume (\alpha_{ext})

		// Ion concentrations
		double KoBl;  // K+ baseline concentration ([K^+_o]_{bl})
		double CaOut; // Ca2+ concentration in the extracellular space ([Ca^{2+}_o])
		double ClOut; // Cl- concentration in the extracellular space ([Cl^-_o]) 
		double ClIn;  // Cl- concentration in the intracellular space ([Cl^-_i])
		double NaOut; // Na+ concentration in the extracellular space ([Cl^-_o])
		double NaIn;  // Na+ concentration in the intracellular space ([Cl^-_i])

		// Electrical parameters
		double Cap;     // Cell capacitance (C_a)
		double GClLeak; // Conductance for Cl- ions
		double GNaLeak; // Conductance for Na+ ions

		// [K+] intake
		double OmegaK;        // K+ Clearance rate (\Omega_K)
		double GKLeak;        // K+ leak conductance (G^K_{leak}) [S.m^-2]
		double GKirMax;       // K+ maximum conductance for Kir channels (G^{max}_{Kir}) [S.m^-2]
		double KNa;           // Threshold Na+ value for the Na/K ATPase
		double KK;            // Threshold Na+ value for the Na/K ATPase
		double JNaKATPaseMax; // Maximum Na/K ATPase flux (mol.s^-1.m^-2)

		// [Ca2+] fluxes with extracellular space
		double GCaLeak;   // Leak Ca2+ conductance (G^{Ca}_{leak}) [S.m^-2]
		double PCaLType;  // Ca2+ permeability of L-type channels (P_{CaL})
		double OMP;       // PMCA pumping rate (O_{MP})

		// GJC Phosphorylation parameters
		double kpact; // Prop factor between OCK and O3K;
		double kpinh; // Prop factor between OPK and \dzeta = kP / kR
		double kPhos; // Rate of unphosphorylations by phosphatases

		// Dependencies between parameter
		bool ensureCaEq; // Ensures that Ca2+ leak is in equilibrium with PMCA outflux.

		//===========================================================||
		// Computed cell parameters (not saved)                      ||
		//===========================================================||
		double VolCyt;    // Volume of cytosol (V_{cyt})
		double VolER;     // Volume of Endoplasmic Reticulum (V_{ER})
		double VolExt;    // Volume of extracellular space (V_{ext})
		double LogForECl; // log for computing the Nernst potential for Chloride ions
		double LogForENa; // log for computing the Nernst potential for sodium ions
		double OCK;       // Phosphorylation rate of CamKII (O_{CK})
		double OPK;       // Phosphorylation rate of PKC (O_{PK})

		//===========================================================||
		// Dynamic values                                            ||
		//===========================================================||
		double KFluxOut; // Potassium flux in the extracellular space (going in) [mol.s-1]
		double KFluxIn; // Potassium flux in the intracellular space (curr Cell to other cell)
	};
}

#endif
