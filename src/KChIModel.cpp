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

#include "KChIModel.h"
#include "KChICell.h"

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_exp.h>

using namespace AstroModel;

// Do not use this in simulation, the model is not correctly calibrated yet
//********************************************************************//
//******************** K C H I  M O D E L ****************************//
//********************************************************************//

double KChIModel::DefaultSij           = 4.5239e-10;
double KChIModel::DefaultT             = 310.15;
KChIModel::GJCCompModel KChIModel::DefaultGJCComp = SimpleEq;
double KChIModel::DefaultKDiffVoltThr  = 1.0e-06;
double KChIModel::DefaultVKirH         = -0.082;
double KChIModel::DefaultVKirS         = 1.0;  
double KChIModel::DefaultVLTmHalf      = -0.050; // Should be -0.04 (Evans2013)
double KChIModel::DefaultVLTmSlope     = -0.005;
double KChIModel::DefaultVLThHalf      = -0.047; // Should be -0.037 (Evans 2013)
double KChIModel::DefaultVLThSlope     = 0.005;
double KChIModel::DefaultKMP           = 5.0e-05; 
double KChIModel::DefaultalphaP        = 50; // Arbitraty
double KChIModel::DefaultalphaM        = 0.5; // Arbitrary
double KChIModel::DefaultIP3BasalPerm  = 0;//0.2 * 1.2209e-04; // F*Cste with F GJC strength in ChI
double KChIModel::DefaultKBasalPerm    = 0.5e-06;//3.5e-06;//0.7e-06;//0.2e-6//1000

std::string KChIModel::ClassName = "KChIModel";

//**********************************************************************
// Default Constructor
//**********************************************************************
KChIModel::KChIModel(ParamHandler & h) : 
	ChIModel(h), Sij(DefaultSij), T(DefaultT), GJCComp(DefaultGJCComp),
	KDiffVoltThr(DefaultKDiffVoltThr), VKirH(DefaultVKirH), VKirS(DefaultVKirS),
	VLTmHalf(DefaultVLTmHalf), VLTmSlope(DefaultVLTmSlope), 
	VLThHalf(DefaultVLThHalf), VLThSlope(DefaultVLThSlope), KMP(DefaultKMP), 
	alphaP(DefaultalphaP), alphaM(DefaultalphaM), IP3BasalPerm(DefaultIP3BasalPerm), 
	KBasalPerm(DefaultKBasalPerm)
{
	TRACE("*** Initializing KChI Model ***")
	SetUpCellsAndODEs(cells.size());
}

//**********************************************************************
// Loading Constructor
//**********************************************************************
KChIModel::KChIModel(std::ifstream & stream, ParamHandler & h) : 
	ChIModel(stream, h)
{

}

//**********************************************************************
// Destructor
//**********************************************************************
KChIModel::~KChIModel()
{
}

//**********************************************************************
// Sets up cells and allocate data for ODEs
//**********************************************************************
void KChIModel::SetUpCellsAndODEs(unsigned int nbCells)
{
	SetFunct(new ODE::KChINetworkFunct(*this), true);
	AllocateMemory(nbCells * KCHIMODEL_NBVALS_PER_CELL);

	for (unsigned int i = 0 ; i < cells.size() ; ++i)
		delete cells[i];
	cells.clear();
	for (unsigned int i = 0 ; i < nbCells ; ++i)
	{
		cells.push_back(new KChICell(this, this->vals + 
			i * KCHIMODEL_NBVALS_PER_CELL, false));
	}

	// Calls SetVals with 0 / 0 arguments (default args)
	// So that it doesn't change the allocated vals.
	SetVals();
}

//**********************************************************************
// Change ODE vals to given pointer (for external use)
//**********************************************************************
void KChIModel::SetVals(double *_vals, unsigned int _nbVals)
{
	ODE::ODEProblem<double, double>::SetVals(_vals, _nbVals);

	for (unsigned int i = 0 ; i < cells.size() ; ++i)
	{
		cells[i]->SetDynVals(this->vals + i * KCHIMODEL_NBVALS_PER_CELL, false);
		cells[i]->SetValNamesPostfix(cells[i]->GetClassName() + "_" +
			StringifyFixed(i));
	}

	ODE::ODEProblem<double, double>::UseCurrentValsAsInitVals();
}

//**********************************************************************
// Initializes the model
//**********************************************************************
void KChIModel::Initialize(ResultSaver saver)
{
	ComputeOtherParameters();
	ChIModel::Initialize(saver);
}

//**********************************************************************
// Compute other parameters
//**********************************************************************
void KChIModel::ComputeOtherParameters()
{
	FoRT = GSL_CONST_MKSA_FARADAY / (GSL_CONST_MKSA_MOLAR_GAS * T);
}

//**********************************************************************
// Computes fluxes across cells
//**********************************************************************
void KChIModel::ComputeFluxes(double )
{
	static double perm;
	static double v;
	KChICell *kcell;

	for (unsigned int i = 0 ;  i < network->size() ; ++i)
	{
		kcell = dynamic_cast<KChICell *>(cells[i]);
		assert(kcell);
		cells[i]->totFlux = 0;
		cells[i]->caSpontLeak = false;
		kcell->KFluxIn = 0;
		kcell->KFluxOut = 0;
		for (unsigned int j = 0 ; j < network->GetNeighbors(i).size() ; ++j)
		{
			// Compute permeability between i and j
			switch (GJCComp) {
				case KChIModel::SimpleEq:
					perm = /*IP3BasalPerm * */ std::min(cells[i]->dynVals[KChICell::Gp], 
						cells[network->GetNeighbors(i)[j]]->dynVals[KChICell::Gp]) * (alphaP - alphaM) + alphaM;
					break;
				case KChIModel::DoubleEq:
					perm = /*IP3BasalPerm * */ (std::min(
						cells[i]->dynVals[KChICell::Gp] * (alphaP - 1.0) + 
						cells[i]->dynVals[KChICell::Gm] * (alphaM - 1.0), 
						cells[network->GetNeighbors(i)[j]]->dynVals[KChICell::Gp] * (alphaP - 1.0) + 
						cells[network->GetNeighbors(i)[j]]->dynVals[KChICell::Gm] * (alphaM - 1.0)) + 1.0);
					break;
				default:
					perm = 1.0;
			}
			cells[i]->totFlux += Sij / kcell->VolCyt * IP3BasalPerm * perm * 
				(*((*network)[i][network->GetNeighbors(i)[j]]))(cells[i]->dynVals[ChICell::IP3] - 
				cells[network->GetNeighbors(i)[j]]->dynVals[ChICell::IP3]);
			
			if (fabs(cells[i]->dynVals[KChICell::Vm] - 
				cells[network->GetNeighbors(i)[j]]->dynVals[KChICell::Vm]) > KDiffVoltThr)
			{
				v = FoRT * (cells[i]->dynVals[KChICell::Vm] - 
					cells[network->GetNeighbors(i)[j]]->dynVals[KChICell::Vm]);
				kcell->KFluxIn += Sij * KBasalPerm * perm * v * 
					(cells[i]->dynVals[KChICell::Ki] - 
					 cells[network->GetNeighbors(i)[j]]->dynVals[KChICell::Ki] * gsl_sf_exp(-v)
					) / (1.0 - gsl_sf_exp(-v));
			}
			else
			{
				kcell->KFluxIn += Sij * KBasalPerm * perm * (cells[i]->dynVals[KChICell::Ki] - 
						cells[network->GetNeighbors(i)[j]]->dynVals[KChICell::Ki]);
			}
		}
	}
}

//**********************************************************************
// Changes the extracellular potassium flux of a cell
//**********************************************************************
void KChIModel::ModifKoutFluxes(unsigned int i, double flux)
{
	KChICell *tmp = dynamic_cast<KChICell*>(cells[i]);
	tmp->KFluxOut += flux;
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool KChIModel::LoadFromStream(std::ifstream & stream)
{
	stream >> Sij;
	stream >> T;
	int tmpVal;
	stream >> tmpVal;
	GJCComp = (GJCCompModel) tmpVal;
	stream >> KDiffVoltThr;
	stream >> VKirH;
	stream >> VKirS;
	stream >> VLTmHalf;
	stream >> VLTmSlope;
	stream >> VLThHalf;
	stream >> VLThSlope;
	stream >> KMP;
	stream >> alphaP;
	stream >> alphaM;
	stream >> IP3BasalPerm;
	stream >> KBasalPerm;

	ComputeOtherParameters();

	bool ok = ChIModel::LoadFromStream(stream);
	return ok and (stream.good() or stream.eof());
}

//**********************************************************************
// Load cells from string and sets function
//**********************************************************************
bool KChIModel::LoadFunctAndCellsFromStream(std::ifstream & stream)
{
	bool ok = true;

	// ODE Function
	SetFunct(new ODE::KChINetworkFunct(*this), true);

	// Cells
	for (unsigned int i = 0 ; i < cells.size() ; ++i)
		delete cells[i];
	cells.clear();

	unsigned int nbCells;
	stream >> nbCells;
	//
	AllocateMemory(nbCells * KCHIMODEL_NBVALS_PER_CELL);
	for (unsigned int i = 0 ; i < nbCells ; ++i)
	{
		cells.push_back(new KChICell(this, this->vals + i * KCHIMODEL_NBVALS_PER_CELL, false));
		ok &= cells.back()->LoadFromStream(stream);
	}

	return ok;
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool KChIModel::SaveToStream(std::ofstream & stream) const
{
	bool ok = ChIModel::SaveToStream(stream);

	stream 
		<< Sij          << std::endl
		<< T            << std::endl
		<< GJCComp      << std::endl
		<< KDiffVoltThr << std::endl
		<< VKirH        << std::endl
		<< VKirS        << std::endl
		<< VLTmHalf     << std::endl
		<< VLTmSlope    << std::endl
		<< VLThHalf     << std::endl
		<< VLThSlope    << std::endl
		<< KMP          << std::endl
		<< alphaP       << std::endl
		<< alphaM       << std::endl
		<< IP3BasalPerm << std::endl
		<< KBasalPerm   << std::endl;

	return ok;
}

//**********************************************************************
// Returns a ParamHandler object with references to internal parameters
//**********************************************************************
ParamHandler KChIModel::BuildModelParamHandler()
{
	ParamHandler params;

	params <= "Sij"         , Sij;
	params <= "T"           , T;
	params <= "GJCComp"     , ((int)GJCComp);
	params <= "KDiffVoltThr", KDiffVoltThr;
	params <= "VKirH"       , VKirH;
	params <= "VKirS"       , VKirS;
	params <= "VLTmHalf"    , VLTmHalf;
	params <= "VLTmSlope"   , VLTmSlope;
	params <= "VLThHalf"    , VLThHalf;
	params <= "VLThSlope"   , VLThSlope;
	params <= "KMP"         , KMP;
	params <= "alphaP"      , alphaP;
	params <= "alphaM"      , alphaM;
	params <= "IP3BasalPerm", IP3BasalPerm;
	params <= "KBasalPerm"  , KBasalPerm;

	params += ChIModel::BuildModelParamHandler();
	return params;
}

//********************************************************************//
// Gives the total number of desired dyn vals (not equivalent to 
// GetNbVals from OPEProblem<double, double>
//********************************************************************//
unsigned int KChIModel::GetTotNbDynVals() const
{
	return cells.size() * KCHIMODEL_NBVALS_PER_CELL;
}

