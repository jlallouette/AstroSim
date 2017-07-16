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

#include "FireDiffuseModel.h"
#include "ODESolvers.h"

#include "ErrorCodes.h"

using namespace AstroModel;
using namespace std;

//********************************************************************//
//*********** F I R E   D I F F U S E   C E L L **********************//
//********************************************************************//

double FireDiffuseCell::DefaultC = 0.0;
double FireDiffuseCell::DefaultThresh = 0.0004;
double FireDiffuseCell::DefaultActQuant = 0.0008;
double FireDiffuseCell::DefaultTauRefr = 7.5;
double FireDiffuseCell::DefaultDegrad = 0.3;
double FireDiffuseCell::DefaultActivDelay = 5.0;


std::string FireDiffuseCell::ClassName = "FireDiffuseCell";
unsigned int FireDiffuseCell::NbValsPerCell = FIREDIFFUSEMODEL_NBVALS_PER_CELL;

//**********************************************************************
// Default Constructor
//**********************************************************************
FireDiffuseCell::FireDiffuseCell(const ODENetworkDynamicsModel<CouplingFunction, 
	FireDiffuseCell> * _model, double *_dv, bool _fv) :
	model(0),
	thresh(DefaultThresh), actQuant(DefaultActQuant), tauRefr(DefaultTauRefr),
	degrad(DefaultDegrad), activDelay(DefaultActivDelay), 
	refractTime(-99999), willActivate(false), actTime(0),
	dynVals(_dv), freeDynVals(_fv), nbDynVals(FIREDIFFUSEMODEL_NBVALS_PER_CELL), 
	totFlux(0)
{
	model = static_cast<const FireDiffuseModel*>(_model);
	funct = new ODE::FireDiffuseCellFunct(*this);
	if (not _dv)
	{
		dynVals = new double[nbDynVals];
		freeDynVals = true;
	}
	FireDiffuseCell::Initialize();
}

//**********************************************************************
// Copy constructor
//**********************************************************************
FireDiffuseCell::FireDiffuseCell(const FireDiffuseCell & c) : model(c.model), 
	thresh(c.thresh), actQuant(c.actQuant), tauRefr(c.tauRefr), degrad(c.degrad), 
	activDelay(DefaultActivDelay),
	refractTime(-99999), willActivate(false), actTime(0),
	dynVals(c.dynVals), freeDynVals(c.freeDynVals), 
	nbDynVals(c.nbDynVals), totFlux(c.totFlux)
{
	funct = new ODE::FireDiffuseCellFunct(*this);
	if (freeDynVals)
	{
		dynVals = new double[nbDynVals];
		for (unsigned int i = 0 ; i < nbDynVals ; ++i)
			dynVals[i] = c.dynVals[i];
	}
	Initialize();
}

//**********************************************************************
// Destructor
//**********************************************************************
FireDiffuseCell::~FireDiffuseCell()
{
	delete funct;
	if (freeDynVals)
		delete[] dynVals;
}

//**********************************************************************
//**********************************************************************
void FireDiffuseCell::Initialize()
{
	dynVals[0] = DefaultC,
	totFlux = 0;
	refractTime = -999999;
	willActivate = false;
	actTime = 0;

	// Parameters
	thresh = DefaultThresh;
	actQuant = DefaultActQuant;
	tauRefr = DefaultTauRefr;
	degrad = DefaultDegrad;
	activDelay = DefaultActivDelay;
}

//**********************************************************************
//**********************************************************************
void FireDiffuseCell::SetToEquilibrium()
{
	dynVals[C] = DefaultC;
	totFlux = 0;
	refractTime = -999999;
	willActivate = false;
	actTime = 0;
}

//**********************************************************************
//**********************************************************************
bool FireDiffuseCell::LoadFromStream(std::ifstream & stream)
{
	stream >> thresh;
	stream >> actQuant;
	stream >> tauRefr;
	stream >> degrad;
	stream >> activDelay;
	return not stream.eof();
}

//**********************************************************************
//**********************************************************************
bool FireDiffuseCell::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< thresh << endl
		<< actQuant << endl
		<< tauRefr << endl
		<< degrad << endl
		<< activDelay << endl;
	return stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler FireDiffuseCell::BuildModelParamHandler()
{
	ParamHandler params;
	params <= "DefaultCVal", DefaultC;
	params <= "thresh", DefaultThresh;   
	params <= "actQuant", DefaultActQuant;   
	params <= "tauRefr", DefaultTauRefr;   
	params <= "degrad", DefaultDegrad;   
	params <= "activDelay", DefaultActivDelay;
	return params;
}

//**********************************************************************
// Change dynamic values to given pointer
//**********************************************************************
void FireDiffuseCell::SetDynVals(double *_dv, bool _f)
{
	if (freeDynVals)
		delete[] dynVals;
	dynVals = _dv;
	freeDynVals = _f;
}

//**********************************************************************
// Return the number of dyn vals per cell
//**********************************************************************
unsigned int FireDiffuseCell::GetNbDynVals()
{
	return nbDynVals;
}

//**********************************************************************
// Update value names in ODEProblem
//**********************************************************************
void FireDiffuseCell::SetValNamesPostfix(std::string pf)
{
	model->AddPostfixToValName(dynVals + C, pf + "_C");
}

//********************************************************************//
//************* F I R E   D I F F U S E   M O D E L ******************//
//********************************************************************//

string FireDiffuseModel::ClassName = "FireDiffuseModel";

//**********************************************************************
// Default Constructor
//**********************************************************************
FireDiffuseModel::FireDiffuseModel(ParamHandler & h) : 
	ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::ODENetworkDynamicsModel(h),
	StimulableCellNetwork::StimulableCellNetwork(h)
{
	TRACE("*** Initializing Fire Diffuse Model ***")
	SetFunct(new ODE::FireDiffuseNetFunct(*this), true);
}

//**********************************************************************
// Loading Constructor
//**********************************************************************
FireDiffuseModel::FireDiffuseModel(std::ifstream & stream, ParamHandler & h) : 
	ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::ODENetworkDynamicsModel(h),
	StimulableCellNetwork::StimulableCellNetwork(h)
{
	if (not LoadFromStream(stream))
		cerr << "Failed to load the model !" << endl;

}

//**********************************************************************
// Destructor
//**********************************************************************
FireDiffuseModel::~FireDiffuseModel()
{
}

//**********************************************************************
// Sets up cells and allocate data for ODEs
//**********************************************************************
void FireDiffuseModel::SetUpCellsAndODEs(unsigned int nbCells)
{
	SetFunct(new ODE::FireDiffuseNetFunct(*this), true);
	ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::SetUpCellsAndODEs(nbCells);
}

//**********************************************************************
// Initializes the model
//**********************************************************************
void FireDiffuseModel::Initialize(ResultSaver saver)
{
	ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::Initialize(saver);
	Stimulable::Initialize();

	// Initialize metrics
	metrics.InitializeMetricsDefault();
}

//**********************************************************************
// Method called before solving the ODE problem
//**********************************************************************
bool FireDiffuseModel::PreSimulationCall(ResultSaver saver)
{
	bool ok = ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::PreSimulationCall(saver);

	if (preRunToEqu)
	{
		// Temporarily deactivate metrics and stimStrats
		std::vector<StimulationStrat *> tmpStimStrats = stimStrats;
		stimStrats.clear();
		isPreRunning = true;

		solver->Solve(*this, tStart, tStart + preRunTime);

		stimStrats = tmpStimStrats;
		isPreRunning = false;
	}
	ODE::ODEProblem<double, double>::UseCurrentValsAsInitVals();

	// Save wanted data
	ok &= metrics.ComputeMetrics<NeedFrequentUpdateMetric>(*this);
	return ok;
}

//**********************************************************************
// Method called before solving the ODE problem
//**********************************************************************
bool FireDiffuseModel::PostSimulationCall(ResultSaver saver)
{
	bool ok = true;

TRACE_UP("*** Computing and saving after simulation metric ***")
	ok &= ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::PostSimulationCall(saver);
TRACE_DOWN("*** After simulation metric data saved ***")

	// Save computed dynamic data and metrics
TRACE_UP("*** Saving dynamic data ***")
	Stimulable::SaveStimStratMetrics(saver);
TRACE_DOWN("*** Finished saving dynamic data ***")

	return ok;
}

//**********************************************************************
// Computes fluxes across cells
//**********************************************************************
void FireDiffuseModel::ComputeFluxes(double )
{
	for (unsigned int i = 0 ;  i < network->size() ; ++i)
	{
		cells[i]->totFlux = 0;
		for (unsigned int j = 0 ; j < network->GetNeighbors(i).size() ; ++j)
			cells[i]->totFlux += (*((*network)[i][network->GetNeighbors(i)[j]]))(
					cells[i]->dynVals[FireDiffuseCell::C] - cells[network->GetNeighbors(i)[j]]->dynVals[FireDiffuseCell::C]);
	}
}

//**********************************************************************
// Changes the flux of a cell
//**********************************************************************
void FireDiffuseModel::ModifFluxes(unsigned int i, double flux)
{
	cells[i]->totFlux += flux;
}

//**********************************************************************
// Notifies the model that all values have been updated for timestep t
//**********************************************************************
void FireDiffuseModel::UpdateVals(double t)
{ 
	tCurr = t; 
	if (not isPreRunning)
	{
		ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::UpdateVals(t);

		for (unsigned int i = 0 ; i < cells.size() ; ++i)
		{
			if (not cells[i]->willActivate and (cells[i]->dynVals[0] >= cells[i]->thresh) and ((tCurr - cells[i]->refractTime) > cells[i]->tauRefr))
			{
				cells[i]->willActivate = true;
				cells[i]->actTime = tCurr;
			}
			if (cells[i]->willActivate and ((tCurr - cells[i]->actTime) > cells[i]->activDelay))
			{
				cells[i]->refractTime = tCurr;
				cells[i]->dynVals[0] = cells[i]->actQuant;
				cells[i]->willActivate = false;
			}
		}
	}
}

//**********************************************************************
// Dynamically dispatch a metric according to its type
//**********************************************************************
bool FireDiffuseModel::AddMetric(Metric *_m, bool _f)
{
	if (ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::AddMetric(_m, _f))
		return true;
	else if (Stimulable::AddMetric(_m, _f))
		return true;
	else
		return metrics.AddMetricAndDependencies(_m, _f, this);
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool FireDiffuseModel::LoadFromStream(std::ifstream & stream)
{
	bool ok = ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::LoadFromStream(stream);
	ok &= Stimulable::LoadFromStream(stream);
	if (ok)
	{
		// metrics
		metrics.FreeAndClean();
		ok &= metrics.LoadFromStream(stream);

		return ok and (stream.good() or stream.eof());
	}
	else 
		return false;
}

//**********************************************************************
// Load funct from stream
//**********************************************************************
bool FireDiffuseModel::LoadFunctFromStream(std::ifstream & )
{
	// ODE Function
	SetFunct(new ODE::FireDiffuseNetFunct(*this), true);
	return true;
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool FireDiffuseModel::SaveToStream(std::ofstream & stream) const
{
	bool ok = ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::SaveToStream(stream);
	ok &= Stimulable::SaveToStream(stream);

	// Metrics
	ok &= metrics.SaveToStream(stream);

	return ok;
}

//**********************************************************************
// Returns a ParamHandler object with references to internal parameters
//**********************************************************************
ParamHandler FireDiffuseModel::BuildModelParamHandler()
{
	ParamHandler params;
	params += ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::BuildModelParamHandler();
	params += Stimulable::BuildModelParamHandler();

	params += metrics.BuildModelParamHandler();
	return params;
}

//**********************************************************************
// Returns all metrics and submetrics
//**********************************************************************
vector<Metric *> FireDiffuseModel::GetAllMetrics() const
{
	vector<Metric *> mTot, mTemp;

	mTemp = ODENetworkDynamicsModel<CouplingFunction, FireDiffuseCell>::GetAllMetrics();
	for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
		mTot.push_back(mTemp[i]);

	mTemp = Stimulable::GetAllMetrics();
	for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
		mTot.push_back(mTemp[i]);

	for (unsigned int i = 0 ; i < metrics.GetMetricsRaw().size() ; ++i)
		mTot.push_back(metrics.GetMetricsRaw()[i]);

	return mTot;
}

//********************************************************************//
// Returns the total fluxes going out of the cell
//********************************************************************//
double FireDiffuseModel::GetTotalFlux(unsigned int i) const
{
	assert(i < cells.size());
	return cells[i]->totFlux;
}

//********************************************************************//
// Return the dynamic value that constitutes the excitable part of the system
//********************************************************************//
double FireDiffuseModel::GetExcDynVal(unsigned int cellNb) const
{
	return GetDynVal(cellNb, FireDiffuseCell::C);
}

