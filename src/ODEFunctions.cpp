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

#include "ODEFunctions.h"
#include "ChICell.h"
#include "KChICell.h"
#include "ChIModel.h"
#include "KChIModel.h"
#include "Synapse.h"
#include "Neuron.h"
#include "NeuronNetModels.h"
#include "AstroNeuroModel.h"
#include "FireDiffuseModel.h"

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include <math.h>

using namespace ODE;
using namespace AstroModel;
using namespace std;

std::string ChICellFunct::ClassName = "ChICellFunct";
std::string KChICellFunct::ClassName = "KChICellFunct";
std::string TMSynapseFunct::ClassName = "TMSynapseFunct";
std::string TMSynapseFunctOptim::ClassName = "TMSynapseFunctOptim";
std::string SFALIFNeuronFunct::ClassName = "SFALIFNeuronFunct";
std::string DummyNeuronFunct::ClassName = "DummyNeuronFunct";
std::string ChINetworkFunct::ClassName = "ChINetworkFunct";
std::string KChINetworkFunct::ClassName = "KChINetworkFunct";
std::string NeuronNetworkFunc::ClassName = "NeuronNetworkFunc";
std::string AstroNeuronNetFunc::ClassName = "AstroNeuronNetFunc";
std::string FireDiffuseCellFunct::ClassName = "FireDiffuseCellFunc";
std::string FireDiffuseNetFunct::ClassName = "FireDiffuseNetFunc";

//********************************************************************//
//********************** C H I  M O D E L ****************************//
//********************************************************************//


//**********************************************************************
// Constructor
//**********************************************************************
ChICellFunct::ChICellFunct(AstroModel::ChICell & _cell) : cell(_cell), model(*_cell.model)
{

}

//**********************************************************************
// Constructor
//**********************************************************************
ChINetworkFunct::ChINetworkFunct(AstroModel::ChIModel & _model) : model(_model)
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void ChICellFunct::CompFunc(const double &, const double *v, double *f) const
{
	static double Ca2;
	static double Ca4;
	static double Q2;
	static double hInf; 
	static double tauH; 
	static double mInf; 
	static double nInf; 
	static double chanProb;
	static double Jchan;
	static double Jleak;
	static double Jpump;
	static double Jspont;
	static double Pplcd;
	static double D5p;  
	static double D3k;  
	static double dIP3i;

	double & diffCa    = f[0];
	double & diffh     = f[1];
	double & diffIP3   = f[2];
	const double & Ca  = v[0];
	const double & h   = v[1];
	const double & IP3 = v[2];

	Q2     = model.d2 * (IP3 + model.d1) / (IP3 + model.d3);

	hInf   = Q2  / (Q2 + Ca);
	tauH   = 1.0 / (cell.a2 * (Q2 + Ca));
	mInf   = IP3  / (IP3 + model.d1);
	nInf   = Ca  / (Ca + model.d5);

	Ca2 = INTPOW2(Ca);
	Ca4 = INTPOW2(Ca2);
	chanProb = mInf * nInf * h;
	Jchan  = cell.rC  * (model.C0 - (1.0 + cell.c1) * Ca) * INTPOW3(chanProb);
	Jleak  = cell.rL  * (model.C0 - (1.0 + cell.c1) * Ca);
	Jpump  = cell.vER * Ca2 / (INTPOW2(cell.Ker) + Ca2);
	// Modified version to match Osama's formula
	Jspont = cell.caSpontLeak ? cell.rL : 0;

	Pplcd  = cell.vd * cell.kd / (cell.kd + IP3) * Ca2 / (Ca2 + INTPOW2(cell.Kplcd));
	D5p    = model.r5p * IP3;
	D3k    = model.v3k * Ca4 / (Ca4 + INTPOW4(model.K3k)) * IP3 / (IP3 + cell.k3);

	dIP3i  = Pplcd - D3k - D5p;

	diffCa  = Jchan + Jleak - Jpump + Jspont;
	diffh   = (hInf - h) / tauH;
	diffIP3 = dIP3i - cell.totFlux + cell.gluIP3Prod;
}

//**********************************************************************
// Function network operator, it returns the derivatives of every dynVal
//**********************************************************************
void ChINetworkFunct::CompFunc(const double & t, const double *v, double *f) const
{
	static unsigned int dynValDim = model.GetNbDynVal(0);

	model.ComputeFluxes(t);
	model.Stimulate(t);

	for (unsigned int i = 0 ; i < model.cells.size() ; ++i)
		model.cells[i]->funct->CompFunc(t, v + i * dynValDim, f + i * dynValDim);
}

//********************************************************************//
//******************** K C H I  M O D E L ****************************//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
KChICellFunct::KChICellFunct(AstroModel::KChICell & _cell) : cell(_cell), 
	model(*dynamic_cast<const KChIModel *>(_cell.model))
{
	
}

//**********************************************************************
// Constructor
//**********************************************************************
KChINetworkFunct::KChINetworkFunct(AstroModel::KChIModel & _model) : model(_model)
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void KChICellFunct::CompFunc(const double & t, const double *v, double *f) const
{
	static double F = GSL_CONST_MKSA_FARADAY;
	static double R = GSL_CONST_MKSA_MOLAR_GAS;
	static double T;
	static double RToF;

	static double Q2;
	static double hInf; 
	static double OmegaH;
	static double mInf; 
	static double Jchan;
	static double Jleak;
	static double Jpump;
	static double Jspont;
	static double Pplcd;
	static double D5p;  
	static double D3k;  

	static double ECa;
	static double LTv;
	static double expLTv;
	static double mLT;
	static double hLT;
	static double ICaLeak;
	static double ICaLType;
	static double ICaPMCA;
	static double ICa;

	static double Ek;
	static double GKir;
	static double Ik;
	static double JNaKATPase;

	static double ICl;
	static double INa;

	static double kCamAct;
	static double kPKCAct;

	const double & Ca  = v[0]; // Cell-averaged Ca2+ concentration
	const double & h   = v[1]; // Fraction of non inactivated IP3R channels on the ER membrane
	const double & IP3 = v[2]; // Cell-averaged concentration of IP3 second messenger
	const double & Ko  = v[3]; // K+ concentration in extracellular space
	const double & Ki  = v[4]; // K+ concentration in intracellular space
	const double & Vm  = v[5]; // Membrane potential
	const double & Cer = v[6]; // Ca2+ Concentration in Endoplasmic Reticulum
	const double & Gp  = v[7]; // GJC enhancing phosphorylation ratio
	const double & Gm  = v[8]; // GJC decreasing phosphorylation ratio

	T = model.T;
	RToF = R * T / F;

	// h
	Q2     = model.d2 * (IP3 + model.d1) / (IP3 + model.d3);
	hInf   = Q2 / (Q2 + Ca);
	OmegaH = (cell.a2 * (Q2 + Ca));
	mInf   = HILL1(IP3, model.d1) * HILL1(Ca, model.d5);

	// Ca
	Jchan    = cell.VolER * cell.rC  * pow(mInf * h, 3) * (Cer - Ca);
	Jleak    = cell.VolER * cell.rL  * (Cer - Ca);
	Jpump    = cell.VolCyt * cell.vER * HILL2(Ca, cell.Ker);
	Jspont   = cell.caSpontLeak ? cell.rL : 0.0;

	ECa = RToF / 2.0 * gsl_sf_log(cell.CaOut / Ca);
	LTv  = 2.0 * Vm / RToF;
	expLTv = gsl_sf_exp(-LTv);
	mLT  = 1.0 / (1.0 + gsl_sf_exp((Vm - model.VLTmHalf) / model.VLTmSlope));
	hLT  = 1.0 / (1.0 + gsl_sf_exp((Vm - model.VLThHalf) / model.VLThSlope));

	ICaLeak  = cell.GCaLeak * (Vm - ECa);
	ICaLType = mLT * hLT * 2.0 * cell.PCaLType * F * LTv * (Ca - cell.CaOut * expLTv) / (1.0 - expLTv);
	ICaPMCA  = 2.0 * F * cell.OMP * HILL2(Ca, model.KMP);
	ICa      = cell.Sa * (ICaLeak + ICaLType + ICaPMCA);

	// IP3
	Pplcd  = cell.vd * cell.kd / (cell.kd + IP3) * HILL2(Ca, cell.Kplcd);
	D5p    = model.r5p * IP3;
	D3k    = model.v3k * HILLn(Ca, model.K3k, 4) * HILL1(IP3, cell.k3);

	// K^+_o
	Ek   = RToF * gsl_sf_log(Ko / Ki);
	GKir = cell.GKirMax / sqrt(Ko * (1.0 + gsl_sf_exp((Vm - model.VKirH - Ek) / model.VKirS)));
	JNaKATPase = cell.Sa * cell.JNaKATPaseMax *
		HILLn(cell.NaIn, cell.KNa, 1.5) * HILL1(Ko, cell.KK);
	Ik   = cell.Sa * (cell.GKLeak + GKir) * (Vm - Ek);

	// Vm
	ICl = cell.Sa * cell.GClLeak * (Vm + RToF * cell.LogForECl);
	INa = cell.Sa * cell.GNaLeak * (Vm - RToF * cell.LogForENa);

	// Gp and Gm
	kCamAct = cell.OCK * HILLn(Ca, cell.k3, 4);
	kPKCAct = cell.OPK * HILL1(Ca, cell.kpi);

	cell.gluIP3Prod = 0;

	f[0] = (Jchan + Jleak - Jpump + Jspont - ICa / (2.0 * F)) / cell.VolCyt;
	f[1] = (hInf - h) * OmegaH;
	f[2] = Pplcd - D3k - D5p - cell.totFlux + cell.gluIP3Prod;

	if (t > 0) {
	TRACE("t = " << t-1.0 << " // Ca = " << v[0] << " // h = " << v[1] << " // I = " << v[2] << " // Ki = " << v[4] << " // Vm = " << v[5] << " // Cer = " << v[6] << " // dKi = " << cell.KFluxIn << " // ICa = " << ICa << " // Jk = " << -Ik/F + 2.0*JNaKATPase) }
	assert(not isnan(f[2]));
	// K+o
	f[3] = (Ik / F - 2.0 * JNaKATPase + cell.KFluxOut)/cell.VolExt - cell.OmegaK * (Ko - cell.KoBl);
	// K+i
	f[4] = (-Ik / F + 2.0 * JNaKATPase - cell.KFluxIn) / cell.VolCyt;
	// Vm
	f[5] = -(Ik + ICa + ICl + INa + F * (cell.KFluxIn + JNaKATPase)) / cell.Cap;
	// Cer
	f[6] = (-Jchan - Jleak + Jpump - Jspont) / cell.VolER;
	switch (model.GJCComp)
	{
		case KChIModel::SimpleEq:
			f[7] = kCamAct * (1.0 - Gp) - kPKCAct * Gp;
			f[8] = 0.0;
			break;
		case KChIModel::DoubleEq:
			f[7] = kCamAct * (1.0 - Gp - Gm) - cell.kPhos * Gp;
			f[8] = kPKCAct * (1.0 - Gp - Gm) - cell.kPhos * Gm;
			break;
		default:
			f[7] = 0.0;
			f[8] = 0.0;
	}
}

//**********************************************************************
// Function network operator, it returns the derivatives of every dynVal
//**********************************************************************
void KChINetworkFunct::CompFunc(const double & t, const double *v, double *f) const
{
	static unsigned int dynValDim = model.GetNbDynVal(0);

	model.ComputeFluxes(t);
	model.Stimulate(t);

	for (unsigned int i = 0 ; i < model.cells.size() ; ++i)
		model.cells[i]->funct->CompFunc((i == 18 and ((int)(t*1000) % 1) == 0) ? (t+1.0) : 0, v + i * dynValDim, f + i * dynValDim);
}

//********************************************************************//
//********* T S O D Y K S   M A R K R A M   S Y N A P S E S **********//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
TMSynapseFunct::TMSynapseFunct(AstroModel::TMSynapse & _syn) : 
	synapse(_syn), model(*_syn.GetModel())
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void TMSynapseFunct::CompFunc(const double & , const double *v, 
	double *f) const
{
	double & diffx       = f[0];
	double & diffu       = f[1];
	double & diffgamma   = f[2];
	const double & x     = v[0];
	const double & u     = v[1];
	const double & gamma = v[2];

	diffx     = synapse.Od * (1.0 - x);
	diffu     = synapse.Of * (- u);
	diffgamma = -synapse.Oc * gamma;
}

//********************************************************************//
//**** O P T I M   T S O D Y K S   M A R K R A M   S Y N A P S E S ***//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
TMSynapseFunctOptim::TMSynapseFunctOptim(AstroModel::TMSynapseOptim & _syn) : 
	TMSynapseFunct(_syn)
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void TMSynapseFunctOptim::CompFunc(const double & , const double *v, 
	double *f) const
{
	double & diffgamma   = f[0];
	const double & gamma = v[0];

	diffgamma = -synapse.Oc * gamma;
}

//********************************************************************//
//************** S F A L I F   N E U R O N S   M O D E L *************//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
SFALIFNeuronFunct::SFALIFNeuronFunct(AstroModel::SFALIFNeuron & _neur) : 
	neuron(_neur), model(*_neur.GetModel())
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void SFALIFNeuronFunct::CompFunc(const double & , const double *v, 
	double *f) const
{
	static double leak;
	static double tmpAMPA;

	double & diffV   = f[0];
	double & diffw   = f[1];
	const double & V = v[0];
	const double & w = v[1];

	leak = neuron.gL / neuron.C * (V - neuron.E0);
	tmpAMPA = 0.0;
	for (unsigned int i = 0 ; i < neuron.dendrSyn.size() ; ++i)
		tmpAMPA += neuron.dendrSyn[i]->GetDynVal(2) * 0.1;
	tmpAMPA = tmpAMPA / neuron.C;

	diffV = - leak - w / neuron.C + neuron.inputCurr / neuron.C;
	diffw = - w / neuron.tauw;
}

//********************************************************************//
//**************** D U M M Y   N E U R O N S   M O D E L *************//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
DummyNeuronFunct::DummyNeuronFunct(AstroModel::DummyNeuron & )
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void DummyNeuronFunct::CompFunc(const double & , const double *, 
	double *f) const
{
	// Do nothing you dummy !
	f[0] = 0.0;
}

//********************************************************************//
//************* N E U R O N   N E T W O R K  M O D E L ***************//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
NeuronNetworkFunc::NeuronNetworkFunc(AstroModel::NeuronNetModel & _model) : model(_model)
{

}

//**********************************************************************
// Function network operator, it returns the derivatives of every dynVal
//**********************************************************************
void NeuronNetworkFunc::CompFunc(const double & t, const double *v, double *f) const
{
	unsigned int nbValDyn = 0;
	for (unsigned int i = 0 ; i < model.neurons.size() ; ++i)
	{
		model.neurons[i]->funct->CompFunc(t, v + nbValDyn, f + nbValDyn);
		nbValDyn += model.neurons[i]->GetNbDynVal();
	}

	for (unsigned int i = 0 ; i < model.synapses.size() ; ++i)
	{
		model.synapses[i]->funct->CompFunc(t, v + nbValDyn, f + nbValDyn);
		nbValDyn += model.synapses[i]->GetNbDynVal();
	}
}

//**********************************************************************
// Constructor
//**********************************************************************
AstroNeuronNetFunc::AstroNeuronNetFunc(AstroModel::AstroNeuroNetModel & _model) : model(_model)
{

}

//**********************************************************************
// Function network operator, it returns the derivatives of every dynVal
//**********************************************************************
void AstroNeuronNetFunc::CompFunc(const double & t, const double *v, double *f) const
{
	model.neuronNet->function->CompFunc(t, v, f);
	unsigned int nbDynVal = model.neuronNet->GetTotNbDynVals();

	double DefaultSpillOvFract = 0.025;
	GlutamatergicSynapse *gluSyn = 0;
	TMSynapse *tmSyn = 0;
	double glu = 0;
	double Calc = 0.0;
	for (unsigned int i = 0 ; i < model.astrToSyn.size() ; ++i)
	{
		model.astroNet->cells[i]->gluIP3Prod = 0;
		Calc = v[model.astroNet->cells[i]->dynVals - model.vals + ChICell::Ca];
		for (unsigned int j = 0 ; j < model.astrToSyn[i].size() ; ++j)
		{
			gluSyn = dynamic_cast<GlutamatergicSynapse *>(model.astrToSyn[i][j]);
			if (gluSyn)
			{
				tmSyn = dynamic_cast<TMSynapse *>(gluSyn);
				glu = (tmSyn ? tmSyn->spillOvFract : DefaultSpillOvFract) * gluSyn->GetGluVal();
			}	
			model.astroNet->cells[i]->gluIP3Prod += 
				model.astroNet->cells[i]->vbeta * (pow(glu, 0.7) / 
				(pow(glu, 0.7) + pow(model.astroNet->cells[i]->kR + 
					model.astroNet->cells[i]->kP*(Calc / (Calc + model.astroNet->cells[i]->kpi)), 0.7)));
		}
	}

	model.astroNet->function->CompFunc(t, v + nbDynVal, f + nbDynVal);
}

//********************************************************************//
//************* F I R E   D I F F U S E   M O D E L ******************//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
FireDiffuseCellFunct::FireDiffuseCellFunct(AstroModel::FireDiffuseCell & _cell) : cell(_cell), model(*_cell.model)
{

}

//**********************************************************************
// Constructor
//**********************************************************************
FireDiffuseNetFunct::FireDiffuseNetFunct(AstroModel::FireDiffuseModel & _model) : model(_model)
{
}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void FireDiffuseCellFunct::CompFunc(const double &, const double *v, double *f) const
{
	double & diffC    = f[0];
	const double & C  = v[0];

	diffC = - cell.degrad * C - cell.totFlux;
}

//**********************************************************************
// Function network operator, it returns the derivatives of every dynVal
//**********************************************************************
void FireDiffuseNetFunct::CompFunc(const double & t, const double *v, double *f) const
{
	static unsigned int dynValDim = model.GetNbDynVal(0);

	model.ComputeFluxes(t);
	model.Stimulate(t);

	for (unsigned int i = 0 ; i < model.cells.size() ; ++i)
		model.cells[i]->funct->CompFunc(t, v + i * dynValDim, f + i * dynValDim);
}
