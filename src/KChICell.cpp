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

#include "KChICell.h"
#include "KChIModel.h"

#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_const_mksa.h>

using namespace AstroModel;

// Do not use this in simulation, model not correctly calibrated yet
//********************************************************************//
//********************* K C H I   C E L L ****************************//
//********************************************************************//

double KChICell::DefaultKo  = 3.5;
double KChICell::DefaultKi  = 113.5;
double KChICell::DefaultVm  = -0.08000;//-0.09298;
double KChICell::DefaultCer = 0.0065;//0.01067; // Attempt to reduce ChI Ca spike
double KChICell::DefaultGp  = 0.5;
double KChICell::DefaultGm  = 0;

double KChICell::DefaultVa            = 6.5450e-14;
double KChICell::DefaultSa            = 7.8540e-09;
double KChICell::DefaultalphExt       = 0.15;
double KChICell::DefaultKoBl          = 3.5;
double KChICell::DefaultCaOut         = 2.0;
double KChICell::DefaultClOut         = 143.5;
double KChICell::DefaultClIn          = 4.8;
double KChICell::DefaultNaOut         = 137.6;
double KChICell::DefaultNaIn          = 30;
double KChICell::DefaultCap           = DefaultSa * 1.0e-02;
double KChICell::DefaultGClLeak       = 0.0;//5.831; // 
double KChICell::DefaultGNaLeak       = 2.75; // 
double KChICell::DefaultOmegaK        = 10;//0.5;
double KChICell::DefaultGKLeak        = 1.0;//10.0;
double KChICell::DefaultGKirMax       = 5.0;//50;//16.96;
double KChICell::DefaultKNa           = 10;// Ostby 2009
double KChICell::DefaultKK            = 1.5;// Ostby 2009
double KChICell::DefaultJNaKATPaseMax = 0; // 
double KChICell::DefaultGCaLeak       = 0;//1.9571e-11; // Depends on other parameters
double KChICell::DefaultPCaLType      = 3.0e-09;//3.0e-08;
double KChICell::DefaultOMP           = 0.010e-7;//2.0e-07;//DefaultVa * 0.0011;
double KChICell::Defaultkpact         = 300.0; // 
double KChICell::Defaultkpinh         = 0.02; // 
double KChICell::DefaultkPhos         = 1.0; // 
bool   KChICell::DefaultEnsureCaEq    = true;


std::string KChICell::ClassName = "KChICell";

//**********************************************************************
// Default Constructor
//**********************************************************************
KChICell::KChICell(const KChIModel * _model, double *_dv, bool _fv) : 
	ChICell(_model, _dv, _fv), 
	Va(DefaultVa), Sa(DefaultSa), alphExt(DefaultalphExt), KoBl(DefaultKoBl), 
	CaOut(DefaultCaOut), ClOut(DefaultClOut), ClIn(DefaultClIn), NaOut(DefaultNaOut), 
	NaIn(DefaultNaIn), Cap(DefaultCap), GClLeak(DefaultGClLeak), 
	GNaLeak(DefaultGNaLeak), OmegaK(DefaultOmegaK),	GKLeak(DefaultGKLeak), 
	GKirMax(DefaultGKirMax), KNa(DefaultKNa), KK(DefaultKK), 
	JNaKATPaseMax(DefaultJNaKATPaseMax), GCaLeak(DefaultGCaLeak), 
	PCaLType(DefaultPCaLType), OMP(DefaultOMP), kpact(Defaultkpact), 
	kpinh(Defaultkpinh), kPhos(DefaultkPhos), ensureCaEq(DefaultEnsureCaEq)
{
	nbDynVals = KCHIMODEL_NBVALS_PER_CELL;
	if (funct)
		delete funct;
	funct = new ODE::KChICellFunct(*this);
	if (not _dv)
	{
		if (dynVals)
			delete [] dynVals;
		dynVals = new double[nbDynVals];
		freeDynVals = true;
	}
	Initialize();
}

//**********************************************************************
// Copy constructor
//**********************************************************************
KChICell::KChICell(const KChICell & c) : ChICell(c), Va(c.Va), Sa(c.Sa), 
	alphExt(c.alphExt), KoBl(c.KoBl), CaOut(c.CaOut), ClOut(c.ClOut), 
	ClIn(c.ClIn), NaOut(c.NaOut), NaIn(c.NaIn), Cap(c.Cap), GClLeak(c.GClLeak), 
	GNaLeak(c.GNaLeak), OmegaK(c.OmegaK), GKLeak(c.GKLeak), GKirMax(c.GKirMax), 
	KNa(c.KNa), KK(c.KK), JNaKATPaseMax(c.JNaKATPaseMax), GCaLeak(c.GCaLeak), 
	PCaLType(c.PCaLType), OMP(c.OMP), kpact(c.kpact), kpinh(c.kpinh), 
	kPhos(DefaultkPhos), ensureCaEq(c.ensureCaEq)
{
	if (funct)
		delete funct;
	funct = new ODE::KChICellFunct(*this);
	Initialize();
}

//**********************************************************************
// Destructor
//**********************************************************************
KChICell::~KChICell()
{
}

//**********************************************************************
//**********************************************************************
void KChICell::Initialize()
{
	ChICell::Initialize();

	KFluxOut = 0.0;
	KFluxIn = 0.0;

	dynVals[Ko]  = DefaultKo;
	dynVals[Ki]  = DefaultKi;
	dynVals[Vm]  = DefaultVm;
	dynVals[Cer] = DefaultCer;
	dynVals[Gp]  = DefaultGp;
	dynVals[Gm]  = DefaultGm;

	if (defaultBiophysParams)
	{
		Va            = DefaultVa;
		Sa            = DefaultSa;       
		alphExt       = DefaultalphExt;   
		KoBl          = DefaultKoBl;      
		CaOut         = DefaultCaOut;     
		ClOut         = DefaultClOut;    
		ClIn          = DefaultClIn;      
		NaOut         = DefaultNaOut;     
		NaIn          = DefaultNaIn;      
		Cap           = DefaultCap;       
		GClLeak       = DefaultGClLeak;
		GNaLeak       = DefaultGNaLeak;
		OmegaK        = DefaultOmegaK;    
		GKLeak        = DefaultGKLeak;    
		GKirMax       = DefaultGKirMax;   
		KNa           = DefaultKNa;
		KK            = DefaultKK;
		JNaKATPaseMax = DefaultJNaKATPaseMax;
		GCaLeak       = DefaultGCaLeak;   
		PCaLType      = DefaultPCaLType;  
		OMP           = DefaultOMP;       
		kpact         = Defaultkpact;     
		kpinh         = Defaultkpinh;     
		kPhos         = DefaultkPhos;
	}

	ComputeOtherParameters();
}

//**********************************************************************
// Compute other parameters
//**********************************************************************
void KChICell::ComputeOtherParameters()
{
	VolCyt = Va / (1.0 + c1);
	VolER = c1 * VolCyt;
	VolExt = alphExt * Va;
	
	LogForECl = gsl_sf_log(ClOut / ClIn);
	LogForENa = gsl_sf_log(NaOut / NaIn);

	OCK = kpact * model->v3k;
	OPK = kpinh * this->kP / this->kR;
	
	if (ensureCaEq)
	{
		const KChIModel *tmp = dynamic_cast<const KChIModel*>(model);
		assert(tmp);
		double ECa = gsl_sf_log(CaOut / DefaultCa) / (2.0 * tmp->FoRT );
		double mInf = 1.0 / (1.0 + gsl_sf_exp((DefaultVm - tmp->VLTmHalf) / tmp->VLTmSlope));
		double hInf = 1.0 / (1.0 + gsl_sf_exp((DefaultVm - tmp->VLThHalf) / tmp->VLThSlope));
		double IPMCA = 2.0 * GSL_CONST_MKSA_FARADAY * OMP * HILL2(DefaultCa, tmp->KMP);
		double LTv = 2.0 * DefaultVm * tmp->FoRT;
		double ILType = mInf * hInf * 2.0 * PCaLType * GSL_CONST_MKSA_FARADAY * LTv * 
			(DefaultCa - CaOut * gsl_sf_exp(-LTv)) / (1.0 - gsl_sf_exp(-LTv));
		GCaLeak = - (IPMCA + ILType) / (DefaultVm - ECa);
		TRACE("GCaLeak : " << GCaLeak)

		double Ek = gsl_sf_log(DefaultKo / DefaultKi) / tmp->FoRT;
		double GKir = GKirMax / sqrt(DefaultKo * (1.0 + 
			gsl_sf_exp((DefaultVm - tmp->VKirH - Ek) / tmp->VKirS)));
		JNaKATPaseMax = (GKLeak + GKir) * (DefaultVm - Ek) / 
			(2.0 * GSL_CONST_MKSA_FARADAY * HILLn(NaIn, KNa, 1.5) * HILL1(DefaultKo, KK));
		TRACE("JNaKATPaseMax : " << JNaKATPaseMax)

		double ENa = gsl_sf_log(DefaultNaOut / DefaultNaIn) / tmp->FoRT;
		GNaLeak = 3.0 * GSL_CONST_MKSA_FARADAY * JNaKATPaseMax *
			HILLn(NaIn, KNa, 1.5) * HILL1(DefaultKo, KK) / (ENa - DefaultVm);
		TRACE("GNaLeak : " << GNaLeak)

		// SERCA pumps activity
		
		mInf  = HILL1(DefaultIP3, model->d1) * HILL1(DefaultCa, model->d5);
		double Jchan = VolER * rC  * pow(mInf * Defaulth, 3) * (DefaultCer - DefaultCa);
		double Jleak = VolER * rL  * (DefaultCer - DefaultCa);
		this->vER = (Jchan + Jleak) / (VolCyt * HILL2(DefaultCa, this->Ker));
	}
}

//**********************************************************************
//**********************************************************************
void KChICell::SetToEquilibrium()
{
	ChICell::SetToEquilibrium();

	KFluxOut = 0.0;
	KFluxIn = 0.0;

	dynVals[Ko]  = DefaultKo;
	dynVals[Ki]  = DefaultKi;
	dynVals[Vm]  = DefaultVm;
	dynVals[Cer] = DefaultCer;
	dynVals[Gp]  = DefaultGp;
	dynVals[Gm]  = DefaultGm;
}

//**********************************************************************
//**********************************************************************
bool KChICell::LoadFromStream(std::ifstream & stream)
{
	bool ok = ChICell::LoadFromStream(stream);
	
	stream >> Va;
	stream >> Sa;      
	stream >> alphExt; 
	stream >> KoBl;    
	stream >> CaOut;   
	stream >> ClOut;   
	stream >> ClIn;    
	stream >> NaOut;   
	stream >> NaIn;    
	stream >> Cap;     
	stream >> GClLeak;
	stream >> GNaLeak;
	stream >> OmegaK;  
	stream >> GKLeak;  
	stream >> GKirMax; 
	stream >> KNa;
	stream >> KK;
	stream >> JNaKATPaseMax;
	stream >> GCaLeak; 
	stream >> PCaLType;
	stream >> OMP;     
	stream >> kpact;   
	stream >> kpinh;   
	stream >> kPhos;
	stream >> ensureCaEq;

	ComputeOtherParameters();

	return ok and stream.good();
}

//**********************************************************************
//**********************************************************************
bool KChICell::SaveToStream(std::ofstream & stream) const
{
	bool ok = ChICell::SaveToStream(stream);

	stream
		<< Va            << std::endl
		<< Sa            << std::endl
		<< alphExt       << std::endl
		<< KoBl          << std::endl
		<< CaOut         << std::endl
		<< ClOut         << std::endl
		<< ClIn          << std::endl
		<< NaOut         << std::endl
		<< NaIn          << std::endl
		<< Cap           << std::endl
		<< GClLeak       << std::endl
		<< GNaLeak       << std::endl
		<< OmegaK        << std::endl
		<< GKLeak        << std::endl
		<< GKirMax       << std::endl
		<< KNa           << std::endl
		<< KK            << std::endl
		<< JNaKATPaseMax << std::endl
		<< GCaLeak       << std::endl
		<< PCaLType      << std::endl
		<< OMP           << std::endl
		<< kpact         << std::endl
		<< kpinh         << std::endl
		<< kPhos         << std::endl
		<< ensureCaEq    << std::endl;

	return ok and stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler KChICell::BuildModelParamHandler()
{
	ParamHandler params;
	params += ChICell::BuildModelParamHandler();

	params <= "DefaultKoVal" , DefaultKo;
	params <= "DefaultKiVal" , DefaultKi;
	params <= "DefaultVmVal" , DefaultVm;
	params <= "DefaultCerval", DefaultCer;
	params <= "DefaultGpVal" , DefaultGp;
	params <= "DefaultGmVal" , DefaultGm;

	params <= "Va"           , DefaultVa;     
	params <= "Sa"           , DefaultSa;     
	params <= "alphExt"      , DefaultalphExt;
	params <= "KoBl"         , DefaultKoBl; 
	params <= "CaOut"        , DefaultCaOut;
	params <= "ClOut"        , DefaultClOut;
	params <= "ClIn"         , DefaultClIn; 
	params <= "NaOut"        , DefaultNaOut;
	params <= "NaIn"         , DefaultNaIn; 
	params <= "Cap"          , DefaultCap;
	params <= "GClLeak"      , DefaultGClLeak;
	params <= "GNaLeak"      , DefaultGNaLeak;
	params <= "OmegaK"       , DefaultOmegaK; 
	params <= "GKLeak"       , DefaultGKLeak; 
	params <= "GKirMax"      , DefaultGKirMax;
	params <= "KNa"          , DefaultKNa;
	params <= "KK"           , DefaultKK;
	params <= "JNaKATPaseMax", DefaultJNaKATPaseMax;
	params <= "GCaLeak"      , DefaultGCaLeak; 
	params <= "PCaLType"     , DefaultPCaLType;
	params <= "OMP"          , DefaultOMP;     
	params <= "kpact"        , Defaultkpact;
	params <= "kpinh"        , Defaultkpinh;
	params <= "kPhos"        , DefaultkPhos;
	params <= "ensureCaEq"   , DefaultEnsureCaEq;

	return params;
}

//**********************************************************************
// Update value names in ODEProblem
//**********************************************************************
void KChICell::SetValNamesPostfix(std::string pf)
{
	ChICell::SetValNamesPostfix(pf);
	model->AddPostfixToValName(dynVals + Ko, pf + "_Ko");
	model->AddPostfixToValName(dynVals + Ki, pf + "_Ki");
	model->AddPostfixToValName(dynVals + Vm, pf + "_Vm");
	model->AddPostfixToValName(dynVals + Cer, pf + "_Cer");
	model->AddPostfixToValName(dynVals + Gp, pf + "_Gp");
	model->AddPostfixToValName(dynVals + Gm, pf + "_Gm");
}

