// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

class mcTransport;
class mcScore;
class mcSource;
class mcMedia;
class XPRNode;

class GeometryParser
{
public:
	static enum mc_particle_t convert_S2T_ptype(const wchar_t* st);
	static enum spectrum_distr_t convert_spec_type(const wchar_t* st);
	static enum profile_distr_t convert_profile_distr_type(const wchar_t* st);

	static mcTransport* ParseTransport(const XPRNode& geometry, const mcMedia* media, int nThreads);
	static mcScore* ParseScore(const XPRNode& item, int nThreads);
	static mcSource* ParseSource(const XPRNode& item, int nThreads);
};
