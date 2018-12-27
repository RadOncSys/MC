// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcParticle.h"
#include "mcObj.h"

class mcRng;
class mcThread;
class mcScoreTrack;

// ������� ����� �������� ��������� ������ ��� 
// ��������������� �������� ������� �����-�����
class mcSource : public mcObj
{
public:
	mcSource();
	mcSource(const char* name, int nThreads);
	virtual ~mcSource();

	void setScoreTrack(double R, double Z1, double Z2, double EMIN, bool doPhotons, bool doElectrons, bool doPositrons);
	bool IsGamma() const { return isGamma_; };
	bool IsStartInside() const { return isStartInside_; };
	void SetIsStartInside(bool b) { isStartInside_ = b; };

	void setModuleName(const char* s);
	const char* getModuleName() const { return attached_module_; }

	virtual void sample(mcParticle& p, mcThread* thread) = 0;

	double etotal() const;

	virtual void dumpVRML(ostream& os) const;

	friend ostream& operator << (ostream& os, const mcSource& s)
	{
		os << (const mcObj&)s;
		os << "Etotal = " << s.etotal() << endl;
		return os;
	}

protected:
	char attached_module_[64];
	// ��������� �� ������, ���������� ����������� ������ ����� ��������
	mcScoreTrack* trackScore_;
	int nThreads_;
	double* etotal_;

	// ���� ��� ����� ��������, �� ������� ��������� ������ �������.
	// �� ��������� ����� ����� ����������� ������� ���������������� �������
	// ������� ���������� �������, ��� �������.
	bool isGamma_;
	
	bool isStartInside_;
};
