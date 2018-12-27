// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcParticle.h"
#include "mcObj.h"
#include "mcScore.h"
#include <vector>

class mcMedia;
class mcRng;
class mcScore;
class mcMediumXE;

#ifndef NNEG
#define  NNEG(a)    ((a) >= 0  ? (a) : -(a))
#endif

#ifndef ZORP
#define  ZORP(a)    ((a) >= 0  ? (a) : 0)
#endif

#ifndef MAX
#define  MAX(a, b)  ((a) > (b) ? (a) :  (b))
#endif

#ifndef MIN
#define  MIN(a, b)  ((a) < (b) ? (a) :  (b))
#endif

#define MC_EPSILON     1e-6 // ������������ � ������ � �������� ����������

enum mc_move_result_t { MCMR_INTERUCT = 0, MCMR_EXIT, MCMR_CONTINUE, MCMR_DISCARGE };

/// <summary>
/// ������� ����� ����������, ���������� �������� � �������
/// ����� ������ ��������� ������� �����-�����.
///
/// ������� ���������� ������������, ��� � ������� ������� ��������� ��� �������
/// ��������� ����� ��� Z, ������, ���������� Z ���� �������� �� ������������.
/// ��� ��������� ������������� ��������� ������� ���� �����.
///
/// ��-������ ����� ��������� �������� ������� �������� / �����������.
/// ��� ��������� �������� ������� ��������� ���������� � ���������� ��� ��������� 
/// ������ � ����������� �� ����������� �������� ����� ��� Z.
///
/// ��-������, ������� ��������� ������� ������������� � ������� ��������� �������,
/// ��� �������� ������� ��������� � scoring.
/// 
/// GG 20150227
/// ������� ����� ������� ���������� ��������, ��������� ���� � ����� ��� ������� ��������.
/// ����� ������������ ��������� ���������� � ������������� ����� ������ �� ��������� � ���������� ������.
/// ���������� �� ������� ������������ ��� ���������� �� ������� ����������� � �� ���������� ������� ���� ������� �������.
/// ����� �������, ������, ������������ ������ ����� ������������� ��� �����, ������ ������� ��������� ������ �����.
/// ������� � �������� � ������ ������ ����������� �� ��������� ����, ���� ����������� ����� ����������.
/// � �������� ����� ������ ��������� �������� ���������� ������ ��� ������ �� ������ ����������� ������������.
/// ������� �� ��������� ������� ������� ����������� � XML ����� ������������ ������.
/// ������, �� ��������� ������� ���������� ������� ������ �������� ������� ��������� ����� ����� ����������� � ���������������� ������.
/// ������� ����� ��������� ������������� ������� ����������� ������� ��� ������������ ���������� ������� ������ ��������� ������.
/// ������� � ���, ��� ������������� ������������� ��������� � ���������� �������.
/// </summary>
class mcTransport : public mcObj
{
public:
	mcTransport();
	mcTransport(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcTransport(void);

	//
	// �������������
	//

	/// <summary>
	/// ���������������� ������� � ������� ������� ���������.
	/// ��� ������������ �����, ��� ��������������� ������� �������������� ���������.
	/// </summary>
	void setPosition(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);

	/// <summary>
	/// ������ �� ������, ���������� ��� �������������� ���� � ������ �� ���������� �������
	/// </summary>
	void setMediaRef(const mcMedia* media);

	/// <summary>
	/// �������, ��������������� ��� ��������� ����� ����������, ���������� �� ������ ����������� ����.
	/// ����������� ������� ����������� ������� ���������� ��������� ����� ��� ����� �������.
	/// ��� ����� ���� �������������� ������. 
	/// ���� ����������� ��������� ������, ����� ��������� ����� ��������� ���� �������.
	/// </summary>
	void setMediumMono(const char* mname);

	/// <summary>
	/// ���������� ��������� ����� �������, �� ��������� � ���������� ���������� � ���������� ��������
	/// </summary>
	void setMediumMono(double d) { defdensity_ = d; }

	/// <summary>
	/// ��������� ������ �������������� ���������, ���������� ��������� ��� ��������� 
	/// ������� ��������� �������� ������������ ������
	/// </summary>
	void MoveToCoordinateSystem(const geomMatrix3D& m);

	//
	// ���������
	//

	/// <summary>
	/// ��������� ������� ����������� ����� ������� ������� ���������� ���������� 
	/// � ������� beginTransport, � ������� ������� ���������� �������������� ���������
	/// �� ������� ������� � ��������� � ����������� ������.
	/// ��������� ������� ����� �������, ������� �� ����� �������� ����������.
	/// </summary>
	virtual void beginTransport(mcParticle& p);

	/// <summary>
	/// ��������� ������� ����������� ������ ������� ������� ���������� ����������. 
	/// � ������� �� ������� beginTransport ������� ����� ����������� ������ � ��� �� ����������� �� �����������.
	/// </summary>
	virtual void beginTransportInside(mcParticle& p);

	/// <summary>
	/// ����� ������� �������� ������, �� ���������� � ����������� ������������� ������� � ������� �������.
	/// �����, ���������� ����, ������ ���������� �������� �������.
	/// ���� ���� � ������� �������� �������.
	/// </summary>
	virtual void endTransport(mcParticle* particle);

	/// <summary>
	// ��������� ������� ������� � ��������� ������. 
	/// ��������� �� ���������� � ����������� �������, ��� ��� �������� �������� �������������.
	/// ���������������� �������������� ����� ����������� ������.
	/// </summary>
	static void simulate(mcThread* thread);

	/// <summary>
	/// ���������� ������� �� �������� ���������� (���������� � ����� �������) ������� ��������.
	/// ���������� ����� ���������� �������� �������� � ����������� ������.
	/// ����������� ������� ���������� ����� �������������� ������, ����������� ������������
	/// ������� (�.�., ���������� �� ������ �������) � ���������� ��������.
	/// � ���������� ������������ �������� ��� � ���������� �� ���� �������
	/// (����������� ������ ���������� ������).
	/// </summary>
	virtual mc_move_result_t moveParticle(mcParticle* particle, double& step, double& edep);

	//
	// Scoring
	//

	/// <summary>
	/// ����� ���������� ���������� ������ score � �����������.
	/// ��� ����� ������� � � ������������� ��������� score 
	/// ������� ��������� �� new � �� �������� delete �� ���������!!!
	/// </summary>
	void setScore(mcScore* score);
	mcScore* getScore() { return score_; }

	/// <summary>
	/// ������������� �������, ������� ������ ����������� �������, ������������ � ������ ����������.
	/// </summary>
	void setDefaultScore(int nThreads);

	/// <summary>
	/// ���������� ������� (�������) � ������.
	/// �������� �������� ����� ������ ����������.
	/// ������ ������ ���� ������� � �������� ��������� ������������ �������� ��������� ������� ������.
	/// ������� ���������� ���� ���������� ����������� �������������� �� ����������� ������� ��������� ������������ �������.
	/// </summary>
	void addRegion(mcTransport* region);

	//
	// ��������������� ������
	//
	void setPreviousTransport(mcTransport* t) { previousTransport_ = t; }
	void setNextTransport(mcTransport* t) { nextTransport_ = t; if (t != nullptr) t->setPreviousTransport(this); }
	const mcTransport* getPreviousTransport() const { return previousTransport_; }
	const mcTransport* getNextTransport() const { return nextTransport_; }

	void setInternalTransport(mcTransport* t);
	void setExternalTransport(mcTransport* t);
	mcTransport* getInternalTransport() const { return internalTransport_; }
	mcTransport* getExternalTransport() const { return externalTransport_; }

	// ������ ���������� ���������� �� �����
	virtual mcTransport* getInternalTransportByName(const char* name);

	void setStamp(int stamp) { stamp_ = stamp; }

	// �������������� ���������
	const geomMatrix3D& MW2T() { return mwtot_; }
	const geomMatrix3D& MT2W() { return mttow_; }

	// ������������ ����
	static double HowManyMFPs(mcRng& rng);

	double etotal() const;

	static void implementException();

	virtual void dump(ostream& os) const;
	virtual void dumpVRML(ostream& os) const;

	void dumpVRMLRing(ostream& os, double r1, double r2, double z, bool normPositive, double x0 = 0, double y0 = 0) const;
	void dumpVRMLSemiCircle(ostream& os, double r, double z, bool normPositive, const geomMatrix3D& M) const;
	void dumpVRMLCylinderSide(ostream& os, double r, double z1, double z2, bool normPositive, double x0 = 0, double y0 = 0) const;
	void dumpVRMLCylinderSemiSide(ostream& os, double r, double h, const geomMatrix3D& M) const;
	void dumpVRMLConicalCylinderSide(ostream& os, double r, double z1, double z2, double f, bool normPositive, double x0 = 0, double y0 = 0) const;
	void dumpVRMLCylinderRing(ostream& os, double r1, double r2, double z1, double z2) const;
	void dumpVRMLConicalRing(ostream& os, double r1, double r2, double z1, double z2, double f) const;
	void dumpVRMLCylinderWithConicalHole(ostream& os, double r1, double r2, double z1, double z2, double f) const;
	void dumpVRMLRectangleRing(ostream& os, double x1, double x2, double y1, double y2, double d, double h) const;
	void dumpVRMLCylinder(ostream& os, double r, double z1, double z2, double x0 = 0, double y0 = 0) const;
	void dumpVRMLPrism(ostream& os, double ax, double ay, double az) const;
	void dumpVRMLPolygonCircle(ostream& os, const std::vector<double>& pz, const std::vector<double>& pr) const;

	//
	// ��������� (������ ����� ��� ��������� ����������� ������� moveParticle)
	//
	virtual double getDistanceInside(mcParticle& p) const;
	virtual double getDistanceOutside(mcParticle& p) const;
	virtual double getDNearInside(const geomVector3D& p) const;

	//
	// Auxilary
	//
	short getDefMedIdx() const { return defmedidx_; }
	double getDefDensity() const { return defdensity_; }

	// GG 20140927 �������� ������������� ��������� ���������� ����� ������ ����������� ������ �� ������ ������� ����������.
	// ������� ����� ������ �� ����� ���������� � ��� ����� ������������ �� ����� ��������������.
	// �� ��������� ������ ��������������� � 0.
	// � ���� ������ ������������ ����� �� Media ��� � ������.
	// ������ ���������� ������ ���������� ������������ �� ������ ���� � ������ ���������� � ���������� � ��� ���� � ������ �������.
	// � ������ ������� ��� ������������ �����, ��������� � ������������� ���� ���������� � ������� ��� �� ����.
	// � ����������� �������� ��� ������ ������� Discharge.
	double transCutoff_phot; // Energy cutoff for photon transport
	double transCutoff_elec; // Energy cutoff for electron transport

protected:
	virtual mcMediumXE* getParticleMedium();
	virtual double getRegionRelDensity();

protected:
	mcScore* score_;
	mcTransport* previousTransport_;
	mcTransport* nextTransport_;
	mcTransport* internalTransport_;
	mcTransport* externalTransport_;

	geomMatrix3D mwtot_;	// �������������� �� ������� ������� � ������� �������
	geomMatrix3D mttow_;	// �������������� �� ������� ������� � ������� �������
	geomMatrix3D mtoe_;		// �������������� �� ������� ������� �� ���������

	const mcMedia* media_;
	short defmedidx_;
	double defdensity_;

	// �������������� ������ ������� ������� ����������
	int stamp_;

	bool isMultiRegions_;
	std::vector<mcTransport*> regions_;
};
