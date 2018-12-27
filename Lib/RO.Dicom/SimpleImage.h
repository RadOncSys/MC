// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include "RO.Dicom.h"
#include <iostream>

namespace RO
{
	namespace Dicom
	{
		// ����� ��� ����������, �� �������������� � �������� ���������� ��������� �����,
		// ������� ����� ������� ����������� ��������, ����� �� ��� �������� ������������ viewer.
		// ��������, ���������� ���������� ����������� � ������� �����-�����.
		class RODICOM_API SimpleImage
		{
		public:
			SimpleImage(void);
			~SimpleImage(void);

			// S E T

			// ����������� ����������� (RTIMAGE ��� DRR)
			void SetModality(const char* modality);

			// ������ ������� �����������
			void SetResolution(int nx, int ny) { nx_ = nx; ny_ = ny; }

			// ������ ������ ����������� � ��
			void SetPixelSize(double psx, double psy) { psx_ = psx; psy_ = psy; }

			// ���������� ������ �������� ���� ����������� � ��
			// (������ ����; ����� ������� ������� ����� ������ �� �������� ������ ������� �� ����� ����)
			void SetImageCorner(double x0, double y0) { x0_ = x0; y0_ = y0; }

			// ��������� ����� �� ��� Z
			void SetPosition(double z) { z0_ = z; }

			// ������� �����������. �� ������ �� �������������� � ��������������
			// ��������������� ��������� ���� �������� �����������.
			// SimpleImage �� ����������� ������ ������� ����� ��������� slope/intercept.
			// � ��������������� ������ �� DICOM �������, �� ����� ������������� �������������.
			void SetPixelData(unsigned short* data, double rs = 1.0, double ri = 0.0);

			// G E T

			unsigned short* GetPixelData() { return pixels_; }
			int Nx() const { return nx_; }
			int Ny() const { return ny_; }
			double Psx() const { return psx_; }
			double Psy() const { return psy_; }
			double PosX() const { return x0_; }
			double PosY() const { return y0_; }
			double PosZ() const { return z0_; }

			// IO

			void ReadFromStream(std::istream& is);
			void ReadFromFile(const char* fname);

			// ������ ������ ������� ����-������
			void WtiteToStream(std::ostream& os) const;
			void WtiteToFile(const char* fname) const;

		private:
			std::string modality_;
			int nx_;
			int ny_;
			double psx_;
			double psy_;
			double x0_;
			double y0_;
			double z0_;
			unsigned short* pixels_;
		};
	}
}
