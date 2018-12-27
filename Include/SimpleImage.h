// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include "RO.Dicom.h"
#include <iostream>
#include <string>

namespace RO
{
	namespace Dicom
	{
		//  ласс дл€ приложений, не чувствительных к большому количеству служебных полей,
		// которым нужно вывести минимальную картинку, чтобы ее мог показать элиментарный viewer.
		// Ќапример, симул€тору портальных изображений в проекте ћонте- арло.
		class RODICOM_API SimpleImage
		{
		public:
			SimpleImage(void);
			~SimpleImage(void);

			// S E T

			// ћодальность изображени€ (RTIMAGE дл€ DRR)
			void SetModality(const char* modality);

			// –азмер матрицы изображени€
			void SetResolution(int nx, int ny) { nx_ = nx; ny_ = ny; }

			// –азмер €чейки изображени€ в см
			void SetPixelSize(double psx, double psy) { psx_ = psx; psy_ = psy; }

			//  оординаты левого верхнего угла изображени€ в см
			// (именно угла; центр первого пиксела будет смещен на половину своего размера от этого угла)
			void SetImageCorner(double x0, double y0) { x0_ = x0; y0_ = y0; }

			// ѕоложение среза по оси Z
			void SetPosition(double z) { z0_ = z; }

			// ћатрица изображени€. ≈е размер не контролируетс€ и предполагаетс€
			// соответствующим указанным выше размерам изображени€.
			void SetPixelData(unsigned short* data);

			// G E T

			unsigned short* GetPixelData() { return pixels_; }
			int Nx() const { return nx_; }
			int Ny() const { return ny_; }
			double PosX() const { return x0_; }
			double PosY() const { return y0_; }
			double PosZ() const { return z0_; }

			// IO

			void ReadFromStream(std::istream& is);
			void ReadFromFile(const char* fname);

			// «апись данных включа€ мета-данные
			void WtiteToStream(std::ostream& os);
			void WtiteToFile(const char* fname);

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
