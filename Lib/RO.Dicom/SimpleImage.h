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
		// Класс для приложений, не чувствительных к большому количеству служебных полей,
		// которым нужно вывести минимальную картинку, чтобы ее мог показать элиментарный viewer.
		// Например, симулятору портальных изображений в проекте Монте-Карло.
		class RODICOM_API SimpleImage
		{
		public:
			SimpleImage(void);
			~SimpleImage(void);

			// S E T

			// Модальность изображения (RTIMAGE для DRR)
			void SetModality(const char* modality);

			// Размер матрицы изображения
			void SetResolution(int nx, int ny) { nx_ = nx; ny_ = ny; }

			// Размер ячейки изображения в см
			void SetPixelSize(double psx, double psy) { psx_ = psx; psy_ = psy; }

			// Координаты левого верхнего угла изображения в см
			// (именно угла; центр первого пиксела будет смещен на половину своего размера от этого угла)
			void SetImageCorner(double x0, double y0) { x0_ = x0; y0_ = y0; }

			// Положение среза по оси Z
			void SetPosition(double z) { z0_ = z; }

			// Матрица изображения. Ее размер не контролируется и предполагается
			// соответствующим указанным выше размерам изображения.
			// SimpleImage по определению хранит матрицу после коррекции slope/intercept.
			// а устанавливается обычно из DICOM объекта, то может потребоваться трансформация.
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

			// Запись данных включая мета-данные
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
