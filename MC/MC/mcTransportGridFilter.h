// Radiation Oncology Monte Carlo open source project
//
// Author: [2024] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransportPrism.h"

// Класс транспорта в гребенчатом фильтре, представляущем 2D матрицу пирамид, 
// состоящих из корпичиков прямоугольного сечения одинаковой высоты, но разного размера.
// Изначально класс разрабатывался для транспорта протонов в проекте ОКО.
// Наследуем Prism с целью использования ее методов для обработки столкновений с данным обектом и выхода из него.
class mcTransportGridFilter : public mcTransportPrism
{
public:
	mcTransportGridFilter(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
		int nx, int ny, int nz, double psx, double psy, double psz);
	virtual ~mcTransportGridFilter(void);

	// В отличие от стандартного начала переносит частицу на поверхность фантома 
	// и вычисляет индекс стартовой ячейки
	void beginTransport(mcParticle& p) override;

	// Старт транспорта может вызываться как в объекте, вложенном во внутреннюю структуру - воздушеый слой.
	void beginTransportInside(mcParticle& p) override;

	// Виртуальная функция перемещения частицы полностью покрывает специфику транспорта в сетке.
	mc_move_result_t moveParticle(mcParticle* particle, double& step, double& edep) override;

	// Специфичные параметры
	void setMedia(short inb, short outb) { inBrickIdx_ = inb; outBrickIdx_ = outb; }
	void setBrickSize(int iz, double x, double y) { bx_[iz] = x; by_[iz] = y; }

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	int getIdxAtPoint(const geomVector3D& p, short* pgidxx, bool& isInBrick) const;
	double getDistanceInsideVoxel(const mcParticle& particle, short* gidxNext, int& idx, bool& isHitCell);

	// Размер матрицы пирамид
	int nx_, ny_, nz_;

	// Шаг между пирамидами
	double psx_, psy_, psz_;

	// Стартовый угол стартовой пирамиды
	double x0_, y0_, z0_;

	// Размеры брикетиков пирамид по слоям 
	std::vector<double> bx_;
	std::vector<double> by_;

	// Индексы сред внутри и снаружи брикетиков
	short inBrickIdx_;
	short outBrickIdx_;
};
