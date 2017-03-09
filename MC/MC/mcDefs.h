// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

#ifndef TWOPI
#define  TWOPI  6.28318530717958648
#endif

#ifndef PI5D2
#define  PI5D2  7.85398163397448310
#endif

#ifndef ABS
#define ABS(x)		  (((x) > 0) ? (x) : -(x))
#endif

#ifndef MIN
#define MIN(x,y)	  (((x) > (y)) ? (y) : (x))
#endif

#ifndef MAX
#define MAX(x,y)	  (((x) > (y)) ? (x) : (y))
#endif

#ifndef ROUND
#define	ROUND(x)		(((x)>=0.0)?int((x)+0.5):int((x)-0.5))
#endif

#ifndef SQUARE
#define  SQUARE(a)  ((a)*(a))
#endif

// Source type mainly correspond to radiation type
enum mcas_type_e { MCA_CO60, MCA_PHOTON, MCA_ELECTRON, MCA_XRAY };

// Note: The source for the following constants is:
// E. Richard Cohen and Barry N. Taylor, "The Fundamental Physical
// Constants," Physics Today, August 1998 Buyers' Guide, pp. BG7-11.
#define  MEV_PER_JOULE (1.60217733e-13)
#define  EMASS         (0.51099906)
#define  TWICE_EMASS   (1.02199812)
#define  EMASS_SQUARED (0.26112004)

// масса протона в ћэ¬ из rpp-2006-book.pdf, стр. 71, и 1/PMASS
#define  PMASS   (938.27203)
#define  PMASS_1	(0.0010657889908537505908600941669337)

#define	NAVOGADRO         (6.0221415E+23)

// —корость света (м/с)
#define	LIGHT_SPEED			299792458