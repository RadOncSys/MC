#include "mcVRMLDumper.h"

void mcVRMLDumper::dumpHead(ostream& os)
{
	os << "#VRML V2.0 utf8" << endl;
}

void mcVRMLDumper::dumpWorldAxis(ostream& os, double a)
{
	double fnt_size = 0.1 * a;

	os << "# World coordinate system axes" << endl;
	for (int i = 0; i < 3; i++)
	{
		// Линия
		os << "Shape {" << endl;
		os << "  appearance Appearance {" << endl;
		os << "    material Material {" << endl;
		//os << "      emissiveColor 1 1 1" << endl;
		os << "      emissiveColor ";
		if (i == 0) os << "1 0 0"; else if (i == 1) os << "0 1 0"; else os << "0 0 1";
		os << endl;
		os << "      transparency 0" << endl;
		os << "    }" << endl;
		os << "  }" << endl;
		os << "  geometry IndexedLineSet {" << endl;
		os << "    coord Coordinate {" << endl;
		os << "      point [" << endl;
		os << "        ";
		if (i == 0) os << -0.5*a << " 0 0 " << endl;
		else if (i == 1) os << "0 " << -0.5*a << " 0 " << endl;
		else os << "0 0 " << -a << endl;
		os << "        ";
		if (i == 0) os << a << " 0 0 " << endl;
		else if (i == 1) os << "0 " << a << " 0 " << endl;
		else os << "0 0 " << a << endl;
		os << "      ]" << endl;
		os << "    }" << endl;
		os << "    coordIndex [" << endl;
		os << "      0 1 " << endl;
		os << "    ]" << endl;
		os << "  }" << endl;
		os << "}" << endl;
		// Текст
		os << "Transform {" << endl;
		os << "  translation ";
		if (i == 0) os << a << " 0 0 " << endl;
		else if (i == 1) os << "0 " << a << " 0 " << endl;
		else os << "0 0 " << a << endl;
		os << "  children Billboard {" << endl;
		os << "    children Shape {" << endl;
		os << "      appearance Appearance {" << endl;
		os << "        material Material {" << endl;
		//os << "          emissiveColor 1 1 1" << endl;
		os << "          emissiveColor ";
		if (i == 0) os << "1 0 0"; else if (i == 1) os << "0 1 0"; else os << "0 0 1";
		os << endl;
		os << "        }" << endl;
		os << "      }" << endl;
		os << "      geometry Text {" << endl;
		os << "        string [\"";
		if (i == 0) os << 'X'; else if (i == 1) os << 'Y'; else os << 'Z';
		os << "\"]" << endl;
		os << "        fontStyle FontStyle {" << endl;
		os << "          size " << fnt_size << endl;
		os << "        }" << endl;
		os << "      }" << endl;
		os << "    }" << endl;
		os << "  }" << endl;
		os << "}" << endl;
	}
}
