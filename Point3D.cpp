#include "mcphoton_lib.h"
// ********************   Point3D class ******************
Point3D::Point3D()
{
	x = 0; y = 0; z = 0;
}
Point3D::Point3D(double x1, double y1, double z1)
{
	SetPoint3D(x1,y1,z1);
}

void Point3D::SetPoint3D(double x1, double y1, double z1){
	x = x1;
	y = y1;
	z = z1;
}


Point3D Point3D ::operator+(const Point3D &coord)
{
	return Point3D(x + coord.x, y + coord.y, z + coord.z);
}
////assigning new result to the vector
Point3D &Point3D ::operator+=(const Point3D &coord)
{
	x += coord.x;
	y += coord.y;
	z += coord.z;
	return *this;
}
//***** substraction//
Point3D Point3D ::operator-(const Point3D &coord)
{
	return Point3D(x - coord.x, y - coord.y, z - coord.z);
}
//assigning new result to the vector
Point3D &Point3D::operator-=(const Point3D &coord)
{
	x -= coord.x;
	y -= coord.y;
	z -= coord.z;
	return *this;
}

Point3D& Point3D::operator *=(double t)
{
	x *= t;
	y *= t;
	z *= t;
	return (*this);
}

Point3D& Point3D::operator /(double t)
{
	double f = 1.0F / t;
	x *= f;
	y *= f;
	z *= f;
	return (*this);
}

double Point3D::distance(const Point3D &coord) {
	return sqrt((x - coord.x)*(x - coord.x) 
		      + (y - coord.y)*(y - coord.y)
		      + (z - coord.z)*(z - coord.z));
}

string Point3D::Print(void)
{
	return to_string(x) + ",	"
		+ to_string(y) + ",	"
		+ to_string(z);
}

