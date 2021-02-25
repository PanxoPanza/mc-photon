#include "MCRT_library.h"

//*****  constructor
Vector3D::Vector3D()
{
	x = 0;
	y = 0;
	z = 0;
}

Vector3D::Vector3D(double x1, double y1, double z1)
{
	x = x1;
	y = y1;
	z = z1;
}

Vector3D::Vector3D(Point3D r)     //initializing object with values.
{
	x = r.x;
	y = r.y;
	z = r.z;
}
Vector3D::Vector3D(const Vector3D &vec)
{
	x = vec.x;
	y = vec.y;
	z = vec.z;
}

void Vector3D::SetVector(double x1, double y1, double z1)
{
	x = x1;
	y = y1;
	z = z1;
}

//***** addition
Vector3D Vector3D ::operator+(const Vector3D &vec)
{
	return Vector3D(x + vec.x, y + vec.y, z + vec.z);
}
////assigning new result to the vector
Vector3D &Vector3D ::operator+=(const Vector3D &vec)
{
	x += vec.x;
	y += vec.y;
	z += vec.z;
	return *this;
}
//***** substraction//
Vector3D Vector3D ::operator-(const Vector3D &vec)
{
	return Vector3D(x - vec.x, y - vec.y, z - vec.z);
}
//assigning new result to the vector
Vector3D &Vector3D::operator-=(const Vector3D &vec)
{
	x -= vec.x;
	y -= vec.y;
	z -= vec.z;
	return *this;
}

//***** scalar multiplication
Vector3D Vector3D ::operator*(double value)
{
	return Vector3D(x*value, y*value, z*value);
}

//assigning new result to the vector.
Vector3D &Vector3D::operator*=(double value)
{
	x *= value;
	y *= value;
	z *= value;
	return *this;
}

//*****  scalar division
Vector3D Vector3D ::operator/(double value)
{
	assert(value != 0);
	return Vector3D(x / value, y / value, z / value);
}
//assigning new result to the vector
Vector3D &Vector3D ::operator/=(double value)
{
	assert(value != 0);
	x /= value;
	y /= value;
	z /= value;
	return *this;
}

Vector3D &Vector3D::operator=(const Vector3D &vec)
{
	x = vec.x;
	y = vec.y;
	z = vec.z;
	return *this;
}

//Dot product
double  Vector3D::dot(const Vector3D &vec)
{
	return x*vec.x + vec.y*y + vec.z*z;
}

//cross product
Vector3D Vector3D::cross(const Vector3D &vec)
{
	double ni = y*vec.z - z*vec.y;
	double nj = z*vec.x - x*vec.z;
	double nk = x*vec.y - y*vec.x;
	return Vector3D(ni, nj, nk);
}

double Vector3D::magnitude()
{
	return sqrt(square());
}

double Vector3D::square()
{
	return x*x + y*y + z*z;
}

Vector3D Vector3D::normalize()
{
	assert(magnitude() != 0);
	*this /= magnitude();
	return *this;
}

double Vector3D::distance(const Vector3D &vec)
{
	Vector3D dist = *this - vec;
	return dist.magnitude();
}


double Vector3D::X() const
{
	return x;
}

double Vector3D::Y() const
{
	return y;
}

double Vector3D::Z() const
{
	return z;
}

string Vector3D::Print(void)
{
	return to_string(x) + ",	" 
		 + to_string(y) + ",	" 
		 + to_string(z);
}