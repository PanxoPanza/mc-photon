#define _USE_MATH_DEFINES
#include "MCRT_library.h"

Panel::Panel() {
	Initialize();
	Allocate();
}

Panel::Panel(Point3D r0, Point3D r1, Point3D r2)
{
	Initialize();
	Allocate();
	SetCoord(r0, r1, r2); //Set coordinates, normal direction, and centroid
}

Panel::~Panel()
{
	CleanUp();
}

void Panel::Initialize(void)
{
	nVtx = 0;
	theta = 0;
	phi = 0;
	vtx = NULL;
	vtx2D = NULL;
	centroid = NULL;
	normal = NULL;
	Tol = 1E-10;
}

void Panel::CleanUp(void)
{
	if (NULL != vtx)		delete[] vtx;
	if (NULL != vtx2D)		delete[] vtx2D;
	if (NULL != centroid)	delete centroid;
	if (NULL != normal)		delete normal;
	Initialize();
}

void Panel::Allocate(void)
{
	nVtx = 3;
	vtx = new Point3D[nVtx];
	vtx2D = new Point3D[nVtx];
	centroid = new Point3D;
	normal = new Vector3D;
}

void Panel::SetCoord(Point3D r1, Point3D r2, Point3D r3) {
	vtx[0] = r1;
	vtx[1] = r2;
	vtx[2] = r3;

	// set panels reference elements
	SetPlaneNormal(); // panel's normal
	SetNormalAngles(); // angles respect to the global coordinates
	SetCentroid(); // panel's centroid

	// Get vertex in local coordinates
	vtx2D[0] = GlobaltoLocal(vtx[0]);
	vtx2D[1] = GlobaltoLocal(vtx[1]);
	vtx2D[2] = GlobaltoLocal(vtx[2]);
	A_panel = ElementArea(vtx2D[0], vtx2D[1], vtx2D[2]); // panel's area
}

void Panel::SetNormalAngles(void) {
	phi = atan2(normal->Y(), normal->X());
	Vector3D Vxy(normal->X(), normal->Y());
	theta = atan2(Vxy.magnitude(), normal->Z());
}

void Panel::SetPlaneNormal(void) //it creates a normal pointing towards de 123 rotation
{
	Vector3D v12 = Vector3D(vtx[1] - vtx[0]); // Get vector V12
	Vector3D v23 = Vector3D(vtx[2] - vtx[1]); // Get vector V23
	*normal = v12.cross(v23).normalize();  // v12 X v23
}

void Panel::SetCentroid(void) {
	*centroid = (vtx[0] + vtx[1] + vtx[2]) / 3.0;
}

double Panel::ElementArea(Point3D r1, Point3D r2, Point3D r3) const
{
	return abs((r1.x*(r2.y - r3.y)
		+ r2.x*(r3.y - r1.y)
		+ r3.x*(r1.y - r2.y)) / 2.0); //Area of a triangle
}

Point3D Panel::GlobaltoLocal(Point3D rx) const
{
	double cosT = normal->Z();
	double sinT = sqrt(1 - cosT*cosT);
	double cosP = sinT < Tol ? cosT / abs(cosT) : normal->X() / sinT;
	double sinP = sinT < Tol ? 0 : normal->Y() / sinT;

	Point3D rx0 = rx - Centroid(), rxx; // Translate
									// Rotate
	rxx.x = +cosP*cosT*rx0.x + sinP*cosT*rx0.y - sinT*rx0.z;
	rxx.y = -sinP*rx0.x + cosP*rx0.y + 0 * rx0.z;
	rxx.z = +cosP*sinT*rx0.x + sinP*sinT*rx0.y + cosT*rx0.z;
	return rxx;
}

Vector3D Panel::GlobaltoLocal(Vector3D Vx) const
{
	double cosT = normal->Z();
	double sinT = sqrt(1 - cosT*cosT);
	double cosP = sinT < Tol ? cosT / abs(cosT) : normal->X() / sinT;
	double sinP = sinT < Tol ? 0 : normal->Y() / sinT;

	return Vector3D(
		+cosP*cosT*Vx.X() + sinP*cosT*Vx.Y() - sinT*Vx.Z(),
		-sinP*Vx.X() + cosP*Vx.Y() + 0 * Vx.Z(),
		+cosP*sinT*Vx.X() + sinP*sinT*Vx.Y() + cosT*Vx.Z());
}

Point3D Panel::LocaltoGlobal(Point3D rx) const
{
	double cosT = normal->Z();
	double sinT = sqrt(1 - cosT*cosT);
	double cosP = sinT < Tol ? cosT/abs(cosT) : normal->X() / sinT;
	double sinP = sinT < Tol ? 0 : normal->Y() / sinT;

	Point3D rxx;
	rxx.x = +cosP*cosT*rx.x - sinP*rx.y + cosP*sinT*rx.z;
	rxx.y = +sinP*cosT*rx.x + cosP*rx.y + sinP*sinT*rx.z;
	rxx.z = -sinT*rx.x + 0 * rx.y + cosT*rx.z;
	return rxx + Centroid();
}

Vector3D Panel::LocaltoGlobal(Vector3D Vx) const
{
	double cosT = normal->Z();
	double sinT = sqrt(1 - cosT*cosT);
	double cosP = sinT < Tol ? cosT / abs(cosT) : normal->X() / sinT;
	double sinP = sinT < Tol ? 0 : normal->Y() / sinT;

	return Vector3D(
		+cosP*cosT*Vx.X() - sinP*Vx.Y() + cosP*sinT*Vx.Z(),
		+sinP*cosT*Vx.X() + cosP*Vx.Y() + sinP*sinT*Vx.Z(),
		-sinT*Vx.X() + 0 * Vx.Y() + cosT*Vx.Z());
}

/*bool Panel::LocalProjection(Point3D &rplane,const Photon &hw) const{
	// Gets photon's position proyected in the panel's local x-y plane
	
	// Transform position and momentum to local coordinates
	rplane = GlobaltoLocal(hw.x0());
	Vector3D hk = GlobaltoLocal(hw.k_hat());
	
	// False if photon is parallel to panel (use tolerance)
	if (abs(hk.Z()) < Tol) return false;

	// False if photon is moving away from the panel
	if (rplane.z > 0 && hk.Z() > 0) return false;
	if (rplane.z < 0 && hk.Z() < 0) return false;

	// Vertical distance to the plane
	double MoveConst = abs(rplane.z / hk.Z());

	// Position of photon when hits infinite plane (z_local = 0)
	rplane.x += MoveConst * hk.X();
	rplane.y += MoveConst * hk.Y();
	rplane.z = 0;
	return true;
}**/

bool Panel::IsInside(const Point3D &r) const
{
	if (abs(r.z) >= Tol) return false; // if point not lying at the surface

	double A1 = ElementArea(r, vtx2D[1], vtx2D[2]);
	double A2 = ElementArea(vtx2D[0], r, vtx2D[2]);
	double A3 = ElementArea(vtx2D[0], vtx2D[1], r);

	if (abs(A_panel - (A1 + A2 + A3))/A_panel <= Tol) return true;
	return false;
}

void Panel::MovePanel(const Point3D &r) {
	for (int i = 0; i < nVtx; i++) {
		vtx[i] += r;
	}
	SetCentroid();
}

// check whether photon hits a panel.
/*bool Panel::PhotonHits(const Photon &hw) const {

	Point3D pos_hw;
	if (LocalProjection(pos_hw, hw)) {
		return IsInside(pos_hw); // Check if point is inside the plane
	}
	else {
		return false;
	}

}*/

// check whether photon hits a panel.
// retrieves the photon's final position and distance to panel
bool Panel::PhotonHits(const Photon &hw, double &dXif, Point3D &Xf) const{
	Vector3D Xic(hw.x0() - Centroid()); // vector from centroid to photon's position
	double Xif_z = normal->dot(Xic); // z-local distance from photon to panel
	double kz_hat = normal->dot(hw.k_hat()); // z-local momentum of photon

	// False if photon is parallel to panel (use tolerance)
	if (abs(kz_hat) < Tol)       { dXif = DBL_MAX; return false; }

	// False if photon is moving away from the panel
	if (Xif_z > 0 && kz_hat > 0) { dXif = DBL_MAX; return false; }
	if (Xif_z < 0 && kz_hat < 0) { dXif = DBL_MAX; return false; }

	dXif = abs(Xif_z) / abs(kz_hat);
	Vector3D k_hat = hw.k_hat();
	Xf = hw.x0();// initial possition
	// move photon to the final position
	Xf.x += dXif*k_hat.X();
	Xf.y += dXif*k_hat.Y();
	Xf.z += dXif*k_hat.Z();

	Xf = GlobaltoLocal(Xf);  // get final position in local coordinates
	return IsInside(Xf); // check if inside of panel
}

double // get vertical distance from panels to photon's position
Panel::DistanceToPhoton(const Photon &hw) const
{
	/*Point3D pos_hw;
	if (LocalProjection(pos_hw, hw)){
		pos_hw = LocaltoGlobal(pos_hw); // pass proyection to local coordinates
		
		// Distance vector between the photon proyected in the plane and the photon origin
		Vector3D DV(pos_hw - hw.x0());
		return DV.magnitude();
	}
	return DBL_MAX;*/

	Vector3D Xic(hw.x0() - Centroid()); // distance from photon's position to centroid
	double Xif_z = normal->dot(Xic); // z-local distance from photon to panel
	double kz_hat = normal->dot(hw.k_hat()); // z-local k of photon

											 // False if photon is parallel to panel (use tolerance)
	if (abs(kz_hat) < Tol) return DBL_MAX;

	// False if photon is moving away from the panel
	if (Xif_z > 0 && kz_hat > 0) return DBL_MAX;
	if (Xif_z < 0 && kz_hat < 0) return DBL_MAX;

	return abs(Xif_z) / abs(kz_hat);
}

Point3D Panel::Vertex(int i) const { return vtx[i];}
Vector3D Panel::Normal(void) const {return *normal;}
Point3D Panel::Centroid(void) const {return *centroid;}
double Panel::NormalTheta(void) const {return theta; }
double Panel::NormalPhi(void) const {return phi; }

string Panel::Print(void) const {
	string output;
	output = "";
	for (int ivtx = 0; ivtx < 3; ivtx++) {
		output += "	Vertex[" + to_string(ivtx) + "] =	"
			+ vtx[ivtx].Print() + "\n";
	}
	output += "	Centroid =	" + centroid->Print() + "\n";
	output += "	Normal =	" + normal->Print() + "\n";
	return output;
}