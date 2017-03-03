
#include "stdafx.h"
#include "cameras/realistic.h"
#include "core/sampler.h"
#include "montecarlo.h"
#include "paramset.h"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;


RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)//specfile identify which lens system is used
	: Camera(cam2world, sopen, sclose, f), filmDistance(filmdistance), apertureDiameter(aperture_diameter), filmDiag(filmdiag)
	// pbrt-v2 doesnot specify hither and yon
{
   // YOUR CODE HERE -- build and store datastructures representing the given lens
   // and film placement.
	ifstream lenfile(specfile.c_str());
	string line;
	assert(lenfile.is_open());
	while (getline(lenfile, line)) {
		if (line[0] == '#') continue;
		stringstream ss(line);
		float radius, thick, n, aperture;
		ss >> radius >> thick >> n >> aperture;
		lenSystem.push_back(SingleLen(radius, thick, n, aperture));
	}
	float diagScale = filmDiag / sqrt(film->xResolution * film->xResolution + film->yResolution * film->yResolution);
	float xScale = film->xResolution * diagScale;
	float yScale = film->yResolution * diagScale;
	RasToScreen = Scale(-xScale, yScale, 1.f) * Translate(Vector(-0.5f, -0.5f, 0.f)) * Scale(1.f / film->xResolution, 1.f / film->yResolution, 1.f);
	//RasToScreen = Translate(Vector(film->xResolution * diagScale / 2, -film->yResolution * diagScale / 2, 0.f)) * Scale(diagScale, diagScale, 1.f) * Scale(-1.f, 1.f, 1.f);
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  // YOUR CODE HERE -- make that ray!
  
  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens
	Point Pras(sample.imageX, sample.imageY, 0.f);
	Point Pscreen;
	RasToScreen(Pras, &Pscreen);
	*ray = Ray(Pscreen, Vector(0, 0, 1), 0.f, INFINITY);
	ray->time = sample.time;
	float lensU, lensV;
	ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
	vector<SingleLen>::const_reverse_iterator rit = lenSystem.rbegin() ;
	lensU *= (*rit).aperture / 2;
	lensV *= (*rit).aperture / 2;
	Point Pdir(lensU, lensV, filmDistance);
	Ray intermediateRay(Pscreen, Normalize(Vector(Pdir - Pscreen)), 0.f);
	for (; rit != lenSystem.rend(); ++rit) {
		ToSphereCoordinate(intermediateRay, rit);
		if (rit->radius == 0) continue;
		// Compute quadratic sphere coefficients
		float A = intermediateRay.d.x*intermediateRay.d.x + intermediateRay.d.y*intermediateRay.d.y + intermediateRay.d.z*intermediateRay.d.z;
		float B = 2 * (intermediateRay.d.x*intermediateRay.o.x + intermediateRay.d.y*intermediateRay.o.y + intermediateRay.d.z*intermediateRay.o.z);
		float C = intermediateRay.o.x*intermediateRay.o.x + intermediateRay.o.y*intermediateRay.o.y +
			intermediateRay.o.z*intermediateRay.o.z - rit->radius*rit->radius;
		// Solve quadratic equation for _t_ values
		float t0, t1;
		if (!Quadratic(A, B, C, &t0, &t1)) {
			return 0.f;
		}
		// Compute intersection distance along ray
		if (t0 > intermediateRay.maxt || t1 < intermediateRay.mint) {
			return 0.f;
		}
		float thit = t0;
		if (t0 < intermediateRay.mint) {
			thit = t1;
			if (thit > intermediateRay.maxt) {
				return 0.f;
			}
		}
		if (intermediateRay(thit).z < 0 && rit->radius > 0)
			thit = t1;
		Point Phit(intermediateRay(thit));
		if (sqrt(Phit.x * Phit.x + Phit.y * Phit.y) * 2 > rit->aperture) {
			return 0.f;
		}
		Vector Vhit_normal(Normalize(Vector(Phit)));
		if (rit->radius > 0) {
			Vhit_normal = -Vhit_normal;
		}
		Vector Vincident(intermediateRay.d);
		vector<SingleLen>::const_reverse_iterator ritNext = ++rit;
		rit--;
		float n2 = 1.f;
		if (ritNext != lenSystem.rend()) {
			n2 = ritNext->n ? ritNext->n : 1;
		}
		float n = rit->n / n2;
		float c1 = -Dot(Vincident, Vhit_normal);
		float c2 = sqrt(1 - n * n * (1 - c1 * c1));
		if (isnan(c2)) {
			return 0.f;
		}
		Vector Vtransmit(n * Vincident + (n * c1 - c2) * Vhit_normal);
		intermediateRay.o = Phit;
		intermediateRay.d = Vtransmit;
	}
	rit--;
	Transform move = Translate(Vector(0.f, 0.f, -rit->radius));
	move(intermediateRay.o, &(ray->o));
	//ray->o = intermediateRay.o;
	ray->d = intermediateRay.d;
	CameraToWorld(*ray, ray);
	float cosTheta = Dot(Vector(0, 0, 1), Normalize(Vector(Pdir - Pscreen)));
	float lenArea = pow(lenSystem.back().aperture / 2.f, 2) * M_PI;
	return pow(cosTheta, 4.f) / pow(fabs(filmDistance), 2.f) * lenArea;
}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}
