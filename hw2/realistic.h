
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"

class SingleLen {
public:
	SingleLen(float _radius, float _thick, float _n, float _aperture)
		:radius(_radius), thick(_thick), n(_n), aperture(_aperture) {
	}
	float radius;
	float thick;
	float n;
	float aperture;
};

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
	void ToSphereCoordinate(Ray& ray, vector<SingleLen>::const_reverse_iterator rit) const {
		float distance = 0.f;
		if (rit == lenSystem.rbegin()) { 
			distance = filmDistance - rit->radius;
		}
		else {
			vector<SingleLen>::const_reverse_iterator ritPrevious = --rit;
			rit++;
			distance = rit->thick - rit->radius + ritPrevious->radius;
		}
		Transform move = Translate(Vector(0.f, 0.f, -distance));
		move(ray.o, &(ray.o));
	} 


	vector<SingleLen> lenSystem;
	float filmDistance;
	float apertureDiameter;
	float filmDiag;
	Transform RasToScreen;
private:
	// RealisticCamera Public Methods

};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H