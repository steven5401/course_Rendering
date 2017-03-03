
/*
pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

This file is part of pbrt.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


// integrators/photonmap.cpp*
#include "stdafx.h"
#include "integrators/progressivephotonmap.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"
#include "film.h"
#include <string>


// PPhotonIntegrator Local Declarations
struct PPhoton {
	PPhoton(const Point &pp, const Spectrum &wt, const Vector &w, BSDF* const _bsdf)
		: p(pp), alpha(wt), wi(w) , bsdf(_bsdf) { }
	Point p;
	Spectrum alpha;
	Vector wi;
	BSDF* bsdf;
};

struct HitPoint {
	HitPoint(Point _p, Vector _n, float _x, float _y, float _alpha, Spectrum _flux, Vector _wo)
		: p(_p), n(_n), imageX(_x), imageY(_y), alpha(_alpha), flux(_flux), wo(_wo) {
		nPhotons = 0;
	}
	HitPoint() {}
	Point p;
	Vector n;
	float imageX, imageY;
	float radius2;
	float alpha;
	Spectrum flux;
	u_int nPhotons;
	Vector wo;
};



struct PPhotonProcess {
	// PPhotonProcess Public Methods
	PPhotonProcess() {};
	void operator()(const Point &p, const PPhoton &photon, float dist2,
		float &maxDistSquared);//can't change
	vector<PPhoton> pphotoninside;
};



inline float kernel(const PPhoton *photon, const Point &p, float maxDist2);


inline void PPhotonProcess::operator()(const Point &p,
	const PPhoton &photon, float distSquared, float &maxDistSquared) {
	pphotoninside.push_back(photon);
}


inline float kernel(const PPhoton *photon, const Point &p,
	float maxDist2) {
	float s = (1.f - DistanceSquared(photon->p, p) / maxDist2);
	return 3.f * INV_PI * s * s;
}


static Spectrum LPPhoton(vector<PPhoton> _pphotons, Vector wo, RNG rng) {
	Spectrum L(0.);
	BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
		BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
	for (int i = 0; i < _pphotons.size(); i++) {
		if (_pphotons[i].bsdf->NumComponents(nonSpecular) == 0)
			continue;
		// Accumulate light from nearby photons
		// Estimate reflected light from photons
		Normal Nf = Dot(_pphotons[i].wi, _pphotons[i].bsdf->dgShading.nn) < 0 ? -_pphotons[i].bsdf->dgShading.nn : _pphotons[i].bsdf->dgShading.nn;
		if (_pphotons[i].bsdf->NumComponents(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY)) > 0) {
			// Compute exitant radiance from photons for glossy surface
			BxDFType flag = Dot(Nf, _pphotons[i].wi) > 0.f ? BSDF_ALL_REFLECTION : BSDF_ALL_TRANSMISSION;
			L += _pphotons[i].bsdf->f(wo, _pphotons[i].wi, flag) * (_pphotons[i].alpha * INV_PI);
		}
		else {
			// Compute exitant radiance from photons for diffuse surface
			Spectrum Lr(0.), Lt(0.);
			if (Dot(Nf, _pphotons[i].wi) > 0.f)
				Lr += _pphotons[i].alpha;
			else
				Lt += _pphotons[i].alpha;
			L += (INV_PI * INV_PI) * (Lr * _pphotons[i].bsdf->rho(wo, rng, BSDF_ALL_REFLECTION) + Lt * _pphotons[i].bsdf->rho(wo, rng, BSDF_ALL_TRANSMISSION));
		}
	}
	return L;
}




// PPhotonIntegrator Method Definitions
PPhotonIntegrator::PPhotonIntegrator(int ncaus, int nind,
	int nl, int mdepth, int mphodepth, float mdist, bool fg,
	int gs, float ga, u_int nphotons, int npass, float rratio, int maxdepth) {
	nCausticPPhotonsWanted = ncaus;
	nIndirectPPhotonsWanted = nind;
	nLookup = nl;
	maxSpecularDepth = mdepth;
	maxPPhotonDepth = mphodepth;
	maxDistSquared = mdist * mdist;
	finalGather = fg;
	cosGatherAngle = cos(Radians(ga));
	gatherSamples = gs;
	nCausticPaths = nIndirectPaths = 0;
	lightSampleOffsets = NULL;
	bsdfSampleOffsets = NULL;
	//euphoria5401//
	nPhotonsPerPass = nphotons;
	totalPhotons = 0;
	nPass = npass;
	reductionRatio = rratio;
	maxDepth = maxdepth;
	//euphoria5401//
}


PPhotonIntegrator::~PPhotonIntegrator() {
	delete[] lightSampleOffsets;
	delete[] bsdfSampleOffsets;
}


void PPhotonIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
	const Scene *scene) {
	// Allocate and request samples for sampling all lights
	uint32_t nLights = scene->lights.size();
	lightSampleOffsets = new LightSampleOffsets[nLights];
	bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
	for (uint32_t i = 0; i < nLights; ++i) {
		const Light *light = scene->lights[i];
		int nSamples = light->nSamples;
		if (sampler) nSamples = sampler->RoundSize(nSamples);
		lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
		bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
	}
}


void PPhotonIntegrator::Preprocess(const Scene *scene,
	const Camera *camera, const Renderer *renderer) {
	//do nothing
	//euphoria5401//
}



Spectrum PPhotonIntegrator::Li(const Scene *scene, const Renderer *renderer,
	const RayDifferential &ray, Intersection &isect,
	const Sample *sample, RNG &rng, MemoryArena &arena) const {
	Spectrum L(0.);
	RayDifferential hpray(ray);
	Vector wo = -hpray.d;
	L += isect.Le(wo);
	// Follow photon path through scene and record intersections  {{{
	bool specularPath = false;
	bool stop = false;
	int nIntersections = 0;
	bool hasIntersect = true;
	float alpha = 1.;
	while (hasIntersect) {
		++nIntersections;
		// Handle photon/surface intersection  {{{
		Vector wo = -hpray.d;
		BSDF *bsdf = isect.GetBSDF(hpray, arena);
		const Point &p = bsdf->dgShading.p;
		const Normal &n = bsdf->dgShading.nn;
		BxDFType specularType = BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR);
		bool hasNonSpecular = (bsdf->NumComponents() > bsdf->NumComponents(specularType));
		if (hasNonSpecular || stop) {
			// Record the information into a HitPoint
			// Compute emitted light if ray hit an area light source
			L += isect.Le(wo);
			// Compute direct lighting for photon map integrator
			L += UniformSampleAllLights(scene, renderer, arena, p, n,
				wo, isect.rayEpsilon, ray.time, bsdf, sample, rng,
				lightSampleOffsets, bsdfSampleOffsets);

			HitPoint *hp = new HitPoint(p, Vector(n), sample->imageX, sample->imageY, alpha, L, wo);
			hitpoints.push_back(hp);
			break;
		}
		// }}}
		// Sample new photon ray direction  {{{
		Vector wi;
		float pdf;
		BxDFType flags;
		// Get random numbers for sampling outgoing photon direction  {{{
		// }}}
		Spectrum fr = bsdf->Sample_f(wo, &wi, BSDFSample(rng), &pdf, BSDF_ALL, &flags);
		if (fr.IsBlack() || pdf == 0.f)
			break;
		specularPath = (nIntersections == 1 || specularPath) &&
			((flags & BSDF_SPECULAR) != 0);
		hpray = RayDifferential(isect.dg.p, wi, hpray,
			isect.rayEpsilon);
		// }}}
		// Possibly terminate ray path  {{{
		if (nIntersections > 3) {
			float continueProbability = .5f;
			if (rng.RandomFloat() > continueProbability)
				stop = true;
			alpha /= continueProbability;
		}
		// }}}
		hasIntersect = scene->Intersect(hpray, &isect);
	}

	return L;
}

void PPhotonIntegrator::Postprocess(const Scene *scene, const Renderer* renderer, const Camera* camera) {
	Point min(hitpoints[0]->p.x, hitpoints[0]->p.y, hitpoints[0]->p.z);
	Point max(hitpoints[0]->p.x, hitpoints[0]->p.y, hitpoints[0]->p.z);
	for (int i = 0; i < hitpoints.size(); i++) {
		Point p = hitpoints[i]->p;
		if (p.x<min.x) min.x = p.x;
		if (p.y<min.y) min.y = p.y;
		if (p.z<min.z) min.z = p.z;
		max.x = std::max(p.x, max.x);
		max.y = std::max(p.y, max.y);
		max.z = std::max(p.z, max.z);
	}
	Vector v = max - min;
	const Film* f = camera->film;
	const int x = f->xResolution;
	const int y = f->yResolution;
	float r = ((v.x + v.y + v.z) / 3.0) / ((x + y) / 2.0) * 2.0;
	for (int i = 0; i < hitpoints.size(); i++) {
		hitpoints[i]->radius2 = r * r;
	}

	PhotonTracing(scene, renderer, camera);
	RadianceEvaluation(scene, camera);
	//fprintf(stderr, "after PhotonTracing()\n");
	//fprintf(stderr, "after RadianceEvaluation()\n");
}

void PPhotonIntegrator::PhotonTracing(const Scene *scene, const Renderer* renderer, const Camera* camera) {
	ProgressReporter progress(nPass, "Photon tracing"); // NOBOOK
	for (int i = 0; i < nPass; i++) {
		//fprintf(stderr, "\nnow in pass %d\r", i);
		PhotonTracingPass(scene, renderer, camera, i);
		progress.Update();
		
	}
	progress.Done();
}

void PPhotonIntegrator::PhotonTracingPass(const Scene *scene, const Renderer* renderer, const Camera* camera, int pass) {
	if (scene->lights.size() == 0) return;
	//ProgressReporter progress(nPhotonsPerPass, "Shooting photons"); // NOBOOK
	// Initialize photon shooting statistics
	MemoryArena arena;
	RNG rng(0);
	vector<PPhoton> pphotons;
	for (u_int nshot = 0; nshot < nPhotonsPerPass; ++nshot) {
		// Trace a photon path and store contribution
		// Choose 4D sample values for photon
		float u[5];
		u[0] = (float)RadicalInverse(nshot + pass + 1, 2);
		u[1] = (float)RadicalInverse(nshot + pass + 1, 3);
		u[2] = (float)RadicalInverse(nshot + pass + 1, 5);
		u[3] = (float)RadicalInverse(nshot + pass + 1, 7);
		u[4] = (float)RadicalInverse(nshot + pass + 1, 11);
		// Choose light to shoot photon from
		int nLights = int(scene->lights.size());
		int lightNum = min(Floor2Int(nLights * (float)RadicalInverse(nshot + pass + 1, 13)), nLights - 1);
		Light *light = scene->lights[lightNum];
		float lightPdf = 1.f / nLights;
		// Generate _photonRay_ from light source and initialize _alpha_
		RayDifferential photonRay;
		float pdf;
		Normal Nl;
		LightSample ls(u[0], u[1], u[2]);
		float time = camera ? camera->shutterOpen : 0.f;
		Spectrum alpha = light->Sample_L(scene, ls, u[3], u[4], time, &photonRay, &Nl, &pdf);
		if (pdf == 0.f || alpha.IsBlack()) continue;
		alpha /= pdf * lightPdf;
		if (!alpha.IsBlack()) {
			// Follow photon path through scene and record intersections  {{{
			bool specularPath = false;
			Intersection photonIsect;
			int nIntersections = 0;
			while (scene->Intersect(photonRay, &photonIsect)) {
				++nIntersections;
				// Handle photon/surface intersection  {{{
				alpha *= renderer->Transmittance(scene, photonRay, NULL, rng, arena);
				Vector wo = -photonRay.d;
				BSDF *photonBSDF = photonIsect.GetBSDF(photonRay, arena);
				BxDFType specularType = BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR);
				bool hasNonSpecular = (photonBSDF->NumComponents() > photonBSDF->NumComponents(specularType));
				if (hasNonSpecular && nIntersections > 1) {
					// Deposit photon at surface
					PPhoton photon(photonIsect.dg.p, alpha, wo, photonBSDF);
					pphotons.push_back(photon);
					//progress.Update(); // NOBOOK
				}
				// }}}
				// Sample new photon ray direction  {{{
				Vector wi;
				float pdf;
				BxDFType flags;
				Spectrum fr = photonBSDF->Sample_f(wo, &wi, BSDFSample(rng),
					&pdf, BSDF_ALL, &flags);
				if (fr.IsBlack() || pdf == 0.f)
					break;
				specularPath = (nIntersections == 1 || specularPath) &&
					((flags & BSDF_SPECULAR) != 0);
				alpha *= fr * AbsDot(wi, photonBSDF->dgShading.nn) / pdf;
				photonRay = RayDifferential(photonIsect.dg.p, wi, photonRay,
					photonIsect.rayEpsilon);
				// }}}
				// Possibly terminate photon path  {{{
				if (nIntersections > 3) {
					float continueProbability = .5f;
					if (rng.RandomFloat() > continueProbability)
						break;
					alpha /= continueProbability;
				}
				// }}}
			}
			// }}}
		}
	}
	if (pphotons.size() == 0) return;
	KdTree<PPhoton>* pphotonmap = new KdTree<PPhoton>(pphotons);
	for (int i = 0; i < hitpoints.size(); i++) {
		PPhotonProcess proc;
		pphotonmap->Lookup(hitpoints[i]->p, proc, hitpoints[i]->radius2);
		if (hitpoints[i]->nPhotons + proc.pphotoninside.size() == 0) continue;
		float factor = (hitpoints[i]->nPhotons + reductionRatio * proc.pphotoninside.size()) / (hitpoints[i]->nPhotons + proc.pphotoninside.size());
		hitpoints[i]->nPhotons += reductionRatio * proc.pphotoninside.size();
		hitpoints[i]->radius2 *= factor;
		hitpoints[i]->flux += LPPhoton(proc.pphotoninside, hitpoints[i]->wo, rng);
		hitpoints[i]->flux *= factor;
		arena.FreeAll();
	}
	totalPhotons += nPhotonsPerPass;
	//progress.Done(); // NOBOOK
}

void PPhotonIntegrator::RadianceEvaluation(const Scene *scene, const Camera* camera) {
	Ray ray;
	Sample *sample = new Sample();
	ProgressReporter progress(hitpoints.size(), "Radiance evaluation"); // NOBOOK
	for (int i = 0; i < hitpoints.size(); i++) {
		HitPoint* hp = hitpoints[i];
		hp->flux *= INV_PI / (hp->radius2 * totalPhotons);

		if (hp->flux.HasNaNs()) {
			Error("Not-a-number radiance value returned "
				"for image sample.  Setting to black.");
			hp->flux = Spectrum(0.f);
		}
		else if (hp->flux.y() < -1e-5) {
			Error("Negative luminance value, %g, returned "
				"for image sample.  Setting to black.", hp->flux.y());
			hp->flux = Spectrum(0.f);
		}
		else if (isinf(hp->flux.y())) {
			Error("Infinite luminance value returned "
				"for image sample.  Setting to black.");
			hp->flux = Spectrum(0.f);
		}

		sample->imageX = hp->imageX;
		sample->imageY = hp->imageY;
		camera->film->AddSample(*sample, hp->flux);
		progress.Update(); // NOBOOK
	}
	progress.Done(); // NOBOOK
	delete sample;
}

PPhotonIntegrator *CreatePPhotonMapSurfaceIntegrator(const ParamSet &params) {
	int nCaustic = params.FindOneInt("causticphotons", 20000);
	int nIndirect = params.FindOneInt("indirectphotons", 100000);
	int nUsed = params.FindOneInt("nused", 50);
	if (PbrtOptions.quickRender) nCaustic = nCaustic / 10;
	if (PbrtOptions.quickRender) nIndirect = nIndirect / 10;
	if (PbrtOptions.quickRender) nUsed = max(1, nUsed / 10);
	int maxSpecularDepth = params.FindOneInt("maxspeculardepth", 5);
	int maxPPhotonDepth = params.FindOneInt("maxphotondepth", 5);
	bool finalGather = params.FindOneBool("finalgather", true);
	int gatherSamples = params.FindOneInt("finalgathersamples", 32);
	if (PbrtOptions.quickRender) gatherSamples = max(1, gatherSamples / 4);
	float maxDist = params.FindOneFloat("maxdist", .1f);
	float gatherAngle = params.FindOneFloat("gatherangle", 10.f);

	int maxDepth = params.FindOneInt("maxdepth", 5);
	int nPhotonsPerPass = params.FindOneInt("nphotonsperpass", 100000);
	int nPass = params.FindOneInt("npass", 64);
	float reductionRation = params.FindOneFloat("alpha", .7f);
	return new PPhotonIntegrator(nCaustic, nIndirect,
		nUsed, maxSpecularDepth, maxPPhotonDepth, maxDist, finalGather, gatherSamples,
		gatherAngle, nPhotonsPerPass, nPass, reductionRation, maxDepth);
}


