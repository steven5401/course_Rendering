
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PPHOTONMAP_H
#define PBRT_INTEGRATORS_PPHOTONMAP_H

// integrators/photonmap.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"

struct PPhoton;
struct RadiancePPhoton;
struct ClosePPhoton;
struct HitPointProcess;
struct RadiancePPhotonProcess;
struct HitPoint;


// PPhotonIntegrator Declarations
class PPhotonIntegrator : public Integrator {
public:
	// PPhotonIntegrator Public Methods
	PPhotonIntegrator(int ncaus, int nindir, int nLookup, int maxspecdepth,
		int maxphotondepth, float maxdist, bool finalGather, int gatherSamples,
		float ga,  u_int nphotons, int npass, float rratio, int maxdepth);
	~PPhotonIntegrator();
	Spectrum Li(const Scene *scene, const Renderer *renderer,
		const RayDifferential &ray, Intersection &isect, const Sample *sample,
		RNG &rng, MemoryArena &arena) const;
	void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
	void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);
	//euphoria5401//
	void Postprocess(const Scene* scene, const Renderer* renderer, const Camera* camera);
	void PhotonTracing(const Scene *scene, const Renderer* renderer, const Camera* camera);
	void PhotonTracingPass(const Scene *scene, const Renderer* renderer, const Camera* camera, int pass);
	void RadianceEvaluation(const Scene *scene, const Camera* camera);
	//euphoria5401//
private:
	// PPhotonIntegrator Private Methods
	friend class PPhotonShootingTask;

	// PPhotonIntegrator Private Data
	uint32_t nCausticPPhotonsWanted, nIndirectPPhotonsWanted, nLookup;
	float maxDistSquared;
	int maxSpecularDepth, maxPPhotonDepth;
	bool finalGather;
	int gatherSamples;
	float cosGatherAngle;

	// Declare sample parameters for light source sampling
	LightSampleOffsets *lightSampleOffsets;
	BSDFSampleOffsets *bsdfSampleOffsets;
	BSDFSampleOffsets bsdfGatherSampleOffsets, indirGatherSampleOffsets;
	int nCausticPaths, nIndirectPaths;
	//euphoria5401//
	u_int nPhotonsPerPass, totalPhotons;
	int nPass;
	float reductionRatio;
	int maxDepth;
	mutable vector<HitPoint*> hitpoints;
	//euphoria5401//
};


PPhotonIntegrator *CreatePPhotonMapSurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_PHOTONMAP_H
