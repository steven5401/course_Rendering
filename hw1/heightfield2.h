
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

#ifndef PBRT_SHAPES_HEIGHTFIELD2_H
#define PBRT_SHAPES_HEIGHTFIELD2_H

// shapes/heightfield2.h*
#include "shape.h"

// Heightfield2 Declarations
class Heightfield2 : public Shape {
public:
	// Heightfield2 Public Methods
	Heightfield2(const Transform *o2, const Transform *w2o, bool ro, int nu, int nv, const float *zs);
	~Heightfield2();
	bool CanIntersect() const;
	BBox ObjectBound() const;
	bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
		DifferentialGeometry *dg) const;
	bool IntersectP(const Ray &ray) const;
	void GetShadingGeometry(const Transform &obj2world,
		const DifferentialGeometry &dg,
		DifferentialGeometry *dgShading) const;
private:
	int posToVoxel(const Point &P, int axis) const {
		Assert(axis < 2);
		int axisnum = axis ? ny : nx;
		BBox bounds = this->ObjectBound();
		Vector invWidth(nx - 1, ny - 1, 1);
		int v = Float2Int((P[axis] - bounds.pMin[axis]) *
			invWidth[axis]);
		return Clamp(v, 0, axisnum - 2);
	}
	float voxelToPos(int p, int axis) const {
		Assert(axis < 2);
		BBox bounds = this->ObjectBound();
		return bounds.pMin[axis] + p * width[axis];
	}
	bool triangleInsert(const Ray &ray, float *tHit, float *rayEpsilon, 
		                DifferentialGeometry *dg, const Point &p1,
		                const Point &p2, const Point &p3) const;
	// Heightfield2 Private Data
	float *z;
	int nx, ny;
	float width[2];
	Normal* normalArray;
};


Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
	bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_HEIGHTFIELD2_H
