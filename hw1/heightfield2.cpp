
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


// shapes/heightfield2.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
	bool ro, int x, int y, const float *zs)
	: Shape(o2w, w2o, ro) {
	nx = x;
	ny = y;
	z = new float[nx*ny];
	normalArray = new Normal[nx*ny];
	width[0] = 1.0 / (nx - 1);
	width[1] = 1.0 / (ny - 1);
	memcpy(z, zs, nx*ny*sizeof(float));
	for (int i = 0; i < nx - 1; i++) {
		for (int j = 0; j < ny - 1; j++) {
			//left-top traingle
			float leftTop = z[(j + 1) * nx + i];
			float leftBot = z[j * nx + i];
			float rightBot = z[j * nx + i + 1];
			float rightTop = z[(j + 1) * nx + i + 1];
			Vector horizontal(width[0], 0, rightTop - leftTop);
			Vector vertical(0, -width[1], leftBot - leftTop);
			Normal normal = Normal(Normalize(Cross(horizontal, vertical)));
			normalArray[j * nx + i] += normal;
			normalArray[(j + 1) * nx + i] += normal;
			normalArray[(j + 1) * nx + i + 1] += normal;
			//right-bot triangle
			Vector horizontal2(-width[0], 0, leftBot - rightBot);
			Vector vertical2(0, width[1], rightTop - rightBot);
			Normal normal2 = Normal(Normalize(Cross(horizontal2, vertical2)));
			normalArray[(j + 1) * nx + i + 1] += normal2;
			normalArray[j * nx + i + 1] += normal2;
			normalArray[j * nx + i] += normal2;
		}
	}
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			normalArray[j * nx + i] = Normalize((*o2w)((Normalize(normalArray[j * nx + i]))));
		}
	}
}


Heightfield2::~Heightfield2() {
	delete[] z;
}


BBox Heightfield2::ObjectBound() const {
	float minz = z[0], maxz = z[0];
	for (int i = 1; i < nx*ny; ++i) {
		if (z[i] < minz) minz = z[i];
		if (z[i] > maxz) maxz = z[i];
	}
	return BBox(Point(0, 0, minz), Point(1, 1, maxz));
}


bool Heightfield2::CanIntersect() const {
	return true;
}


bool Heightfield2::triangleInsert(const Ray &ray, float *tHit, float *rayEpsilon,
	DifferentialGeometry *dg, const Point &p1, const Point &p2, const Point &p3) const {
	//p1 is 90 degree point
	//p2 is horizonal point
	//p3 is vertical point
	Vector e1 = p2 - p1;
	Vector e2 = p3 - p1;
	Vector s1 = Cross(ray.d, e2);
	float divisor = Dot(s1, e1);

	if (divisor == 0.)
		return false;
	float invDivisor = 1.f / divisor;

	// Compute first barycentric coordinate
	Vector s = ray.o - p1;
	float b1 = Dot(s, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
		return false;

	// Compute second barycentric coordinate
	Vector s2 = Cross(s, e1);
	float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b1 + b2 > 1.)
		return false;

	// Compute _t_ to intersection point
	float t = Dot(e2, s2) * invDivisor;
	if (t < ray.mint || t > ray.maxt)
		return false;
	
	// Compute triangle partial derivatives
	Vector dpdu, dpdv;
	dpdu = p2.x > p1.x ? (p2 - p1) : (p1 - p2);
	dpdu /= dpdu.x;
	dpdv = p3.y > p1.y ? (p3 - p1) : (p1 - p3);
	dpdv /= dpdv.y;

	// Fill in _DifferentialGeometry_ from triangle hit
	const Transform &o2w = *ObjectToWorld;
	Point o2w_ray_t = o2w(ray(t));
	*dg = DifferentialGeometry(o2w_ray_t, o2w(dpdu), o2w(dpdv),
		o2w(Normal(0, 0, 0)), o2w(Normal(0, 0, 0)),
		ray(t).x, ray(t).y, this);
	*tHit = t;
	*rayEpsilon = 1e-3f * *tHit;
	//PBRT_RAY_TRIANGLE_INTERSECTION_HIT(const_cast<Ray *>(&ray), t);
	return true;
}

bool Heightfield2::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
	DifferentialGeometry *dg) const {
	float rayT;
	Ray ray;
	(*WorldToObject)(r, &ray);
	BBox bounds = this->ObjectBound();
	if (bounds.Inside(ray(ray.mint))) {
		rayT = ray.mint;
	} else if (!bounds.IntersectP(ray, &rayT)) {
		//PBRT_GRID_RAY_MISSED_BOUNDS();
		return false;
	}
	Point gridIntersect = ray(rayT);
	// Set up 2D DDA for ray
	float NextCrossingT[2], DeltaT[2];
	int Step[2], Out[2], Pos[2];
	for (int axis = 0; axis < 2; ++axis) {
		// Compute current voxel for axis
		Pos[axis] = posToVoxel(gridIntersect, axis);
		float invRayD = 1 / ray.d[axis];
		if (ray.d[axis] >= 0) {
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT +
				(voxelToPos(Pos[axis] + 1, axis) - gridIntersect[axis]) * invRayD;
			DeltaT[axis] = width[axis] * invRayD;
			Step[axis] = 1;
			Out[axis] = axis ? ny - 1: nx - 1;
		}
		else {
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT +
				(voxelToPos(Pos[axis], axis) - gridIntersect[axis]) * invRayD;
			DeltaT[axis] = -width[axis] * invRayD;
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}
	bool hitSomething = false;
	bool first = true;
	for (;;) {
		// Check for intersection in current voxel and advance to next
		/*
		p4----p3
		| 1 /  |
		|  /   |
		| / 0  |
		p1----p2

								  y+        point @ hf4x4.pbrt mapping to coordinate@hftest.exr
			-------x+      		 /
			|					/
			|			==>	   /
			|				   | \
			y+				   |  \
							   |   \
							   z+   x+
		*/

		Point p1(Pos[0] * width[0], Pos[1] * width[1], z[Pos[1] * nx + Pos[0]]);
		Point p2((Pos[0] + 1) * width[0], Pos[1] * width[1], z[Pos[1] * nx + Pos[0] + 1]);
		Point p3((Pos[0] + 1) * width[0], (Pos[1] + 1) * width[1], z[(Pos[1] + 1) * nx + Pos[0] + 1]);
		Point p4(Pos[0] * width[0], (Pos[1] + 1) * width[1], z[(Pos[1] + 1) * nx + Pos[0]]);
		bool stepAxis;
		if (first) {
			first = false;
			float x = ray.d.x, y = ray.d.y;
			float m = y / x;
			bool firstmeet = (m < 1 && x > 0 && y > 0) ||
				(m < 0 && x > 0 && y < 0) ||
				(m > 1 && x < 0 && y < 0);
			if (firstmeet) {
				if (triangleInsert(ray, tHit, rayEpsilon, dg, p4, p3, p1)) {
					return true;
				}
				else if (triangleInsert(ray, tHit, rayEpsilon, dg, p2, p1, p3)) {
					return true;
				}
			}
			else {
				if (triangleInsert(ray, tHit, rayEpsilon, dg, p2, p1, p3)) {
					return true;
				}
				else if (triangleInsert(ray, tHit, rayEpsilon, dg, p4, p3, p1)) {
					return true;
				}
			}
		}
		else {
			/*bool firstmeet = (m < 1 && x > 0 && y > 0) ||
							 (m < 0 && x > 0 && y < 0) ||
							 (m > 1 && x < 0 && y < 0);*/
			if (/*firstmeet*/!stepAxis) {
				if (triangleInsert(ray, tHit, rayEpsilon, dg, p4, p3, p1)) {
					return true;
				}
				else if (triangleInsert(ray, tHit, rayEpsilon, dg, p2, p1, p3)) {
					return true;
				}
			}
			else {
				if (triangleInsert(ray, tHit, rayEpsilon, dg, p2, p1, p3)) {
					return true;
				}
				else if (triangleInsert(ray, tHit, rayEpsilon, dg, p4, p3, p1)) {
					return true;
				}
			}
		}
		// Advance to next voxel

		// Find _stepAxis_ for stepping to next voxel
		stepAxis = (NextCrossingT[0] > NextCrossingT[1]);
		if (ray.maxt < NextCrossingT[stepAxis])
			break;
		Pos[stepAxis] += Step[stepAxis];
		if (Pos[stepAxis] == Out[stepAxis])
			break;
		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}
	return hitSomething;
}
bool Heightfield2::IntersectP(const Ray &ray) const {
	DifferentialGeometry dg;
	float tHit;
	float rayEpsilon;
	return Intersect(ray, &tHit, &rayEpsilon, &dg);
}

void Heightfield2::GetShadingGeometry(const Transform &obj2world,
	const DifferentialGeometry &dg,
	DifferentialGeometry *dgShading) const {

	// Compute barycentric coordinates for point
	float b[3];

	// Initialize _A_ and _C_ matrices for barycentrics
	float uv[3][2];
	const Transform* w2o = dg.shape->WorldToObject;
	Point hitPoint = (*w2o)(dg.p);//local
	int pos[2];
	pos[0] = posToVoxel(hitPoint, 0);
	pos[1] = posToVoxel(hitPoint, 1);
	Point leftBottomPoint(voxelToPos(pos[0], 0), voxelToPos(pos[1], 1), 0);//local
	Vector m = hitPoint - leftBottomPoint;
	bool slopeBigOne = m.y > m.x;
	Normal leftBot = normalArray[pos[1] * nx + pos[0]];
	Normal leftTop = normalArray[(pos[1] + 1) * nx + pos[0]];
	Normal rightTop = normalArray[(pos[1] + 1) * nx + pos[0] + 1];
	Normal rightBot = normalArray[pos[1] * nx + pos[0] + 1];//world
	if (slopeBigOne) {
		uv[0][0] = voxelToPos(pos[0], 0);
		uv[0][1] = voxelToPos(pos[1], 1);
		uv[1][0] = voxelToPos(pos[0] + 1, 0);
		uv[1][1] = voxelToPos(pos[1] + 1, 1);
		uv[2][0] = voxelToPos(pos[0], 0);
		uv[2][1] = voxelToPos(pos[1] + 1, 1);
	} else {
		uv[0][0] = voxelToPos(pos[0], 0);
		uv[0][1] = voxelToPos(pos[1], 1);
		uv[1][0] = voxelToPos(pos[0] + 1, 0);
		uv[1][1] = voxelToPos(pos[1], 1);
		uv[2][0] = voxelToPos(pos[0] + 1, 0);
		uv[2][1] = voxelToPos(pos[1] + 1, 1);//local
	}

	float A[2][2] =
	{ { uv[1][0] - uv[0][0], uv[2][0] - uv[0][0] },
	{ uv[1][1] - uv[0][1], uv[2][1] - uv[0][1] } };
	float C[2] = { dg.u - uv[0][0], dg.v - uv[0][1] };
	if (!SolveLinearSystem2x2(A, C, &b[1], &b[2])) {
		// Handle degenerate parametric mapping
		b[0] = b[1] = b[2] = 1.f / 3.f;
	}
	else
		b[0] = 1.f - b[1] - b[2];

	// Use _n_ and _s_ to compute shading tangents for triangle, _ss_ and _ts_
	Normal ns;//world
	Vector ss, ts;//world
	if (slopeBigOne) {
		/*
		2---1
		|  /
		| /
		|/
		0
		*/
		ns = Normalize(b[0] * leftBot +
			b[1] * rightTop + 
			b[2] * leftTop);
	} else {
		/*
		    2
		   /|
		  / |
		 /  |
		0---1
		*/
		ns = Normalize(b[0] * leftBot +
			b[1] * rightBot +
			b[2] * rightTop);//worldCoordinate
	}
	ss = Normalize(dg.dpdu);//worldCoordinate

	ts = Cross(ss, ns);//worldCoordinate
	if (ts.LengthSquared() > 0.f) {
		ts = Normalize(ts);
		ss = Cross(ts, ns);
	}
	else
		CoordinateSystem((Vector)ns, &ss, &ts);
	Normal dndu, dndv;
	/*
	// Compute $\dndu$ and $\dndv$ for triangle shading geometry
		// Compute deltas for triangle partial derivatives of normal
		float du1 = uv[0][0] - uv[2][0];
		float du2 = uv[1][0] - uv[2][0];
		float dv1 = uv[0][1] - uv[2][1];
		float dv2 = uv[1][1] - uv[2][1];
		Normal dn1, dn2;
		if (slopeBigOne) {
			dn1 = leftBot - leftTop;
			dn2 = rightTop - leftTop;
		}
		else {
			dn1 = leftBot - rightTop;
			dn2 = rightBot - rightTop;
		}
		float determinant = du1 * dv2 - dv1 * du2;
		if (determinant == 0.f)
			dndu = dndv = Normal(0, 0, 0);
		else {
			float invdet = 1.f / determinant;
			dndu = (dv2 * dn1 - dv1 * dn2) * invdet;
			dndv = (-du2 * dn1 + du1 * dn2) * invdet;
		}*/
		//dndu /= dndu.x;
		//dndv /= dndv.y;
	//	dndu = dndv = Normal(0, 0, 0);/
	if (slopeBigOne) {
		dndu = (leftTop - rightTop) / width[0];
		dndv = (leftTop - leftBot) / width[1];
	}
	else {
		dndu = (rightBot - leftBot) / width[0];
		dndv = (rightTop - rightBot) / width[1];
	}
	*dgShading = DifferentialGeometry(dg.p, ss, ts,
		obj2world(dndu), obj2world(dndv),
		dg.u, dg.v, dg.shape);
	dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx;
	dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy;
	dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy;
//*dgShading = dg;
}

Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
	bool reverseOrientation, const ParamSet &params) {
	int nu = params.FindOneInt("nu", -1);
	int nv = params.FindOneInt("nv", -1);
	int nitems;
	const float *Pz = params.FindFloat("Pz", &nitems);
	Assert(nitems == nu*nv);
	Assert(nu != -1 && nv != -1 && Pz != NULL);
	return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}


