
#pragma once

#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/vertex/component_ocf.h>
#include <vcg/simplex/edge/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/simplex/face/component_ocf.h>

#include <vcg/complex/used_types.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/allocate.h>

#include <vcg/simplex/face/topology.h>

#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/color.h>
#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/update/position.h>
#include <vcg/complex/trimesh/update/quality.h>
#include <vcg/complex/trimesh/update/selection.h>
#include <vcg/complex/trimesh/update/topology.h>

#include <vcg\space\point3.h>

#include <cstdlib> //for rand()
#include <ctime> //for time()

#include <vector>
using std::vector;
using namespace vcg;
//用vcg库定义三维网格结构

class CVertex;
class CFace;

//class CUsedTypes: public vcg::UsedTypes< vcg::Use<CVertex>::AsVertexType, vcg::Use<CFace>::AsFaceType>{};
class CEdge;
class CUsedTypes : public vcg::UsedTypes < vcg::Use<CVertex>::AsVertexType, vcg::Use<CEdge>::AsEdgeType, vcg::Use<CFace>::AsFaceType > {};

class CVertex : public vcg::Vertex<CUsedTypes, vcg::vertex::InfoOcf, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b, vcg::vertex::BitFlags,
	vcg::vertex::VFAdjOcf,          /*  0b */
	vcg::vertex::MarkOcf,           /*  0b */
	vcg::vertex::TexCoordfOcf,      /*  0b */
	vcg::vertex::CurvaturefOcf,     /*  0b */
	vcg::vertex::CurvatureDirfOcf,  /*  0b */
	vcg::vertex::RadiusfOcf         /*  0b */ >
{
public:
	vector<int> neighbors;
	vector<int> original_neighbors;

	bool bIsOriginal;
	int m_index;

  int dual_index;
	int target_index;

	bool is_fixed_sample; //feature points (blue color) 
	bool is_skel_ignore;

	/* for skeletonization */
	double eigen_confidence;
	Point3f eigen_vector0; //Associate with the biggest eigen value
	Point3f eigen_vector1; // Also use for remember last better virtual point
	Point3f eigen_vector2; //The smallest eigen value : should be PCA normal N()
	
	double eigen_value0;
	double eigen_value1;
	double eigen_value2;



  bool is_skel_virtual; //in our papaer, we said bridge point instead of virtual point
	bool is_skel_branch;
  bool is_fixed_original; 
  
  double skel_radius; // remember radius for branches
	float nearest_neighbor_dist;

  bool is_dual_sample;
  bool is_boundary;
	
public:
	operator Point3f &()
	{
		return P();
	}

	operator const Point3f &() const
	{
		return cP();
	}

	float & operator[](unsigned int i)
	{
		return P()[i];
	}

	CVertex() :
		m_index(0),
		dual_index(0),
		bIsOriginal(false),
		is_fixed_sample(false),
		eigen_confidence(0.),
		is_skel_branch(false),
		is_skel_ignore(false),
		is_skel_virtual(false),
		is_fixed_original(false),
		is_dual_sample(false),
		is_boundary(false),
		nearest_neighbor_dist(0.0),
		eigen_vector0(Point3f(1, 0, 0)),
		eigen_vector1(Point3f(0, 1, 0)),
		eigen_vector2(Point3f(0, 0, 1)),
		target_index(-1),
		eigen_value0(0.),
		eigen_value1(0.),
		eigen_value2(0.),
    skel_radius(-1.0)
		{
			N() = Point3f(0,0,0);

		}

	/* for skeletonization */
	void remove() //important, some time we don't want to earse points, just remove them
	{
		neighbors.clear();
		original_neighbors.clear();
		is_skel_ignore = true;
    N() = P();
		P() = Point3f(88888888888.8, 88888888888.8, 88888888888.8);
	}

	bool isSample_Moving()
	{
		return (!is_skel_ignore && !is_fixed_sample && !is_skel_branch);
	}

	bool isSample_JustMoving()
	{
		return (!is_skel_ignore && !is_fixed_sample && !is_skel_virtual && !is_skel_branch);
	}

	bool isSample_MovingAndVirtual()
	{
		return (!is_skel_ignore && !is_fixed_sample && is_skel_virtual && !is_skel_branch);
	}

	bool isSample_JustFixed()
	{
		return (!is_skel_ignore && is_fixed_sample && !is_skel_virtual && !is_skel_branch);
	}

	bool isSample_FixedAndBranched()
	{
		return (!is_skel_ignore && is_fixed_sample && !is_skel_virtual && is_skel_branch);
	}


	void setSample_JustMoving()
	{
		is_fixed_sample = false;
		is_skel_virtual = false;
		is_skel_branch  = false;
	}

	void setSample_MovingAndVirtual()
	{
		is_fixed_sample = false;
		is_skel_virtual = true;
		is_skel_branch  = false;
	}

	void setSample_JustFixed()
	{
		is_fixed_sample = true;
		is_skel_virtual = false;
		is_skel_branch  = false;
	}

	void setSample_FixedAndBranched()
	{
		is_fixed_sample = true;
		is_skel_virtual = false;
		is_skel_branch  = true;
	}

	void recompute_m_render()
	{
		srand(time(NULL));
		int x = rand()%1000;
		int y = rand()%1000;
		int z = rand()%1000;

		Point3f normal = N();
		normal.Normalize();

		Point3f helper(x/1000.0, y/1000.0, z/1000.0);
		Point3f new_m3_to_m5 = normal ^ helper;
		new_m3_to_m5.Normalize();
		Point3f new_m6_to_m8 = normal ^ new_m3_to_m5;

		eigen_vector0 = new_m3_to_m5;
		eigen_vector1 = new_m6_to_m8;

	}
};

// class CFace : public vcg::Face<CUsedTypes, vcg::face::VertexRef> {};
// class CMesh : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace> > {};

class CEdge : public vcg::Edge < CUsedTypes, vcg::edge::EVAdj >
{
public:
	inline CEdge(){};
	inline CEdge(CVertex * v0, CVertex * v1){ V(0) = v0; V(1) = v1; };
	static inline CEdge OrderedEdge(CVertex* v0, CVertex* v1){
		if (v0 < v1) return CEdge(v0, v1);
		else return CEdge(v1, v0);
	}
};

class CFace : public vcg::Face < CUsedTypes,
	vcg::face::InfoOcf,              /* 4b */
	vcg::face::VertexRef,            /*12b */
	vcg::face::BitFlags,             /* 4b */
	vcg::face::Normal3f,             /*12b */
	vcg::face::QualityfOcf,          /* 0b */
	vcg::face::MarkOcf,              /* 0b */
	vcg::face::Color4bOcf,           /* 0b */
	vcg::face::FFAdjOcf,             /* 0b */
	vcg::face::VFAdjOcf,             /* 0b */
	vcg::face::WedgeTexCoordfOcf     /* 0b */
	/*vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef, vcg::face::BitFlags*/ >
{};

class CMesh : public vcg::tri::TriMesh < vcg::vertex::vector_ocf<CVertex>, vcg::face::vector_ocf<CFace>/*, std::vector<CEdge>*/ >
{

};

