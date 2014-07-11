#pragma once
#include "GlobalFunction.h"
#include "PointCloudAlgorithm.h"
#include "normal_extrapolation.h"
#include <iostream>

//#include <CGAL/wlop_simplify_and_regularize_point_set_test_AABB_tree.h>
//#include <CGAL/tags.h>

using namespace std;
using namespace vcg;


// better code is going to be in CGAL 
class WLOP : public PointCloudAlgorithm
{
public:
	WLOP(RichParameterSet* para);
	~WLOP(void);
public:

	void run();
	void setInput(DataMgr* pData);
	RichParameterSet* getParameterSet(){ return para; }
	void setParameterSet(RichParameterSet* _para){ para = _para; }
	void clear();

	void setFirstIterate();
  int getIterateNum(){ return nTimeIterated; }
	double getErrorX(){ return error_x; }

  void computeInitialSampleNeighbor();

  void runSkelWlop();
  void runDragWlop();
  void runRegularizeSamples();

  void computeJointNeighborhood();

protected:
	WLOP(){}

private:
	void input(CMesh* _samples, CMesh* _original);
	void initVertexes();

	double iterate();
	void computeAverageTerm(CMesh* samples, CMesh* original);
  void computeAverageAddSampleTerm(CMesh* samples, CMesh* original);

  void computeSampleAverageTerm(CMesh* samples);
	void computeRepulsionTerm(CMesh* samples);

	void computeDensity(bool isOriginal, double radius);
	void recomputePCA_Normal();

  void stepForward();

  void updateAdaptiveNeighbor();


private:
	RichParameterSet* para;

private:
	CMesh* samples;
	CMesh* original;
  CMesh* dual_samples;

	Box3f box;
	int nTimeIterated;
	double error_x;

	vector<double> samples_density;
	vector<double> original_density;

	vector<Point3f> repulsion;
	vector<double>  repulsion_weight_sum;

  vector<Point3f> samples_average;
  vector<double>  samples_average_weight_sum;

	vector<Point3f> average;
	vector<double>  average_weight_sum;

	vector<CVertex> mesh_temp;

  bool use_adaptive_mu;
  vector<bool> is_sample_close_to_original;
};
