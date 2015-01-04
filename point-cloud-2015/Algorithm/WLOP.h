#pragma once
#include "GlobalFunction.h"
#include "PointCloudAlgorithm.h"
#include "normal_extrapolation.h"
#include <iostream>

//#include <CGAL/wlop_simplify_and_regularize_point_set_test_AABB_tree.h>
//#include <CGAL/tags.h>
#include "Algorithm/pointcloud_normal.h"
#include "eigenlib/Eigen/Eigen"

using namespace std;
using namespace vcg;

//typedef Eigen::Matrix<float, 3, 3> Matrix33f;

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

	void computeInitialNeighborSize();
  void computeInitialSampleNeighbor();
	void computeNearestNeighborDist();
	
  void runSkelWlop();
  void runDragWlop();
  void runRegularizeSamples();
  void runRegularizeNormals();
	void runDetectKitePoitns();


	void runComputeDistribution();
	void runEllipsoidFitting();


	void runComputeCorrespondence();
	vector<SphereSlots> computeDistributions(CMesh* samples, CMesh* dual_samples);


	void updateSphereSlots(SphereSlots& sphere_slots, Point3f dir, double dist_weight);
	double getSphereSlotsConfidence(SphereSlots& sphere_slots);
	

  void computeJointNeighborhood();
	void runShowPickDistribution();

	void runProgressiveNeighborhood();

protected:
	WLOP(){}

private:
	void input(CMesh* _samples, CMesh* _original);
	void initVertexes(bool clear_neighbor);

	double iterate();

	vector<Point3f> computeNewSamplePositions(int& error_x);

	void computeAverageTerm(CMesh* samples, CMesh* original);
  void computeAverageAddSampleTerm(CMesh* samples, CMesh* original);

  void computeSampleAverageTerm(CMesh* samples);
	void computeSampleSimilarityTerm(CMesh* samples);
	void computeRepulsionTerm(CMesh* samples);

	void computeDensity(bool isOriginal, double radius);
	void recomputePCA_Normal();

  void stepForward();
	void smoothSkelDistance();

  void updateAdaptiveNeighbor();

  void runProjection();
	void runComputeConfidence();
	void runComputeInnerClusering();

	void runComputeHoleConfidence();

	void runMatLOP();


	void addSamplesToOriginalTemporary();
	void removeSamplesFromOriginal();

	void innerpointsClassification();


	void runInnerPointsRegularization();
	void runSearchNeighborhood();
	void runSmoothNeighborhood();

	void runMoveBackward();
	void runSelfWLOP();
	void runNormalSmoothing();
	void runSelfPCA();
	void runSelfProjection();

	void runMoveSample();
	void runMoeveSkel();

	void runComputeEigenDirections(CMesh* dual_samples, CMesh* samples);
	void runComputeEigenNeighborhood(CMesh* dual_samples, CMesh* samples);



private:
	RichParameterSet* para;

private:
	CMesh* samples;
	CMesh* original;
  CMesh* dual_samples;
	CMesh* target_samples;
	CMesh* target_dual_samples;

	Box3f box;
	int nTimeIterated;
	double error_x;

	double repulsion_radius;

	vector<double> samples_density;
	vector<double> original_density;

	vector<Point3f> repulsion;
	vector<double>  repulsion_weight_sum;

// 	vector<Point3f> repulsion_x;
// 	vector<Point3f> repulsion_y;
// 	vector<Point3f> repulsion_z;

	vector<double> repulsion_x_length;
	vector<double> repulsion_y_length;
	vector<double> repulsion_z_length;


	//Matrix33f matA;
	//vector< Eigen::Matrix3f> repulsion_matA_set;
	vector<Matrix33f> repulsion_matA_set;

  vector<Point3f> samples_average;
  vector<double>  samples_average_weight_sum;

	vector<Point3f> samples_similarity;
	//vector<double>  samples_similarity_weight_sum;

	vector<Point3f> average;
	vector<double>  average_weight_sum;

	vector<Point3f> average_low_confidence;
	vector<double>  average_weight_sum_low_confidence;

	vector<CVertex> mesh_temp;

  bool use_adaptive_mu;
  vector<bool> is_sample_close_to_original;

  vector<int> pioneer_points_id;
  vector<Point3f> pioneer_points_position;
  vector< vector<Point3f>> pioneer_points_origininal;

	SphereSlots default_sphere;
	bool use_closest_dual;
	int added_sample_num;
	bool use_kite_points;
	bool use_eigen_neighborhood;

	bool use_ellipsoid_weight;
	bool use_ellipsoid_repulsion;

public:
	typedef enum { Bilateral, DistanceDiff, NormalDiff, Constant}WeightType;
	vector< vector<double>> neighbor_weights;
	vector< vector<int>> const_surface_neighbors;

	void computeConstNeighborhoodUsingKNN(int knn);
	void computeConstNeighborhoodUsingRadius(double radius);


	void compute_neighbor_weights(vector<CVertex>& samples, 
		                            vector<CVertex>& target,
																vector< vector<double>>& neighbors_weights,
																double radius,
																double sigma,
																WeightType type = Bilateral);

	void compute_neighbor_weights(vector<CVertex>& samples,
		                            vector<CVertex>& target,
		                            vector< vector<int>>& neighbors_indexes,
																vector< vector<double>>& neighbors_weights,
		                            double radius,
		                            double sigma,
																WeightType type = Bilateral);

};
