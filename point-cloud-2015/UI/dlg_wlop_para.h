#pragma once


#include "glarea.h"
#include <QtGui>
#include <QtWidgets/QFrame>
#include <QtWidgets/QWidget>
#include <iostream>

#include "ui_para_wlop.h"

#include "ParameterMgr.h"
#include "calculationthread.h"
#include "Algorithm/pointcloud_normal.h"
#include <vcg/complex/trimesh/update/normal.h>

using namespace std;

class WlopParaDlg : public QFrame
{
	Q_OBJECT
	public:
		WlopParaDlg(QWidget *p, ParameterMgr * _paras, GLArea * _area);
		~WlopParaDlg();
		void initConnects();
		void setFrameConent();
	signals:
		void parameterChanged();

	private slots:
		bool initWidgets();
		void getRadiusValues(double _val);
		void getRepPow(double _val);
		void getFitPow(double _val);
		void getIter(int _val);
		void getMu(double _val);
    void getMu3(double _val);
		void getOriginalAverageKNN(double _val);
		void getDualRadius(double _val);

		void get_increasing_step_size(double _val);
		void get_local_neighbor_size(double _val);

		void get_local_neighbor_size_for_surface_points(double _val);
		void get_cooling_parameter(double _val);

		void get_local_angle_threshold(double _val);

		void get_eigen_neighbor_para1(double _val);
		void get_eigen_neighbor_para2(double _val);

		void Average_Closest_Dist(double _val);
		void search_dual_index_para(double _val);
		void Average_Dist_To_Input_Threshold(double _val);
		void Choose_ADT_Threshold_Percentage(double _val);
		void Original_Confidence_KNN(double _val);
		void Similarity_Term_Neighbor_Para(double _val);
		void Similarity_Length_Outlier_Threshold(double _val);
		void Density_Confidence_Segment_Threshold(double _val);
		void Eigen_Directional_Threshold(double _val);
		void Save_Move_Dist_Along_Normal_Para(double _val);
		void Big_Repulsion_Power(double _val);

		void Protect_Small_Tubular_Para(double _val);
		void Protect_High_Confidence_Para(double _val);


		void isDensity(bool _val);
		void isPca(bool _val);
    void needSampleAverage(bool _val);
    void needOriginalCombineSample(bool _val);
    void needAverageMovement(bool _val);

    void useEllipticalOriginalNeighbor(bool _val);
    void useAdaptiveSampleNeighbor(bool _val);
    void useAdaptiveMu(bool _val);
    void useTangentVector(bool _val);
		void useConfidence(bool _val);
		void useOriginalKNN(bool _val);
		void useKitePoints(bool _val);

		void useBackwardFirst(bool _val);
		void useEigenNeighborhood(bool _val);

		void useSeparateNeighborhood(bool _val);
		void useEllipsoidWeight(bool _val);
		void useEllipsoidRepulsion(bool _val);



		void use_Average_Dist_Threshold(bool _val);
		void use_Confidence_To_Combine_Normal(bool _val);
		void use_Only_Do_Repuslion(bool _val);
		void use_Only_Do_Avergage(bool _val);
		void use_Use_Confidence_To_Merge(bool _val);

		//
		void applyWlop();
    void applyDualWlop();
		void applyAnisotropicLop();
    void applyStepForward();
    void applyDualConnection();
		void applyMatLOP();

		void applyRunEstimateAverageDistThreshold();


		void copySamplesToDualSamples();

    void applySkelWlop();
    void applyDragWlop();
    void applyRegularizeSamples();
    void applyRegularizeNormals();
		void applyDetectKitePoints();


    void applyProjection();
    void applyNormalReform();
		void applyComputeConfidence();
		void applyComputeDistribution();
		void applyComputeInnerClustering();
		void applyComputeCorrespondence();
		void applyShowPickDistribution();

		void applyProgressiveNeighborhood();

		void applyInnerPointsClassification();

		void applyEllipsoidFitting();


		void applyInnerPointsRegularization();
		void applySearchNeighborhood();
		void applySmoothNeighborhood();

		void applyMoveBackward();
		void applySelfWLOP();
		void applyNormalSmoothing();

		void applySelfPCA();
		void applySelfPorjection();

		void applyMoveSample();
		void applyMoveSkel();

		void applyComputeEigenDirections();
		void applyComputeEigenNeighbor();

		void applyComputeInitialNeighborhood();
		void oneKEY();

	private:
		bool run_backward_first;

		

		Ui::para_wlop * ui;
		ParameterMgr * m_paras;
		GLArea * area;

		CalculationThread calculation_thread;

};