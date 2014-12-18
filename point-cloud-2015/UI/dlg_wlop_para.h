#pragma once


#include "glarea.h"
#include <QtGui>
#include <QtWidgets/QFrame>
#include <QtWidgets/QWidget>
#include <iostream>

#include "ui_para_wlop.h"

#include "ParameterMgr.h"
#include "calculationthread.h"

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
		void get_local_angle_threshold(double _val);


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

		//
		void applyWlop();
    void applyDualWlop();
		void applyAnisotropicLop();
    void applyStepForward();
    void applyDualConnection();
		void applyMatLOP();

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

	private:
		Ui::para_wlop * ui;
		ParameterMgr * m_paras;
		GLArea * area;

		CalculationThread calculation_thread;

};