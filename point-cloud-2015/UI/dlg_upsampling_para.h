#pragma once

#include "glarea.h"
#include <QtGui>
#include <QtWidgets/QFrame>
#include <QtWidgets/QWidget>
#include <iostream>
#include <QFileDialog>

#include "..//GeneratedFiles/ui_upsampling_para.h"
#include "ParameterMgr.h"

using namespace std;

class UpsamplingParaDlg : public QFrame
{
	Q_OBJECT
	public:
		UpsamplingParaDlg(QWidget *p, ParameterMgr * _paras, GLArea * _area);
		~UpsamplingParaDlg();
		void initConnects();
		
		void setFrameConent();
	//signals:
	//	void parameterChanged();

	private slots:
		//Common
		bool initWidgets();
		void runAddPts();
		void runProjection();
		void getRadiusValues(double _val);
		void getSigma(double _val);
		void setEdgeParameter(double _val);
		void setNum(int _val);
		
		//Threshold Method
		void setThreshold(double _val);
		void setUsingThresholdProcess(bool _val);
    void setUseConstantThreshold(bool _val);

		void needSnapFiles(bool _val);
    void useAdaptiveUpsampling(bool _val);
    void useUpsampleSkeletalPoints(bool _val);


		void getSnapShotResolution(double _val);
		void getSnapShotIndex(double _val);
		void getBeginIndex(double _val);
		void getEndIndex(double _val);
		void getSpeed(double _val);
		void applyPlayVideo();

    void runPointsExtrapolation();

		void loadVideoFiles();


		void getRotateCenterX(double _val);
		void getRotateCenterY(double _val);
		void getRotateCenterZ(double _val);
		void getRotateNormalX(double _val);
		void getRotateNormalY(double _val);
		void getRotateNormalZ(double _val);
		void getRotateStep(double _val);
		void getRotateAngle(double _val);

    void getDensityThreshold(double _val);



		void rotateStep();
		void rotateAnimation();

	private:
		Ui::Upsampling_para * ui;
		ParameterMgr * m_paras;
		GLArea * area;

		int video_begin_index;
		int video_end_index;
		double video_speed;
};