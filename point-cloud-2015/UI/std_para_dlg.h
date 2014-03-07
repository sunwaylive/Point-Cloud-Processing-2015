#pragma once

#include "glarea.h"
#include <iostream>
#include <QtWidgets/QDialog>
#include <QtGui>
#include <QtWidgets/QFrame>
#include <QtWidgets/QWidget>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QGridLayout>


#include "ParameterMgr.h"
#include "UI/dlg_wlop_para.h"
#include "UI/dlg_normal_para.h"
#include "UI/dlg_skeleton_para.h"
#include "UI/dlg_upsampling_para.h"


using namespace std;

class StdParaDlg : public QDockWidget
{
	Q_OBJECT
	public:
		StdParaDlg(ParameterMgr* _paras, GLArea * _area, QWidget* parent = 0);
		~StdParaDlg();

		bool showWlopParaDialog();
		bool showNormalParaDlg();
		bool showSkeletonParaDlg();
		bool showUpsamplingParaDlg();

	private:
		void init();
		void createFrame();
		void loadWlopFrame();
		void loadNormalFrame();
		void loadSkeletonFrame();
		void loadUpsamplingFrame();

	private slots:
		 void closeClick();
		 
	private:
		WlopParaDlg * para_wlop;
		NormalParaDlg * para_normal;
		SkeletonParaDlg* para_skeleton;
		UpsamplingParaDlg* para_upsampling;
		 
		ParameterMgr * paras;

		QFrame * mainFrame;
		GLArea * gla;
};
