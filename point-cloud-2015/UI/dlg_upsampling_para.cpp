#include "UI/dlg_upsampling_para.h"

UpsamplingParaDlg::UpsamplingParaDlg(QWidget *p, ParameterMgr * _paras, GLArea * _area) : QFrame(p)
{
	ui = new Ui::Upsampling_para;
	UpsamplingParaDlg::ui->setupUi(this);
	m_paras = _paras;
	area = _area;

	video_begin_index = 0;
	video_end_index = 0;
	video_speed = 1.1;

	if(!initWidgets())
	{
		cout << "Warning:  UpsamplingParaDlg::initWidgets failed!" << endl;
		return ;
	}
	initConnects();
}

void UpsamplingParaDlg::initConnects()
{
	if (!connect(area,SIGNAL(needUpdateStatus()),this,SLOT(initWidgets())))
	{
		cout << "can not connect signal" << endl;
	}

	connect(ui->radius,SIGNAL(valueChanged(double)),this,SLOT(getRadiusValues(double)));
	connect(ui->sigma,SIGNAL(valueChanged(double)),this,SLOT(getSigma(double)));

	connect(ui->add_num,SIGNAL(valueChanged(int)),this,SLOT(setNum(int)));
	connect(ui->using_threshol_process,SIGNAL(clicked(bool)),this,SLOT(setUsingThresholdProcess(bool)));
  connect(ui->use_constant_threshold,SIGNAL(clicked(bool)),this,SLOT(setUseConstantThreshold(bool)));
  connect(ui->use_upsample_skeletal_points, SIGNAL(clicked(bool)), this, SLOT(useUpsampleSkeletalPoints(bool)));

	connect(ui->need_snap_files, SIGNAL(clicked(bool)), this, SLOT(needSnapFiles(bool)));
  connect(ui->use_adaptive_upsampling, SIGNAL(clicked(bool)), this, SLOT(useAdaptiveUpsampling(bool)));


  connect(ui->threshold,SIGNAL(valueChanged(double)),this,SLOT(setThreshold(double)));
	connect(ui->apply_add_point,SIGNAL(clicked()),this,SLOT(runAddPts()));
	connect(ui->pushButton_Projection,SIGNAL(clicked()),this,SLOT(runProjection()));
  connect(ui->Points_extrapolation,SIGNAL(clicked()),this,SLOT(runPointsExtrapolation()));

	connect(ui->edge_paramete,SIGNAL(valueChanged(double)),this,SLOT(setEdgeParameter(double)));

	connect(ui->begin_index,SIGNAL(valueChanged(double)),this,SLOT(getBeginIndex(double)));
	connect(ui->end_index,SIGNAL(valueChanged(double)),this,SLOT(getEndIndex(double)));
	connect(ui->video_speed,SIGNAL(valueChanged(double)),this,SLOT(getSpeed(double)));
	connect(ui->pushButton_play_video,SIGNAL(clicked()),this,SLOT(applyPlayVideo()));
	connect(ui->wlop_snapshot_resolution,SIGNAL(valueChanged(double)),this,SLOT(getSnapShotResolution(double)));
	connect(ui->wlop_snapshot_index,SIGNAL(valueChanged(double)),this,SLOT(getSnapShotIndex(double)));
  connect(ui->density_threshold, SIGNAL(valueChanged(double)), this, SLOT(getDensityThreshold(double)));

	connect(ui->pushButton_load_video_files, SIGNAL(clicked()), this, SLOT(loadVideoFiles()));


	connect(ui->rotate_center_X, SIGNAL(valueChanged(double)), this, SLOT(getRotateCenterX(double)));
	connect(ui->rotate_center_Y, SIGNAL(valueChanged(double)), this, SLOT(getRotateCenterY(double)));
	connect(ui->rotate_center_Z, SIGNAL(valueChanged(double)), this, SLOT(getRotateCenterZ(double)));
	connect(ui->rotate_normal_X, SIGNAL(valueChanged(double)), this, SLOT(getRotateNormalX(double)));
	connect(ui->rotate_normal_Y, SIGNAL(valueChanged(double)), this, SLOT(getRotateNormalY(double)));
	connect(ui->rotate_normal_Z, SIGNAL(valueChanged(double)), this, SLOT(getRotateNormalZ(double)));
	connect(ui->rotate_step, SIGNAL(valueChanged(double)), this, SLOT(getRotateStep(double)));
	connect(ui->rotate_angle, SIGNAL(valueChanged(double)), this, SLOT(getRotateAngle(double)));

	connect(ui->pushButton_rotate, SIGNAL(clicked()), this, SLOT(rotateStep()));
	connect(ui->pushButton_rotate_around, SIGNAL(clicked()), this, SLOT(rotateAnimation()));

}

bool UpsamplingParaDlg::initWidgets()
{
	cout << "init" << endl;
	// 
	ui->radius->setValue(m_paras->upsampling.getDouble("CGrid Radius"));
	ui->sigma->setValue(m_paras->upsampling.getDouble("Feature Sigma"));
	ui->add_num->setValue(m_paras->upsampling.getInt("Number of Add Point"));
	
	Qt::CheckState state = m_paras->upsampling.getBool("Using Threshold Process") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);
	ui->using_threshol_process->setCheckState(state);
	
	//state = m_paras->upsampling.getBool("Auto Recompute Radius For Dist") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);

  state = m_paras->upsampling.getBool("Use Constant Threshold") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);
  ui->use_constant_threshold->setCheckState(state);

	state = m_paras->glarea.getBool("Need Snap Files") ? (Qt::CheckState::Checked) : (Qt::CheckState::Unchecked);
	ui->need_snap_files->setCheckState(state);

  state = m_paras->upsampling.getBool("Use Adaptive Upsampling") ? (Qt::CheckState::Checked) : (Qt::CheckState::Unchecked);
  ui->use_adaptive_upsampling->setCheckState(state);

  state = m_paras->upsampling.getBool("Use Upsample On Skeletal Points") ? (Qt::CheckState::Checked) : (Qt::CheckState::Unchecked);
  ui->use_upsample_skeletal_points->setCheckState(state);

	ui->threshold->setValue(m_paras->upsampling.getDouble("Dist Threshold"));
	ui->edge_paramete->setValue(m_paras->upsampling.getDouble("Edge Parameter"));

	ui->wlop_snapshot_resolution->setValue(m_paras->glarea.getDouble("Snapshot Resolution"));
	ui->wlop_snapshot_index->setValue(m_paras->glarea.getDouble("Snapshot Index"));

  ui->density_threshold->setValue(m_paras->upsampling.getDouble("Density Threshold Dist"));



	ui->rotate_center_X->setValue(area->rotate_pos.X());
	ui->rotate_center_Y->setValue(area->rotate_pos.Y());
	ui->rotate_center_Z->setValue(area->rotate_pos.Z());
	ui->rotate_normal_X->setValue(area->rotate_normal.X());
	ui->rotate_normal_Y->setValue(area->rotate_normal.Y());
	ui->rotate_normal_Z->setValue(area->rotate_normal.Z());
	ui->rotate_step->setValue(area->rotate_delta);
	ui->rotate_angle->setValue(area->rotate_angle);
	//update();
	//repaint();
	return true;
}


void UpsamplingParaDlg::getRadiusValues(double _val)
{
	m_paras->upsampling.setValue("CGrid Radius",DoubleValue(_val));
	m_paras->setGlobalParameter("CGrid Radius", DoubleValue(_val));
	// for debug
	cout << m_paras->upsampling.getDouble("CGrid Radius") << endl;

	area->updateGL();
}


void UpsamplingParaDlg::getSigma(double _val)
{
	m_paras->upsampling.setValue("Feature Sigma",DoubleValue(_val));
}


void UpsamplingParaDlg::runAddPts()
{
	area->runUpsampling();
}

void UpsamplingParaDlg::runProjection()
{
	cout << "void UpsamplingParaDlg::runProjection()" << endl;
	m_paras->upsampling.setValue("Run Projection", BoolValue(true));
	area->runUpsampling();
	m_paras->upsampling.setValue("Run Projection", BoolValue(false));
	area->dataMgr.recomputeQuad();
	area->updateGL();
}


void UpsamplingParaDlg::runPointsExtrapolation()
{
  m_paras->upsampling.setValue("Run Points Extrapolation", BoolValue(true));
  area->runUpsampling();
  m_paras->upsampling.setValue("Run Points Extrapolation", BoolValue(false));
  area->dataMgr.recomputeQuad();javascript:;
  area->updateGL();
}


void UpsamplingParaDlg::setNum(int _val)
{
	m_paras->upsampling.setValue("Number of Add Point",IntValue(_val));
}


//
void UpsamplingParaDlg::setThreshold(double _val)
{
	m_paras->upsampling.setValue("Dist Threshold",DoubleValue(_val));
}

void UpsamplingParaDlg::setUsingThresholdProcess(bool _val)
{
	m_paras->upsampling.setValue("Using Threshold Process", BoolValue(_val));
}

void UpsamplingParaDlg::needSnapFiles(bool _val)
{
	m_paras->glarea.setValue("Need Snap Files", BoolValue(_val));
}

void UpsamplingParaDlg::useUpsampleSkeletalPoints(bool _val)
{
  m_paras->upsampling.setValue("Use Upsample On Skeletal Points", BoolValue(_val));
}


void UpsamplingParaDlg::useAdaptiveUpsampling(bool _val)
{
  m_paras->upsampling.setValue("Use Adaptive Upsampling", BoolValue(_val));

  if (_val)
  {
    CMesh* samples = area->dataMgr.getCurrentSamples();
    CMesh* original = area->dataMgr.getCurrentOriginal();

    GlobalFun::computeAverageDistToInput(samples, original, 50);

    double threshold = m_paras->upsampling.getDouble("Density Threshold Dist");

    for (int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];
      if (v.nearest_neighbor_dist < threshold)
      {
        v.is_fixed_sample = true;
      }
      else
      {
        v.is_fixed_sample = false;
      }
    }
  }
}


void UpsamplingParaDlg::setUseConstantThreshold(bool _val)
{
  m_paras->upsampling.setValue("Use Constant Threshold", BoolValue(_val));

  if (_val)
  {
    m_paras->upsampling.setValue("Run Predict Constant Threshold", BoolValue(true));
    area->runUpsampling();
    m_paras->upsampling.setValue("Run Predict Constant Threshold", BoolValue(false));
  }
}



void UpsamplingParaDlg::setEdgeParameter(double _val)
{
	m_paras->upsampling.setValue("Edge Parameter",DoubleValue(_val));
	
}


void UpsamplingParaDlg::setFrameConent()
{
	if(layout()) delete layout();
	QGridLayout * vLayout = new QGridLayout(this);
	vLayout->setAlignment(Qt::AlignTop);
	setLayout(vLayout);

	showNormal();
	adjustSize();
}

UpsamplingParaDlg::~UpsamplingParaDlg()
{
	cout << "De-construct UpsamplingParaDlg Frame." << endl;
	delete ui;
	ui = NULL;
	area = NULL;
	m_paras = NULL;
}


void UpsamplingParaDlg::getSnapShotResolution(double _val)
{
	m_paras->glarea.setValue("Snapshot Resolution", DoubleValue(_val));
	area->updateGL();
}

void UpsamplingParaDlg::getSnapShotIndex(double _val)
{
	m_paras->glarea.setValue("Snapshot Index", DoubleValue(_val));
	area->updateGL();
}

void UpsamplingParaDlg::getBeginIndex(double _val)
{
	video_begin_index = _val;

	if (_val > 5)
	{
		CMesh* samples = area->dataMgr.getCurrentSamples();

		for (int i = 0; i < samples->vert.size(); i++)
		{
			samples->vert[i].is_skel_ignore = false;
		}

		for (int i = video_begin_index+1; i < samples->vert.size(); i++)
		{
			samples->vert[i].is_skel_ignore = true;
		}
	}


	area->updateGL();
}

void UpsamplingParaDlg::getEndIndex(double _val)
{
	video_end_index = _val;
}

void UpsamplingParaDlg::getSpeed(double _val)
{
	video_speed = _val;
}


void UpsamplingParaDlg::applyPlayVideo()
{
	m_paras->glarea.setValue("SnapShot Each Iteration", BoolValue(true));
	m_paras->glarea.setValue("No Snap Radius",BoolValue(true));

	CMesh* samples = area->dataMgr.getCurrentSamples();

	int begin_index = video_begin_index;
	int end_index = video_end_index;

	if (begin_index > end_index)
	{
		end_index = samples->vn - 1;
	}

	if (end_index >= samples->vn)
	{
		end_index = samples->vn - 1;
	}

	double speed = video_speed;
	double step_size = 2;

	int current_index = begin_index;
	int last_index = current_index;

	for (int i = current_index+1; i < samples->vert.size(); i++)
	{
		samples->vert[i].is_skel_ignore = true;
	}



	while(current_index < end_index)
	{
		area->saveSnapshot();
		area->updateGL();

		step_size *= speed;
		step_size += 1;
		current_index += step_size;

		if (current_index >= samples->vert.size())
		{
			break;
		}

		for (int i = last_index; i < current_index; i++)
		{
			samples->vert[i].is_skel_ignore = false;
		}

		last_index = current_index;
	}

	for (int i = 0; i < end_index; i++)
	{
		samples->vert[i].is_skel_ignore = false;
	}
	area->saveSnapshot();
	area->updateGL();

}


struct naturalSortCompare {

	inline bool isNumber(QChar c) {
		return c >= '0' && c <= '9';
	}

	inline bool operator() (const QString& s1, const QString& s2) {
		if (s1 == "" || s2 == "") return s1 < s2;

		// Move to the first difference between the strings
		int startIndex = -1;
		int length = s1.length() > s2.length() ? s2.length() : s1.length();
		for (int i = 0; i < length; i++) {
			QChar c1 = s1[i];
			QChar c2 = s2[i];
			if (c1 != c2) {
				startIndex = i;
				break;
			}
		}

		// If the strings are the same, exit now.
		if (startIndex < 0) return s1 < s2;

		// Now extract the numbers, if any, from the two strings.
		QString sn1;
		QString sn2;
		bool done1 = false;
		bool done2 = false;
		length = s1.length() < s2.length() ? s2.length() : s1.length();

		for (int i = startIndex; i < length; i++) {
			if (!done1 && i < s1.length()) {
				if (isNumber(s1[i])) {
					sn1 += QString(s1[i]);
				}
				else {
					done1 = true;
				}
			}

			if (!done2 && i < s2.length()) {
				if (isNumber(s2[i])) {
					sn2 += QString(s2[i]);
				}
				else {
					done2 = true;
				}
			}

			if (done1 && done2) break;
		}

		// If none of the strings contain a number, use a regular comparison.
		if (sn1 == "" && sn2 == "") return s1 < s2;

		// If one of the strings doesn't contain a number at that position,
		// we put the string without number first so that, for example,
		// "example.bin" is before "example1.bin"
		if (sn1 == "" && sn2 != "") return true;
		if (sn1 != "" && sn2 == "") return false;

		return sn1.toInt() < sn2.toInt();
	}

};

void UpsamplingParaDlg::loadVideoFiles()
{
	QString file_location = QFileDialog::getExistingDirectory(this, "choose a directory...", "", QFileDialog::ShowDirsOnly);
	if (!file_location.size())
		return;

	QDir dir(file_location);
	if (!dir.exists())
		return;

	dir.setFilter(QDir::Files);
	//dir.setSorting(QDir::Name);
	dir.setSorting(QDir::LocaleAware);

	global_paraMgr.glarea.setValue("Snapshot Index", DoubleValue(1.0));
	global_paraMgr.glarea.setValue("SnapShot Each Iteration", BoolValue(true));


	QFileInfoList list = dir.entryInfoList();
	QStringList str_list;

	for (int i = 0; i < list.size(); ++i)
	{
		QFileInfo fileInfo = list.at(i);
		QString f_name = fileInfo.fileName();
		str_list.push_back(f_name);

	}


	std::sort(str_list.begin(), str_list.end(), naturalSortCompare());

	for (int i = 0; i < str_list.size(); ++i)
	{
		//QFileInfo fileInfo = list.at(i);
		QString f_name = str_list.at(i);

		if (!f_name.endsWith(".skel"))
			continue;

		f_name = file_location + "\\" + f_name;

		cout << f_name.toStdString().c_str() << endl;

		area->dataMgr.loadSkeletonFromSkel(f_name);

		global_paraMgr.glarea.setValue("GLarea Busying", BoolValue(true));

		area->dataMgr.recomputeQuad();

		area->updateUI();
		update();
		area->updateGL();
		emit area->needUpdateStatus();

		area->saveSnapshot();

		//Sleep(1500);
	}
}



void UpsamplingParaDlg::getRotateCenterX(double _val)
{
	area->rotate_pos.X() = _val;
	update(); area->update(); //area->updateGL();
}

void UpsamplingParaDlg::getRotateCenterY(double _val)
{
	area->rotate_pos.Y() = _val;
	update(); area->update(); //area->updateGL();
}

void UpsamplingParaDlg::getRotateCenterZ(double _val)
{
	area->rotate_pos.Z() = _val;
	update(); area->update(); //area->updateGL();
}

void UpsamplingParaDlg::getRotateNormalX(double _val)
{
	area->rotate_normal.X() = _val;
	update(); area->update(); //area->updateGL();
}

void UpsamplingParaDlg::getRotateNormalY(double _val)
{
	area->rotate_normal.Y() = _val;
	update(); area->update(); //area->updateGL();
}

void UpsamplingParaDlg::getRotateNormalZ(double _val)
{
	area->rotate_normal.Z() = _val;
	update(); area->update(); //area->updateGL();
}

void UpsamplingParaDlg::getRotateStep(double _val)
{
	area->rotate_delta = _val;
	update(); area->update(); //area->updateGL();
}

void UpsamplingParaDlg::getRotateAngle(double _val)
{
	area->rotate_angle = _val;
	update(); area->update(); //area->updateGL();
}

void UpsamplingParaDlg::getDensityThreshold(double _val)
{
  m_paras->upsampling.setValue("Density Threshold Dist", DoubleValue(_val));
  update(); area->update(); //area->updateGL();
  area->updateUI();
  area->needUpdateStatus();
}


void UpsamplingParaDlg::rotateStep()
{
	area->rotate_angle += area->rotate_delta;
	if (area->rotate_angle >= 360)
	{
		area->rotate_angle -= 360;
	}
	initWidgets();
	if (m_paras->glarea.getBool("SnapShot Each Iteration"))
	{
		area->saveSnapshot();
	}
	update(); area->update(); area->updateGL();
}

void UpsamplingParaDlg::rotateAnimation()
{
	if (m_paras->glarea.getBool("SnapShot Each Iteration"))
	{
		area->saveSnapshot();
	}

	int rotate_time = 360 / area->rotate_delta;
	for (int i = 0; i < rotate_time; i++)
	{
		rotateStep();
	}
	update(); area->update(); area->updateGL();


}