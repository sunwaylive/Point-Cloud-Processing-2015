#include "UI/dlg_wlop_para.h"

WlopParaDlg::WlopParaDlg(QWidget *p, ParameterMgr * _paras, GLArea * _area) : QFrame(p)
{
	ui = new Ui::para_wlop;
	WlopParaDlg::ui->setupUi(this);
	m_paras = _paras;
	area = _area;

	if(!initWidgets())
	{
		cerr << "Warning:  WlopParaDlg::initWidgets failed!" << endl;
		 return ;
	}

	initConnects();
}


void WlopParaDlg::initConnects()
{
	if (!connect(area,SIGNAL(needUpdateStatus()),this,SLOT(initWidgets())))
	{
		cout << "can not connect signal" << endl;
	}

	if(!connect(ui->radius,SIGNAL(valueChanged(double)),this,SLOT(getRadiusValues(double))))
	{
		cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
	}

	//if(!connect(ui->rep_pow,SIGNAL(valueChanged(double)),this,SLOT(getRepPow(double))))
	//{
	//	cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
	//}
	//if(!connect(ui->fit_pow,SIGNAL(valueChanged(double)),this,SLOT(getFitPow(double))))
	//{
	//	cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
	//}
	if(!connect(ui->mu,SIGNAL(valueChanged(double)),this,SLOT(getMu(double))))
	{
		cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
	}
  if(!connect(ui->sample_average_mu3,SIGNAL(valueChanged(double)),this,SLOT(getMu3(double))))
  {
    cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
  }
	if(!connect(ui->iter,SIGNAL(valueChanged(int)),this,SLOT(getIter(int))))
	{
		cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
	}
	if(!connect(ui->compute_density,SIGNAL(clicked(bool)),this,SLOT(isDensity(bool))))
	{
		cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
	}
	if(!connect(ui->compute_pca,SIGNAL(clicked(bool)),this,SLOT(isPca(bool))))
	{
		cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
	}
  if(!connect(ui->Need_sample_average,SIGNAL(clicked(bool)),this,SLOT(needSampleAverage(bool))))
  {
    cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
  }
  if(!connect(ui->Need_original_combine_sample,SIGNAL(clicked(bool)),this,SLOT(needOriginalCombineSample(bool))))
  {
    cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
  }
  if(!connect(ui->Need_averaging_movement,SIGNAL(clicked(bool)),this,SLOT(needAverageMovement(bool))))
  {
    cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
  }

  if(!connect(ui->Use_Elliptical_Original_Neighbor,SIGNAL(clicked(bool)),this,SLOT(useEllipticalOriginalNeighbor(bool))))
  {
    cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
  }
  if(!connect(ui->Use_KNN_Sample_Neighbor,SIGNAL(clicked(bool)),this,SLOT(useAdaptiveSampleNeighbor(bool))))
  {
    cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
  }
  if(!connect(ui->Use_Adaptive_Mu,SIGNAL(clicked(bool)),this,SLOT(useAdaptiveMu(bool))))
  {
    cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
  }

	//
	if(!connect(ui->wlop_apply,SIGNAL(clicked()),this,SLOT(applyWlop())))
	{
		cerr << "cannot connect WlopParaDlg::applyWlop()." << endl;
	}
  if(!connect(ui->dual_wlop_apply,SIGNAL(clicked()),this,SLOT(applyDualWlop())))
  {
    cerr << "cannot connect WlopParaDlg::applyDualWlop()." << endl;
  }

	connect(ui->anisotropic_lop_apply,SIGNAL(clicked()),this,SLOT(applyAnisotropicLop()));
  connect(ui->step_forward,SIGNAL(clicked()),this,SLOT(applyStepForward()));
  connect(ui->dual_connection,SIGNAL(clicked()),this,SLOT(applyDualConnection()));


}


bool WlopParaDlg::initWidgets()
{
	ui->radius->setValue(m_paras->wLop.getDouble("CGrid Radius"));
	ui->mu->setValue(m_paras->wLop.getDouble("Repulsion Mu"));
	/*ui->rep_pow->setValue(m_paras->wLop.getDouble("Repulsion Power"));
	ui->fit_pow->setValue(m_paras->wLop.getDouble("Average Power"));*/
	ui->iter->setValue(m_paras->wLop.getDouble("Num Of Iterate Time"));
  ui->sample_average_mu3->setValue(m_paras->wLop.getDouble("Sample Average Mu3"));

	
	Qt::CheckState state = m_paras->wLop.getBool("Need Compute Density") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);
	ui->compute_density->setCheckState(state);

	state = m_paras->wLop.getBool("Need Compute PCA") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);
	ui->compute_pca->setCheckState(state);

  state = m_paras->wLop.getBool("Need Sample Average") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);
  ui->Need_sample_average->setCheckState(state);

  state = m_paras->wLop.getBool("Original Combine Sample") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);
  ui->Need_original_combine_sample->setCheckState(state);

  state = m_paras->wLop.getBool("Need Averaging Movement") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);
  ui->Need_averaging_movement->setCheckState(state);

  state = m_paras->wLop.getBool("Use Elliptical Original Neighbor") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);
  ui->Use_Elliptical_Original_Neighbor->setCheckState(state);

  state = m_paras->wLop.getBool("Use Adaptive Sample Neighbor") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);
  ui->Use_KNN_Sample_Neighbor->setCheckState(state);

  state = m_paras->wLop.getBool("Use Adaptive Mu") ? (Qt::CheckState::Checked): (Qt::CheckState::Unchecked);
  ui->Use_Adaptive_Mu->setCheckState(state);

	update();
	repaint();
	return true;
}

void WlopParaDlg::setFrameConent()
{
	if(layout()) delete layout();
	QGridLayout * vLayout = new QGridLayout(this);
	vLayout->setAlignment(Qt::AlignTop);
	setLayout(vLayout);

	showNormal();
	adjustSize();
}

void WlopParaDlg::getRadiusValues(double _val)
{
	m_paras->setGlobalParameter("CGrid Radius",DoubleValue(_val));
	area->updateGL();
}


void WlopParaDlg::getRepPow(double _val)
{
	m_paras->wLop.setValue("Repulsion Power",DoubleValue(_val));
}

void WlopParaDlg::getFitPow(double _val)
{
	m_paras->wLop.setValue("Average Power",DoubleValue(_val));
}

void WlopParaDlg::getIter(int _val)
{
	m_paras->wLop.setValue("Num Of Iterate Time",DoubleValue(_val));
}

void WlopParaDlg::getMu(double _val)
{
	m_paras->wLop.setValue("Repulsion Mu", DoubleValue(_val));
}

void WlopParaDlg::getMu3(double _val)
{
  m_paras->wLop.setValue("Sample Average Mu3", DoubleValue(_val));
}

void WlopParaDlg::isDensity(bool _val)
{
	m_paras->wLop.setValue("Need Compute Density",BoolValue(_val));
	if (_val)
	{
		area->wlop.setFirstIterate();
	}
}

 void WlopParaDlg::needSampleAverage(bool _val)
 {
   m_paras->wLop.setValue("Need Sample Average", BoolValue(_val));

 }

void WlopParaDlg::isPca(bool _val)
{
	m_paras->wLop.setValue("Need Compute PCA",BoolValue(_val));

}

void WlopParaDlg::needOriginalCombineSample(bool _val)
{
  m_paras->wLop.setValue("Original Combine Sample", BoolValue(_val));
}

void WlopParaDlg::needAverageMovement(bool _val)
{
  m_paras->wLop.setValue("Need Averaging Movement", BoolValue(_val));
}

void WlopParaDlg::useEllipticalOriginalNeighbor(bool _val)
{
  m_paras->wLop.setValue("Use Elliptical Original Neighbor", BoolValue(_val));
}

void WlopParaDlg::useAdaptiveSampleNeighbor(bool _val)
{
  m_paras->wLop.setValue("Use Adaptive Sample Neighbor", BoolValue(_val));

  if (_val)
  {
    m_paras->wLop.setValue("Run Compute Initial Sample Neighbor", BoolValue(true));
    area->runWlop();
    m_paras->wLop.setValue("Run Compute Initial Sample Neighbor", BoolValue(false));
  }
}

void WlopParaDlg::useAdaptiveMu(bool _val)
{
  m_paras->wLop.setValue("Use Adaptive Mu", BoolValue(_val));

  if (!_val)
  {
    CMesh* samples = area->dataMgr.getCurrentSamples();

    for (int i = 0; i < samples->vert.size(); i++)
    {
      samples->vert[i].is_fixed_sample = false;
    }
  }
}


// apply
void WlopParaDlg::applyWlop()
{
	//Timer timer;
	//timer.start("WWWWLLLLLOOOOOPPPP Time");
	//area->runWlop();
	//timer.end();

  if (global_paraMgr.glarea.getBool("SnapShot Each Iteration"))
  {
    area->runWlop();
  }
  else
  {
    global_paraMgr.glarea.setValue("Running Algorithm Name", StringValue("WLOP"));
    calculation_thread.setArea(area);
    calculation_thread.start();
  }

}

//void WlopParaDlg::applyDualConnection()
//{
//  CMesh* samples = area->dataMgr.getCurrentSamples();
//  CMesh* dual_samples = area->dataMgr.getCurrentDualSamples();
//
//  GlobalFun::computeAnnNeigbhors(samples->vert, dual_samples->vert, 5, false, "WlopParaDlg::applyDualConnection()");
//  //GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 5, false, "WlopParaDlg::applyDualConnection()");
//
//  cout << "finished ANN " << endl;
//
//  for (int i = 0; i < dual_samples->vert.size(); i++)
//  {
//    int dual_index = dual_samples->vert[i].neighbors[0];
//    if (dual_index < 0 && dual_index >= samples->vert.size())
//    {
//      dual_index = 0;
//    }
//    dual_samples->vert[i].dual_index = dual_index;
//  }
//
//  return;
//}


void WlopParaDlg::applyDualConnection()
{
  CMesh* samples = area->dataMgr.getCurrentSamples();
  CMesh* dual_samples = area->dataMgr.getCurrentDualSamples();

  //GlobalFun::computeAnnNeigbhors(samples->vert, dual_samples->vert, 5, false, "WlopParaDlg::applyDualConnection()");
  GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 5, false, "WlopParaDlg::applyDualConnection()");

  cout << "finished ANN " << endl;

  for (int i = 0; i < samples->vert.size(); i++)
  {
    int dual_index = samples->vert[i].neighbors[0];
    if (dual_index < 0 && dual_index >= dual_samples->vert.size())
    {
      dual_index = 0;
    }
    samples->vert[i].dual_index = dual_index;
  }

  return;
}


 void WlopParaDlg::applyDualWlop()
 {
   double temp_radius = global_paraMgr.wLop.getDouble("CGrid Radius");
   global_paraMgr.setGlobalParameter("CGrid Radius", DoubleValue(temp_radius));

   m_paras->wLop.setValue("Run Dual WLOP", BoolValue(true));
   area->runWlop();
   m_paras->wLop.setValue("Run Dual WLOP", BoolValue(false));

   global_paraMgr.setGlobalParameter("CGrid Radius", DoubleValue(temp_radius));

   //if (global_paraMgr.glarea.getBool("SnapShot Each Iteration"))
   //{
   //  m_paras->wLop.setValue("Run Dual WLOP", BoolValue(true));
   //  area->runWlop();
   //  m_paras->wLop.setValue("Run Dual WLOP", BoolValue(false));
   //}
   //else
   //{ 
   //  m_paras->wLop.setValue("Run Dual WLOP", BoolValue(true));
   //  global_paraMgr.glarea.setValue("Running Algorithm Name", StringValue("WLOP"));
   //  calculation_thread.setArea(area);
   //  calculation_thread.start();
   //}
 }

void WlopParaDlg::applyAnisotropicLop()
{
	if (global_paraMgr.glarea.getBool("SnapShot Each Iteration"))
	{
		m_paras->wLop.setValue("Run Anisotropic LOP", BoolValue(true));
		area->runWlop();
		m_paras->wLop.setValue("Run Anisotropic LOP", BoolValue(false));
	}
	else
	{
		m_paras->wLop.setValue("Run Anisotropic LOP", BoolValue(true));
		global_paraMgr.glarea.setValue("Running Algorithm Name", StringValue("WLOP"));
		calculation_thread.setArea(area);
		calculation_thread.start();
	}

}

void WlopParaDlg::applyStepForward()
{
  m_paras->wLop.setValue("Run Step Forward", BoolValue(true));
  area->runWlop();
  m_paras->wLop.setValue("Run Step Forward", BoolValue(false));

  
}


WlopParaDlg::~WlopParaDlg()
{
	cout << "De-construct WlopParaDlg Frame." << endl;
	delete ui;
	ui = NULL;
}