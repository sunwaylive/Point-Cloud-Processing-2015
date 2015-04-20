#include "mainwindow.h"

MainWindow::MainWindow(QWidget *parent, Qt::WindowFlags flags)
: QMainWindow(parent, flags)
{
	cout << "MainWindow constructed" << endl;
	ui.setupUi(this);
	area = new GLArea(this);
	setCentralWidget(area);

	init();
	initWidgets();
	createActionGroups();
	iniStatusBar();
	initConnect();

	area->loadDefaultModel();
}

MainWindow::~MainWindow()
{
	if(area) delete area;
	area = NULL;

}

void MainWindow::initWidgets()
{

	ui.actionShow_Samples->setChecked(global_paraMgr.glarea.getBool("Show Samples"));
	ui.actionShow_Original->setChecked(global_paraMgr.glarea.getBool("Show Original"));
	ui.actionShow_Normals->setChecked(global_paraMgr.glarea.getBool("Show Normal"));
	ui.actionShow_Neighborhood_Ball->setChecked(global_paraMgr.glarea.getBool("Show Radius"));
	ui.actionShow_All_Raidus->setChecked(global_paraMgr.glarea.getBool("Show All Radius"));	
	ui.actionCull_Points->setChecked(global_paraMgr.drawer.getBool("Need Cull Points"));
	ui.actionShow_Individual_Color->setChecked(global_paraMgr.drawer.getBool("Show Individual Color"));
	ui.actionShow_Confidence_Color->setChecked(global_paraMgr.drawer.getBool("Show Confidence Color"));

	ui.actionShow_Sample_Quads->setChecked(paras->glarea.getBool("Show Samples Quad"));
	ui.actionShow_Sample_Dot->setChecked(paras->glarea.getBool("Show Samples Dot"));
	ui.actionShow_Sample_Circle->setChecked(paras->glarea.getBool("Show Samples Circle"));
	ui.actionShow_Sample_Sphere->setChecked(paras->glarea.getBool("Show Samples Sphere"));

	ui.action_Dual_Quad->setChecked(paras->glarea.getBool("Show Dual Samples Quad"));
	ui.action_Dual_Dot->setChecked(paras->glarea.getBool("Show Dual Samples Dot"));
	ui.action_Dual_Circle->setChecked(paras->glarea.getBool("Show Dual Samples Circle"));
	ui.action_Dual_Sphere->setChecked(paras->glarea.getBool("Show Dual Samples Sphere"));
	
  ui.actionShow_Dual->setChecked(global_paraMgr.glarea.getBool("Show Dual Samples"));
  ui.actionShow_Connection->setChecked(global_paraMgr.glarea.getBool("Show Dual Connection"));
  ui.actionPick_Dual->setChecked(global_paraMgr.glarea.getBool("Pick Dual Point"));

	ui.actionShow_Original_Quad->setChecked(paras->glarea.getBool("Show Original Quad"));
	ui.actionShow_Original_Dot->setChecked(paras->glarea.getBool("Show Original Dot"));
	ui.actionShow_Original_Circle->setChecked(paras->glarea.getBool("Show Original Circle"));
	ui.actionShow_Original_Sphere->setChecked(paras->glarea.getBool("Show Original Sphere"));
	
	ui.actionShow_Normal_Color->setChecked(paras->drawer.getBool("Use Color From Normal"));
	ui.actionSnap_Each_Iteration->setChecked(paras->glarea.getBool("SnapShot Each Iteration"));
	ui.actionNo_Snap_Radius->setChecked(paras->glarea.getBool("No Snap Radius"));
  ui.actionShow_Skeleton->setChecked(paras->glarea.getBool("Show Skeleton"));
  ui.actionShow_colorful_branches->setChecked(paras->drawer.getBool("Use Differ Branch Color"));
  ui.actionShow_Picked_Neighbor->setChecked(paras->drawer.getBool("Draw Picked Point Neighbor"));
	ui.actionShow_Cloest_Dual->setChecked(paras->glarea.getBool("Show Cloest Dual Connection"));
	ui.actionShow_Target->setChecked(paras->glarea.getBool("Show Target Samples"));
	ui.actionShow_Eigens->setChecked(paras->glarea.getBool("Show Eigen Directions"));
	ui.actionMultiple_Pick->setChecked(paras->glarea.getBool("Multiply Pick Point"));
	ui.actionShow_Skeletal_Points->setChecked(paras->glarea.getBool("Show Skeltal Points"));
  ui.actionShow_Neighbor->setChecked(paras->drawer.getBool("Draw Picked Point Neighbor"));

	ui.actionRandom_Erase->setChecked(paras->glarea.getBool("Random Erase"));
	ui.actionShow_Segment_Color->setChecked(paras->drawer.getBool("Show Segmentation Color"));

  update();
  repaint();

}

void MainWindow::initConnect()
{
	if (!connect(area,SIGNAL(needUpdateStatus()),this,SLOT(updateStatusBar())))
	{
		cout << "can not connect signal" << endl;
	}

	if (!connect(area, SIGNAL(needUpdateStatus()), this, SLOT(initWidgets())))
	{
		cout << "can not connect signal" << endl;
	}

	connect(ui.actionImport_Ply, SIGNAL(triggered()), this, SLOT(openFile()));
	connect(ui.actionSave_Ply, SIGNAL(triggered()), this, SLOT(saveFile()));
	connect(ui.actionDownSample, SIGNAL(triggered()), this, SLOT(downSample()));
	connect(ui.actionSubSample, SIGNAL(triggered()), this, SLOT(subSample()));
	connect(ui.actionNormalize, SIGNAL(triggered()), this, SLOT(normalizeData()));
	connect(ui.actionClear_Data, SIGNAL(triggered()), this, SLOT(clearData()));
	connect(ui.actionImport_Image, SIGNAL(triggered()), this, SLOT(openImage()));
	connect(ui.actionSave_View, SIGNAL(triggered()), this, SLOT(saveView()));
	connect(ui.actionSave_Skel, SIGNAL(triggered()), this, SLOT(saveSkel()));
  connect(ui.actionQianSample, SIGNAL(triggered()), this, SLOT(getQianSample()));
	
	connect(ui.actionSnapShot, SIGNAL(triggered()), this, SLOT(saveSnapshot()));
	connect(ui.actionRun_Wlop, SIGNAL(triggered()), this, SLOT(runWLop()));
	connect(ui.actionRun_PCA, SIGNAL(triggered()), this, SLOT(runPCA_Normal()));
	connect(ui.actionReorientate, SIGNAL(triggered()), this, SLOT(reorientateNormal()));
	connect(ui.actionWLOP_Setting, SIGNAL(triggered()), this, SLOT(showWLopDlg()));
	connect(ui.actionNormal_Setting, SIGNAL(triggered()), this, SLOT(showNormalDlg()));
	connect(ui.actionUpsample_Setting, SIGNAL(triggered()), this, SLOT(showUpsampleDlg()));

	connect(ui.actionInitial_Sampling_2, SIGNAL(triggered()), this, SLOT(initialSampling()));
	connect(ui.actionAuto_Play, SIGNAL(triggered()), this, SLOT(autoPlaySkeleton()));
	connect(ui.actionStop, SIGNAL(triggered()), this, SLOT(setStop()));
	connect(ui.actionStep, SIGNAL(triggered()), this, SLOT(stepPlaySkeleton()));
	connect(ui.actionJump, SIGNAL(triggered()), this, SLOT(jumpPlaySkeleton()));
	connect(ui.actionSkeleton_Setting, SIGNAL(triggered()), this, SLOT(showSkeletonDlg()));	
	connect(ui.actionRecompute_Quad,SIGNAL(triggered()),this,SLOT(recomputeQuad()));
	
	connect(ui.actionLight_On_Off, SIGNAL(toggled(bool)), this, SLOT(lightOnOff(bool)));
	connect(ui.actionShow_Individual_Color, SIGNAL(toggled(bool)), this, SLOT(showIndividualColor(bool)));
	connect(ui.actionShow_Samples,SIGNAL(toggled(bool)),this,SLOT(showSamples(bool)));
	connect(ui.actionShow_Original,SIGNAL(toggled(bool)),this,SLOT(showOriginal(bool)));
	connect(ui.actionShow_Normals,SIGNAL(toggled(bool)),this,SLOT(showNormals(bool)));
	connect(ui.actionShow_Neighborhood_Ball,SIGNAL(toggled(bool)),this,SLOT(showNeighborhoodBall(bool)));
	connect(ui.actionShow_All_Raidus,SIGNAL(toggled(bool)),this,SLOT(showAllNeighborhoodBall(bool)));
	connect(ui.actionCull_Points,SIGNAL(toggled(bool)),this,SLOT(cullPoints(bool)));
	connect(ui.actionShow_Normal_Color,SIGNAL(toggled(bool)),this,SLOT(showNormalColor(bool)));
	connect(ui.actionSnap_Each_Iteration,SIGNAL(toggled(bool)),this,SLOT(setSnapshotEachIteration(bool)));
	connect(ui.actionNo_Snap_Radius,SIGNAL(toggled(bool)),this,SLOT(setNoSnapshotWithRadius(bool)));
  connect(ui.actionShow_Skeleton,SIGNAL(toggled(bool)),this,SLOT(showSkeleton(bool)));
  connect(ui.actionShow_colorful_branches,SIGNAL(toggled(bool)),this,SLOT(showColorfulBranches(bool)));
  connect(ui.actionShow_Picked_Neighbor,SIGNAL(toggled(bool)),this,SLOT(showPickPointNeighbor(bool)));
	connect(ui.actionShow_Confidence_Color, SIGNAL(toggled(bool)), this, SLOT(showConfidenceColor(bool)));

  connect(ui.actionShow_Dual,SIGNAL(toggled(bool)),this,SLOT(showDualPoints(bool)));
  connect(ui.actionShow_Connection,SIGNAL(toggled(bool)),this,SLOT(showConnection(bool)));
  connect(ui.actionPick_Dual,SIGNAL(toggled(bool)),this,SLOT(pickDualPoints(bool)));
	connect(ui.actionShow_Eigens, SIGNAL(toggled(bool)), this, SLOT(showEigenDirections(bool)));
	connect(ui.actionMultiple_Pick, SIGNAL(toggled(bool)), this, SLOT(multiplPick(bool)));

	connect(ui.actionRandom_Erase, SIGNAL(toggled(bool)), this, SLOT(randomErasePick(bool)));
	connect(ui.actionShow_Segment_Color, SIGNAL(toggled(bool)), this, SLOT(ShowSegmentColor(bool)));
	connect(ui.actionShow_Skeletal_Points, SIGNAL(toggled(bool)), this, SLOT(showSkeletalPoints(bool)));
  connect(ui.actionShow_Neighbor, SIGNAL(toggled(bool)), this, SLOT(showPickedNeighbors(bool)));

	connect(sample_draw_type,SIGNAL(triggered(QAction *)),this,SLOT(setSmapleType(QAction *)));
	connect(original_draw_type,SIGNAL(triggered(QAction *)),this,SLOT(setOriginalType(QAction *)));
	connect(dual_sample_draw_type, SIGNAL(triggered(QAction *)), this, SLOT(setDualSmapleType(QAction *)));

	connect(ui.actionSample_Color,SIGNAL(triggered()),this,SLOT(sampleColor()));
	connect(ui.actionOriginal_Color,SIGNAL(triggered()),this,SLOT(originalColor()));
	connect(ui.actionBackground_Color,SIGNAL(triggered()),this,SLOT(backGroundColor()));
	connect(ui.actionAmbient_Color,SIGNAL(triggered()),this,SLOT(ambientColor()));
	connect(ui.actionDiffuse_Color,SIGNAL(triggered()),this,SLOT(diffuseColor()));
	connect(ui.actionSpecular_Color,SIGNAL(triggered()),this,SLOT(specularColor()));
	connect(ui.actionNormal_Color,SIGNAL(triggered()),this,SLOT(normalColor()));
	connect(ui.actionFeature_Color,SIGNAL(triggered()),this,SLOT(featureColor()));

	connect(ui.actionErase_Pick,SIGNAL(triggered()),this,SLOT(removePickPoints()));
  connect(ui.actionSwitch_Sample_DualSample,SIGNAL(triggered()),this,SLOT(switchSampleDualSample()));
	connect(ui.actionSwitch_Sample_Original, SIGNAL(triggered()), this, SLOT(switchSampleOriginal()));


  connect(ui.actionReorient_Pick, SIGNAL(triggered()), this, SLOT(reorientPick()));
  connect(ui.actionClean_Pick, SIGNAL(triggered()), this, SLOT(cleanPick()));

	connect(ui.actionShow_Confidence_Color, SIGNAL(toggled(bool)), this, SLOT(pickDualPoints(bool)));
	connect(ui.actionShow_Cloest_Dual, SIGNAL(toggled(bool)), this, SLOT(showClosestDualConnection(bool)));
	connect(ui.actionShow_Target, SIGNAL(toggled(bool)), this, SLOT(showTargets(bool)));

	connect(ui.actionSpray_Erase, SIGNAL(triggered()), this, SLOT(sprayErasePick()));
	connect(ui.actionAdd_Samples_To_Original, SIGNAL(triggered()), this, SLOT(addSamplesToOriginal()));

	QTimer *timer = new QTimer(this);
	connect(timer, SIGNAL(timeout()), this->area, SLOT(update()));
	timer->start(100);
}

void MainWindow::iniStatusBar()
{
	QStatusBar *status_bar = statusBar();

	original_size_label = new QLabel;
	original_size_label->setMinimumSize(200,30);
	original_size_label->setIndent(2);
	original_size_label->setFrameShape(QFrame::NoFrame);
	original_size_label->setFrameShadow(QFrame::Plain);

	downSample_num_label = new QLabel;
	downSample_num_label->setMinimumSize(200,30);
	downSample_num_label->setFrameShape(QFrame::NoFrame);
	downSample_num_label->setFrameShadow(QFrame::Plain);

	radius_label = new QLabel;
	radius_label->setMinimumSize(200,30);
	radius_label->setFrameShape(QFrame::NoFrame);
	radius_label->setFrameShadow(QFrame::Plain);
	
	sample_size_lable = new QLabel;
	sample_size_lable->setMinimumSize(200,30);
	sample_size_lable->setFrameShape(QFrame::NoFrame);
	sample_size_lable->setFrameShadow(QFrame::Plain);

	error_label = new QLabel;
	error_label->setMinimumSize(200,30);
	error_label->setFrameShape(QFrame::NoFrame);
	error_label->setFrameShadow(QFrame::Plain);

  iteration_label = new QLabel;
  iteration_label->setMinimumSize(200,30);
  iteration_label->setFrameShape(QFrame::NoFrame);
  iteration_label->setFrameShadow(QFrame::Plain);

	updateStatusBar();

	status_bar->addWidget(downSample_num_label);
  status_bar->addWidget(iteration_label);
	status_bar->addWidget(error_label);
	status_bar->addWidget(radius_label);
	status_bar->addWidget(original_size_label);
	status_bar->addWidget(sample_size_lable);
}

void MainWindow::updateStatusBar()
{
	QString title = strTitle +  " -- " + area->dataMgr.curr_file_name;
	setWindowTitle(title);

	int o_size = 0;
	if(!area->dataMgr.isOriginalEmpty())
		o_size = area->dataMgr.getCurrentOriginal()->vert.size();
	QString strOriginal = "Original points: " + QString::number(o_size); 
	original_size_label->setText(strOriginal);

	int s_size = 0;
	if(!area->dataMgr.isSamplesEmpty())
		s_size = area->dataMgr.getCurrentSamples()->vert.size();

	QString str = "Sample points: " + QString::number(s_size);
	sample_size_lable->setText(str);

	double nDown = global_paraMgr.data.getDouble("Down Sample Num");
	QString strSub = "Down: " + QString::number(nDown);
	downSample_num_label->setText(strSub);

	double radius  = global_paraMgr.data.getDouble("CGrid Radius");
	QString strRadius = "Radius: " + QString::number(radius);
	radius_label->setText(strRadius);

  int iteration_unm = 0;
  QString running_name = global_paraMgr.glarea.getString("Running Algorithm Name");
  if (running_name == QString("WLOP"))
  {
    iteration_unm = area->wlop.getIterateNum();
  }
  else if (running_name == QString("Skeletonization"))
  {
    iteration_unm = area->skeletonization.getIterateNum();
  }
  QString iterateNum= "Iterate: " + QString::number(iteration_unm);
  iteration_label->setText(iterateNum);


	double error = 0;
	//QString running_name = global_paraMgr.glarea.getString("Running Algorithm Name");
	if (running_name == QString("WLOP"))
	{
		error = global_paraMgr.wLop.getDouble("Current Movement Error");
	}
	else if (running_name == QString("Skeletonization"))
	{
		error = global_paraMgr.skeleton.getDouble("Current Movement Error");
	}
	
	QString strError = "Movement: " + QString::number(error);
	error_label->setText(strError);

  initWidgets();

	update();
	repaint();
}

void MainWindow::createActionGroups()
{
	sample_draw_type = new QActionGroup(this);
	ui.actionShow_Sample_Quads->setActionGroup(sample_draw_type);
	ui.actionShow_Sample_Dot->setActionGroup(sample_draw_type);
	ui.actionShow_Sample_Circle->setActionGroup(sample_draw_type);
	ui.actionShow_Sample_Sphere->setActionGroup(sample_draw_type);

	original_draw_type = new QActionGroup(this);
	ui.actionShow_Original_Quad->setActionGroup(original_draw_type);
	ui.actionShow_Original_Dot->setActionGroup(original_draw_type);
	ui.actionShow_Original_Circle->setActionGroup(original_draw_type);
	ui.actionShow_Original_Sphere->setActionGroup(original_draw_type);

	dual_sample_draw_type = new QActionGroup(this);
	ui.action_Dual_Quad->setActionGroup(dual_sample_draw_type);
	ui.action_Dual_Dot->setActionGroup(dual_sample_draw_type);
	ui.action_Dual_Circle->setActionGroup(dual_sample_draw_type);
	ui.action_Dual_Sphere->setActionGroup(dual_sample_draw_type);

	QString str = strTitle + " -- Welcome!";
	setWindowTitle(str);
}



void MainWindow::init()
{
	strTitle = "Point Cloud";
	paraDlg_Skeleton = NULL;
	paraDlg_Upsample = NULL;
	paraDlg_WLOP = NULL;
	paraDlg_Normal = NULL;
	paras = &global_paraMgr;
}

void MainWindow::showWLopDlg()
{

	if(paraDlg_WLOP != 0)
	{
		paraDlg_WLOP->close();
		delete paraDlg_WLOP;
	}

	paraDlg_WLOP = new StdParaDlg(paras, area, this);
	paraDlg_WLOP->setAllowedAreas(Qt::LeftDockWidgetArea
		| Qt::RightDockWidgetArea);
	addDockWidget(Qt::RightDockWidgetArea,paraDlg_WLOP);
	paraDlg_WLOP->setFloating(false);
	paraDlg_WLOP->hide();
	paraDlg_WLOP->showWlopParaDialog();
}

void MainWindow::showNormalDlg()
{
	if(paraDlg_Normal != 0)
	{
		paraDlg_Normal->close();
		delete paraDlg_Normal;
	}

	paraDlg_Normal = new StdParaDlg(paras, area, this);
	paraDlg_Normal->setAllowedAreas(Qt::LeftDockWidgetArea
		| Qt::RightDockWidgetArea);
	addDockWidget(Qt::RightDockWidgetArea,paraDlg_Normal);
	paraDlg_Normal->setFloating(false);
	paraDlg_Normal->hide();
	paraDlg_Normal->showNormalParaDlg();
}

void MainWindow::showUpsampleDlg()
{
	if(paraDlg_Upsample != 0)
	{
		paraDlg_Upsample->close();
		delete paraDlg_Upsample;
	}

	paraDlg_Upsample = new StdParaDlg(paras, area, this);
	paraDlg_Upsample->setAllowedAreas(Qt::LeftDockWidgetArea
		| Qt::RightDockWidgetArea);
	addDockWidget(Qt::LeftDockWidgetArea,paraDlg_Upsample);
	paraDlg_Upsample->setFloating(false);
	paraDlg_Upsample->hide();
	paraDlg_Upsample->showUpsamplingParaDlg();

}

void MainWindow::showSkeletonDlg()
{
	if(paraDlg_Skeleton != 0)
	{
		paraDlg_Skeleton->close();
		delete paraDlg_Skeleton;
	}

	paraDlg_Skeleton = new StdParaDlg(paras, area, this);
	paraDlg_Skeleton->setAllowedAreas(Qt::LeftDockWidgetArea
		| Qt::RightDockWidgetArea);
	addDockWidget(Qt::LeftDockWidgetArea,paraDlg_Skeleton);
	paraDlg_Skeleton->setFloating(false);
	paraDlg_Skeleton->hide();
	paraDlg_Skeleton->showSkeletonParaDlg();	
}

void MainWindow::autoPlaySkeleton()
{

}

void MainWindow::stepPlaySkeleton()
{

}

void MainWindow::jumpPlaySkeleton()
{

}

void MainWindow::initialSampling()
{

}

void MainWindow::setStop()
{

}

void MainWindow::openFile()
{
	QString file = QFileDialog::getOpenFileName(this, "Select a ply file", "", "*.ply");
	if(!file.size()) return;

	area->dataMgr.loadPlyToSample(file);
	area->initAfterOpenFile();
	area->updateGL();
}

void MainWindow::openImage()
{
	QString file = QFileDialog::getOpenFileName(this, "Select a ply file", "", "");
	if(!file.size()) return;

	area->dataMgr.loadImage(file);
	area->initAfterOpenFile();
	area->updateGL();
}

void MainWindow::saveFile()
{
	QString file = QFileDialog::getSaveFileName(this, "Save samples as", "", "*.ply");
	if(!file.size()) return;

	if (!file.endsWith(".ply"))
	{
		file += QString(".ply");
	}
	area->cleanPickPoints();

	if (global_paraMgr.glarea.getBool("Show Original"))
	{
		if (global_paraMgr.glarea.getBool("Show Samples"))
		{
			area->dataMgr.eraseRemovedSamples();
			area->dataMgr.savePly(file, *area->dataMgr.getCurrentSamples());
		}

		if (global_paraMgr.glarea.getBool("Show Dual Samples"))
		{
			//area->dataMgr.eraseRemovedSamples();
			file.replace(".ply", "_dual.ply");
			area->dataMgr.savePly(file, *area->dataMgr.getCurrentDualSamples());
			file.replace("_dual.ply", ".ply");
		}

		if (file.endsWith("ply") && !area->dataMgr.isOriginalEmpty())
		{
			file.replace(".ply", "_original.ply");
			area->dataMgr.savePly(file, *area->dataMgr.getCurrentOriginal());
		}

		return;
	}

	if (global_paraMgr.glarea.getBool("Show Dual Samples"))
	{
		//area->dataMgr.eraseRemovedSamples();
		file.replace(".ply", "_dual.ply");
		area->dataMgr.savePly(file, *area->dataMgr.getCurrentDualSamples());
		file.replace("_dual.ply", ".ply");
	}

	if (file.endsWith("ply") && !area->dataMgr.isSamplesEmpty())
	{
		area->dataMgr.savePly(file, *area->dataMgr.getCurrentSamples());
	}
}

void MainWindow::downSample()
{
  if (global_paraMgr.glarea.getBool("GLarea Busying"))
  {
    global_paraMgr.glarea.setValue("Algorithom Stop", BoolValue(true));
    global_paraMgr.glarea.setValue("GLarea Busying", BoolValue(false));
    return;
  }

  
	area->dataMgr.downSamplesByNum();
  area->dataMgr.skeleton.clear();
	area->initSetting();
	area->updateGL();
}

void MainWindow::getQianSample()
{
 /* CMesh* samples = area->dataMgr.getCurrentSamples();
  CMesh* original = area->dataMgr.getCurrentOriginal();

  if (!original->vert.empty())
  {
    original->vert.clear();
  }
  for (int i = 0; i < samples->vert.size(); i++)
  {
    CVertex v = samples->vert[i];
    v.bIsOriginal = true;
    original->vert.push_back(v);
  }
  original->vn = samples->vn;*/

  area->dataMgr.downSamplesByNum(false);
  area->initAfterOpenFile();
  //area->initSetting();
  area->updateGL();

}

void MainWindow::subSample()
{
	area->dataMgr.subSamples();
	area->initSetting();
	area->updateGL();
}

void MainWindow::normalizeData()
{
	area->dataMgr.normalizeAllMesh();


  QString file = QFileDialog::getSaveFileName(this, "Select a ply file", "", "*.normalize");
  if (!file.size()) return;

  if (!file.endsWith(".normalize"))
  {
    file += QString(".normalize");
  }

  area->dataMgr.saveNomalization(file);

	area->initView();
	area->updateGL();
}

void MainWindow::clearData()
{
	area->dataMgr.clearData();
	area->updateUI();
	area->updateGL();
}


void MainWindow::saveSnapshot()
{
//   Timer timer;
//   timer.start("walk");
//   CMesh* dual_samples = area->dataMgr.getCurrentDualSamples();
//   GlobalFun::computeAnnNeigbhors(dual_samples->vert, dual_samples->vert, 350, false, "test random walk");
//   timer.end();
//   return;

	area->saveSnapshot();
}

void MainWindow::saveView()
{
   Timer timer;
   timer.start("walk");
   CMesh* dual_samples = area->dataMgr.getCurrentDualSamples();
   GlobalFun::computeRandomwalkNeighborhood(dual_samples, 6, 250);
   timer.end();

  return;

  CMesh* samples = area->dataMgr.getCurrentSamples();
  double percentage = 0.02;
  GlobalFun::removeOutliersBaseOnDistance(samples, 10, percentage);
  //GlobalFun::removeOutliersBaseOnNormal(samples, 50, percentage);


// 	ofstream outpara("outpara.paras");
// 	global_paraMgr.outputAllParameters(outpara);
// 	outpara.close();
// 
//  	QString file = QFileDialog::getSaveFileName(this, "Select a ply file", "", "*.VPoint");
//  	if(!file.size()) return;
//  	area->saveView(file);
// 
// 	file.replace(".VPoint", ".paras");
// 	ifstream inpara(file.toStdString().c_str());
// 	global_paraMgr.inputAllParameters(inpara);
// 	inpara.close();
}

void MainWindow::saveSkel()
{
	QString file = QFileDialog::getSaveFileName(this, "Save samples as", "", "*.skel");
	if(!file.size()) return;

	if (!file.endsWith(".skel"))
	{
		file += QString(".skel");
	}

	if (global_paraMgr.glarea.getBool("Show Skeleton"))
	{
		//area->dataMgr.saveTargetSkeletonAsSkel(file);
		area->dataMgr.saveSkeletonAsSkel(file);
	}
	else
	{
		area->dataMgr.saveSkeletonAsSkel(file);
	}
// 
// 	file.replace(".skel", ".View");
// 	area->saveView(file);

	file.replace(".skel", ".paras");
	ofstream outpara(file.toStdString().c_str());
	global_paraMgr.outputAllParameters(outpara);
	outpara.close();

	file.replace(".paras", ".VPoint");
	area->saveVPoint(file);

	area->updateGL();
}


void MainWindow::setSnapshotEachIteration(bool _val)
{
	paras->glarea.setValue("SnapShot Each Iteration",BoolValue(_val));
}

void MainWindow::setNoSnapshotWithRadius(bool _val)
{
	paras->glarea.setValue("No Snap Radius", BoolValue(_val));
}

void MainWindow::showPickPointNeighbor(bool _val)
{
  paras->drawer.setValue("Draw Picked Point Neighbor", BoolValue(_val));
}

void MainWindow::showEigenDirections(bool _val)
{
	paras->glarea.setValue("Show Eigen Directions", BoolValue(_val));
}

void MainWindow::showSkeletalPoints(bool _val)
{
	paras->glarea.setValue("Show Skeltal Points", BoolValue(_val));
}

void MainWindow::showPickedNeighbors(bool _val)
{
  paras->drawer.setValue("Draw Picked Point Neighbor", BoolValue(_val));
}


void MainWindow::showColorfulBranches(bool _val)
{
  //if (area->paraMgr.glarea.getBool("GLarea Busying"))
  //{
  //  return;
  //}

  if (_val)
  {
    area->glDrawer.generateRandomColorList();
  }
  paras->drawer.setValue("Use Differ Branch Color", BoolValue(_val));
  area->updateGL();
}


void MainWindow::showTargets(bool _val)
{
	paras->glarea.setValue("Show Target Samples", BoolValue(_val));
	paras->glarea.setValue("Show Target Dual Samples", BoolValue(_val));

	paras->drawer.setValue("Show Feature Color", BoolValue(_val));

	area->updateGL();
}


void MainWindow::showClosestDualConnection(bool _val)
{
	global_paraMgr.wLop.setValue("Dual Samples Represent Skeltal Points", BoolValue(true));
  global_paraMgr.glarea.setValue("Show Cloest Dual Connection", BoolValue(_val));

	if (_val && !area->dataMgr.isSkeletalPointsEmpty())
	{
		global_paraMgr.wLop.setValue("Run Compute Dual Index", BoolValue(true));
		area->runWlop();
		global_paraMgr.wLop.setValue("Run Compute Dual Index", BoolValue(false));
	}
	else
	{
		CMesh* samples = area->dataMgr.getCurrentSamples();

		for (int i = 0; i < samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];
			v.dual_index = i;
		}
	}

	emit area->needUpdateStatus();

// 	CMesh* samples = area->dataMgr.getCurrentSamples();
// 	CMesh* dual_samples = area->dataMgr.getCurrentDualSamples();
// 
// 	paras->glarea.setValue("Show Cloest Dual Connection", BoolValue(_val));
// 
// 	GlobalFun::computeBallNeighbors(dual_samples, NULL, global_paraMgr.wLop.getDouble("CGrid Radius"), samples->bbox);
// 
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		CVertex& v = samples->vert[i];
// 		CVertex& dual_v = dual_samples->vert[i];
// 
// 		double min_dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_v.P());
// 		int dual_idx = i;
// 		for (int j = 0; j < dual_v.neighbors.size(); j++)
// 		{
// 			int index = dual_v.neighbors[j];
// 			CVertex& dual_t = dual_samples->vert[index];
// 
// 			double dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_t.P());
// 			if (dist2 < min_dist2)
// 			{
// 				min_dist2 = dist2;
// 				dual_idx = index;
// 			}
// 		}
// 
// 		v.dual_index = dual_idx;
// 	}
// 
//  	if (_val)
//  	{
//  		//GlobalFun::computeAnnNeigbhors(samples->vert, dual_samples->vert, 1, false, "showClosestDualConnection");
//  		GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 1, false, "showClosestDualConnection");
//  
//  		for (int i = 0; i < samples->vert.size(); i++)
//  		{
//  			CVertex& v = samples->vert[i];
//  			v.dual_index = v.neighbors[0];
//  		}
//  	}

}

void MainWindow::lightOnOff(bool _val)
{
	paras->glarea.setValue("Light On or Off",BoolValue(_val));
	if (_val)
	{
		glEnable(GL_LIGHTING);
	}
	else
	{
		glDisable(GL_LIGHTING);
	}
	area->updateGL();
}

void MainWindow::showOriginal(bool _val)
{
	global_paraMgr.glarea.setValue("Show Original", BoolValue(_val));
	area->updateGL();
}

void MainWindow::showSamples(bool _val)
{
	global_paraMgr.glarea.setValue("Show Samples", BoolValue(_val));
	area->updateGL();
}

void MainWindow::showDualPoints(bool _val)
{
  global_paraMgr.glarea.setValue("Show Dual Samples", BoolValue(_val));
  area->updateGL();
}

void MainWindow::showConnection(bool _val)
{
  global_paraMgr.glarea.setValue("Show Dual Connection", BoolValue(_val));
  area->updateGL();
}

void MainWindow::pickDualPoints(bool _val)
{
  global_paraMgr.glarea.setValue("Pick Dual Point", BoolValue(_val));
  area->updateGL();
}

 void MainWindow::showSkeleton(bool _val)
 {
   global_paraMgr.glarea.setValue("Show Skeleton", BoolValue(_val));
   area->updateGL();
 }

void MainWindow::showNormals(bool _val)
{
	global_paraMgr.glarea.setValue("Show Normal", BoolValue(_val));
	area->updateGL();
}

void MainWindow::cullPoints(bool _val)
{
	global_paraMgr.drawer.setValue("Need Cull Points", BoolValue(_val));
	area->updateGL();
}

void MainWindow::showNormalColor(bool _val)
{
	cout << "show normal" << endl;
	global_paraMgr.drawer.setValue("Use Color From Normal", BoolValue(_val));
	area->updateGL();
}


void MainWindow::showNeighborhoodBall(bool _val)
{
	global_paraMgr.glarea.setValue("Show Radius", BoolValue(_val));
	area->updateGL();
}

void MainWindow::showAllNeighborhoodBall(bool _val)
{
	global_paraMgr.glarea.setValue("Show All Radius", BoolValue(_val));
	area->updateGL();
}

void MainWindow::showIndividualColor(bool _val)
{
	global_paraMgr.drawer.setValue("Show Individual Color", BoolValue(_val));
	area->updateGL();
}

void MainWindow::setSmapleType(QAction * action)
{
	if(action == ui.actionShow_Sample_Quads)
	{
		paras->glarea.setValue("Show Samples Quad",BoolValue(true));
		paras->glarea.setValue("Show Samples Dot",BoolValue(false));
		paras->glarea.setValue("Show Samples Circle",BoolValue(false));
		paras->glarea.setValue("Show Samples Sphere",BoolValue(false));
	}
	else if (action == ui.actionShow_Sample_Dot)
	{
		paras->glarea.setValue("Show Samples Quad",BoolValue(false));
		paras->glarea.setValue("Show Samples Dot",BoolValue(true));
		paras->glarea.setValue("Show Samples Circle",BoolValue(false));
		paras->glarea.setValue("Show Samples Sphere",BoolValue(false));
	}
	else if(action == ui.actionShow_Sample_Circle)
	{
		paras->glarea.setValue("Show Samples Quad",BoolValue(false));
		paras->glarea.setValue("Show Samples Dot",BoolValue(false));
		paras->glarea.setValue("Show Samples Circle",BoolValue(true));
		paras->glarea.setValue("Show Samples Sphere",BoolValue(false));
	}
	else if (action == ui.actionShow_Sample_Sphere)
	{
		paras->glarea.setValue("Show Samples Quad",BoolValue(false));
		paras->glarea.setValue("Show Samples Dot",BoolValue(false));
		paras->glarea.setValue("Show Samples Circle",BoolValue(false));
		paras->glarea.setValue("Show Samples Sphere",BoolValue(true));
	}
	area->updateGL();
}


void MainWindow::setDualSmapleType(QAction * action)
{
	if (action == ui.action_Dual_Quad)
	{
		paras->glarea.setValue("Show Dual Samples Quad", BoolValue(true));
		paras->glarea.setValue("Show Dual Samples Dot", BoolValue(false));
		paras->glarea.setValue("Show Dual Samples Circle", BoolValue(false));
		paras->glarea.setValue("Show Dual Samples Sphere", BoolValue(false));
	}
	else if (action == ui.action_Dual_Dot)
	{
		paras->glarea.setValue("Show Dual Samples Quad", BoolValue(false));
		paras->glarea.setValue("Show Dual Samples Dot", BoolValue(true));
		paras->glarea.setValue("Show Dual Samples Circle", BoolValue(false));
		paras->glarea.setValue("Show Dual Samples Sphere", BoolValue(false));
	}
	else if (action == ui.action_Dual_Circle)
	{
		paras->glarea.setValue("Show Dual Samples Quad", BoolValue(false));
		paras->glarea.setValue("Show Dual Samples Dot", BoolValue(false));
		paras->glarea.setValue("Show Dual Samples Circle", BoolValue(true));
		paras->glarea.setValue("Show Dual Samples Sphere", BoolValue(false));
	}
	else if (action == ui.action_Dual_Sphere)
	{
		paras->glarea.setValue("Show Dual Samples Quad", BoolValue(false));
		paras->glarea.setValue("Show Dual Samples Dot", BoolValue(false));
		paras->glarea.setValue("Show Dual Samples Circle", BoolValue(false));
		paras->glarea.setValue("Show Dual Samples Sphere", BoolValue(true));
	}
	area->updateGL();
}


void MainWindow::setOriginalType(QAction * action)
{
	if(action == ui.actionShow_Original_Quad)
	{
		paras->glarea.setValue("Show Original Quad",BoolValue(true));
		paras->glarea.setValue("Show Original Dot",BoolValue(false));
		paras->glarea.setValue("Show Original Circle",BoolValue(false));
		paras->glarea.setValue("Show Original Sphere",BoolValue(false));
	}
	else if (action == ui.actionShow_Original_Dot)
	{
		paras->glarea.setValue("Show Original Quad",BoolValue(false));
		paras->glarea.setValue("Show Original Dot",BoolValue(true));
		paras->glarea.setValue("Show Original Circle",BoolValue(false));
		paras->glarea.setValue("Show Original Sphere",BoolValue(false));
	}
	else if(action == ui.actionShow_Original_Circle)
	{
		paras->glarea.setValue("Show Original Quad",BoolValue(false));
		paras->glarea.setValue("Show Original Dot",BoolValue(false));
		paras->glarea.setValue("Show Original Circle",BoolValue(true));
		paras->glarea.setValue("Show Original Sphere",BoolValue(false));
	}
	else if (action == ui.actionShow_Original_Sphere)
	{
		paras->glarea.setValue("Show Original Quad",BoolValue(false));
		paras->glarea.setValue("Show Original Dot",BoolValue(false));
		paras->glarea.setValue("Show Original Circle",BoolValue(false));
		paras->glarea.setValue("Show Original Sphere",BoolValue(true));
	}
	area->updateGL();
}

void MainWindow::runWLop()
{
	global_paraMgr.glarea.setValue("Running Algorithm Name", StringValue("WLOP"));
	calculation_thread.setArea(area);
	calculation_thread.start();

}

void MainWindow::runPCA_Normal()
{
	int knn = global_paraMgr.norSmooth.getInt("PCA KNN");
	CMesh* samples = area->dataMgr.getCurrentSamples();
	vcg::NormalExtrapolation<vector<CVertex> >::ExtrapolateNormals(samples->vert.begin(), samples->vert.end(), knn, -1);
}

void MainWindow::reorientateNormal()
{
	if (area->dataMgr.isSamplesEmpty())
	{
		return;
	}

	CMesh* samples = area->dataMgr.getCurrentSamples();
	for (int i = 0; i < samples->vert.size(); i++)
	{
		samples->vert[i].N() *= -1;
	}
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
	event->accept();
}

void MainWindow::dropEvent ( QDropEvent * event )
{
	const QMimeData * data = event->mimeData();

	if (data->hasUrls())
	{
		QList< QUrl > url_list = data->urls();
		for (int i=0, size=url_list.size(); i<size; i++)
		{
			QString path = url_list.at(i).toLocalFile();
			area->openByDrop(path);
			cout << "open file: "<< path.toStdString() << endl;
		}
	}
}

void MainWindow::sampleColor()
{
	area->changeColor("Sample Point Color");
	area->updateGL();
}

void MainWindow::originalColor()
{
	area->changeColor("Original Point Color");
	area->updateGL();
}

void MainWindow::backGroundColor()
{
	//area->changeColor("Background Color");
	area->changeColor("Backface Color");

	area->updateGL();
}

void MainWindow::ambientColor()
{
	area->changeColor("Light Ambient Color");
	area->updateGL();
}

void MainWindow::diffuseColor()
{
	area->changeColor("Light Diffuse Color");
	area->updateGL();
}

void MainWindow::specularColor()
{
	area->changeColor("Light Specular Color");
	area->updateGL();
}

void MainWindow::normalColor()
{
	area->changeColor("Normal Line Color");
	area->updateGL();
}

void MainWindow::featureColor()
{
	//area->changeColor("Feature Color");
	area->changeColor("DLink Color");

	area->updateGL();
}

void MainWindow::recomputeQuad()
{
	//cout << "recompute quad" << endl;
// 	if (area->dataMgr.isSamplesEmpty())
// 	{
// 		return;
// 	}

	area->dataMgr.recomputeQuad();
	area->updateGL();
}


void MainWindow::removePickPoints()
{
	area->removePickPoint();
	area->updateUI();
	area->updateGL();
}

void MainWindow::switchSampleDualSample()
{
  area->cleanPickPoints();
  area->dataMgr.switchSampleDualSample();
  area->updateUI();
}

void MainWindow::switchSampleOriginal()
{
	area->cleanPickPoints();
	area->dataMgr.switchSampleOriginal();
	area->updateUI();
}


void MainWindow::reorientPick()
{
  area->reorientPick();
  area->updateUI();
  area->updateGL();
}

void MainWindow::cleanPick()
{
  area->cleanPick();
  area->updateUI();
  area->updateGL();
}

void MainWindow::showConfidenceColor(bool _val)
{
	global_paraMgr.drawer.setValue("Show Confidence Color", BoolValue(_val));
	area->updateGL();
}

void MainWindow::multiplPick(bool _val)
{
	global_paraMgr.glarea.setValue("Multiply Pick Point", BoolValue(_val));
	area->updateGL();
}

void MainWindow::randomErasePick(bool _val)
{
	global_paraMgr.glarea.setValue("Random Erase", BoolValue(_val));
	area->updateGL();
}

void MainWindow::ShowSegmentColor(bool _val)
{
	global_paraMgr.drawer.setValue("Show Segmentation Color", BoolValue(_val));
	area->updateGL();
}

void MainWindow::sprayErasePick()
{
	area->sprayErasePick();
	area->updateGL();
}

void MainWindow::addSamplesToOriginal()
{
	CMesh* samples = area->dataMgr.getCurrentSamples();
	CMesh* original = area->dataMgr.getCurrentOriginal();

	samples->face.clear();
	original->face.clear();

	int idx = original->vert.back().m_index + 1;

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex t = samples->vert[i];
		t.bIsOriginal = true;
		t.m_index = idx++;

		original->vert.push_back(t);
		original->bbox.Add(t.P());
	}
	original->vn = original->vert.size();
}


