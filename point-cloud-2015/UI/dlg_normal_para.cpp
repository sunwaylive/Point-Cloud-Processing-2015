#include "UI/dlg_normal_para.h"

NormalParaDlg::NormalParaDlg(QWidget *p, ParameterMgr * _paras, GLArea * _area) : QFrame(p)
{
	ui = new Ui::normal_paras;
	NormalParaDlg::ui->setupUi(this);
	area = _area;
	m_paras = _paras;

	if(!initWidgets())
	{
		cerr << " NormalParaDlg::initWidgets failed." << endl;
		return;
	}
	initConnects();
}

// 
void NormalParaDlg::initConnects()
{
	connect(ui->radius,SIGNAL(valueChanged(double)),this,SLOT(getRadiusValues(double)));
	connect(ui->sigma,SIGNAL(valueChanged(double)),this,SLOT(getSigma(double)));
	connect(ui->KNN,SIGNAL(valueChanged(int)),this,SLOT(getKNN(int)));
	connect(ui->normal_smoothing,SIGNAL(clicked()),this,SLOT(applyNormalSmoothing()));
	connect(ui->normal_PCA,SIGNAL(clicked()),this,SLOT(applyPCANormal()));
	connect(ui->pushButton_reorientate_normal,SIGNAL(clicked()),this,SLOT(reorientateNormal()));
	connect(ui->use_anistropic_PCA, SIGNAL(clicked(bool)),this,SLOT(isAPCA(bool)));
	connect(ui->smooth_normal_interate_num,SIGNAL(valueChanged(int)),this,SLOT(getIterateNum(int)));
}

// 
bool NormalParaDlg::initWidgets()
{
	ui->radius->setValue(m_paras->norSmooth.getDouble("CGrid Radius"));
	ui->sigma->setValue(m_paras->norSmooth.getDouble("Sharpe Feature Bandwidth Sigma"));
	ui->KNN->setValue(m_paras->norSmooth.getInt("PCA KNN"));
	ui->KNN->setRange(0, 100000);
	ui->smooth_normal_interate_num->setValue(m_paras->norSmooth.getInt("Number Of Iterate"));
	Qt::CheckState state = m_paras->norSmooth.getBool("Run Anistropic PCA") ? (Qt::CheckState::Checked) : (Qt::CheckState::Unchecked);
	ui->use_anistropic_PCA->setCheckState(state);

	return true;
}

void NormalParaDlg::getRadiusValues(double _val)
{
	m_paras->setGlobalParameter("CGrid Radius",DoubleValue(_val));
}

void NormalParaDlg::getPcaThreshold(double _val)
{
	m_paras->norSmooth.setValue("PCA Threshold", DoubleValue(_val));
}

void NormalParaDlg::getIterateNum(int _val)
{
	m_paras->norSmooth.setValue("Number Of Iterate", IntValue(_val));
}

void NormalParaDlg::getSigma(double _val)
{
	m_paras->norSmooth.setValue("Sharpe Feature Bandwidth Sigma",DoubleValue(_val));
}

void NormalParaDlg::getKNN(int _val)
{
	m_paras->norSmooth.setValue("PCA KNN",IntValue(_val));
}



void NormalParaDlg::isAPCA(bool _val)
{
	m_paras->norSmooth.setValue("Run Anistropic PCA", BoolValue(_val));
}


void NormalParaDlg::reorientateNormal()
{
	if (area->dataMgr.isSamplesEmpty())
	{
		return;
	}

	CMesh* samples = area->dataMgr.getCurrentSamples();
	if (global_paraMgr.glarea.getBool("Show Samples"))
	{
		samples = area->dataMgr.getCurrentSamples();
	}
	else if (global_paraMgr.glarea.getBool("Show Dual Samples"))
	{
		samples = area->dataMgr.getCurrentDualSamples();
	}
	else if (global_paraMgr.glarea.getBool("Show Original")
		&& !area->dataMgr.isOriginalEmpty())
	{
		samples = area->dataMgr.getCurrentOriginal();
	}
	else
	{
		samples = area->dataMgr.getCurrentSamples();
	}
	
	for (int i = 0; i < samples->vert.size(); i++)
	{
		samples->vert[i].N() *= -1;
	}
}

void NormalParaDlg::applyNormalSmoothing()
{
	if (global_paraMgr.glarea.getBool("Show Normal"))
	{
		area->runNormalSmoothing();
		area->dataMgr.recomputeQuad();
		area->updateGL();
	}
	else
	{
		int knn = global_paraMgr.norSmooth.getInt("PCA KNN");
		CMesh* samples;

		if (global_paraMgr.glarea.getBool("Show Samples"))
		{
			samples = area->dataMgr.getCurrentSamples();
		}

		vector<Point3f> normals_save(samples->vert.size());
		for (int i = 0; i < samples->vert.size(); i++)
		{
			samples->vert[i].N() = samples->vert[i].N().Normalize();
			normals_save[i] = samples->vert[i].N();
		}

		vcg::tri::PointCloudNormal<CMesh>::Param pca_para;
		pca_para.fittingAdjNum = knn;

		vcg::tri::PointCloudNormal<CMesh>::Compute(*samples, pca_para, NULL);

		for (int i = 0; i < samples->vert.size(); i++)
		{
			if (normals_save[i] * samples->vert[i].N() < 0)
			{
				samples->vert[i].N() *= -1;
			}
		}

	}
}

void NormalParaDlg::applyPCANormal()
{
	bool use_previous_orientation = global_paraMgr.norSmooth.getBool("Run Anistropic PCA");
	if (m_paras->norSmooth.getBool("Run Anistropic PCA"))
	{
		area->runNormalSmoothing();
	}
	else
	{
    int knn = global_paraMgr.norSmooth.getInt("PCA KNN");
    CMesh* samples;

	if (global_paraMgr.glarea.getBool("Show Samples"))
	{
	  samples = area->dataMgr.getCurrentSamples();
	}
    else if (global_paraMgr.glarea.getBool("Show Dual Samples"))
    {
      samples = area->dataMgr.getCurrentDualSamples();
    }
    else if(global_paraMgr.glarea.getBool("Show Original")
            && !area->dataMgr.isOriginalEmpty())
    {
      samples = area->dataMgr.getCurrentOriginal();
    }
    else
    {
      samples = area->dataMgr.getCurrentSamples();
    }
    //vcg::NormalExtrapolation<vector<CVertex> >::ExtrapolateNormals(samples->vert.begin(), samples->vert.end(), knn, -1);
    cout << "new normal estimation: " << endl;

		vector<Point3f> remember_normal(samples->vert.size());
		for (int i = 0; i < samples->vert.size(); i++)
		{
			remember_normal[i] = samples->vert[i].N();
		}
     vcg::tri::PointCloudNormal<CMesh>::Param pca_para;
     pca_para.fittingAdjNum = knn;
 
     vcg::tri::PointCloudNormal<CMesh>::Compute(*samples, pca_para, NULL);

		 if (use_previous_orientation)
		 {
			 for (int i = 0; i < samples->vert.size(); i++)
			 {
				 CVertex& v = samples->vert[i];
				 if (v.N() * remember_normal[i] < 0)
				 {
					 v.N() *= -1;
				 }
			 }
		 }

	}

	cout << "new normal estimation end " << endl;

	area->dataMgr.recomputeQuad();
	area->updateGL();
}

NormalParaDlg::~NormalParaDlg()
{
	delete ui;
	ui = NULL;
	area = NULL;
	m_paras = NULL;
}