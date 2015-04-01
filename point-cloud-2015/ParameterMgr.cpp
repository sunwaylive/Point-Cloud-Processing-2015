#include "ParameterMgr.h"
#include <iostream>

int ParameterMgr::init_time = 0;
ParameterMgr global_paraMgr;

ParameterMgr::ParameterMgr(void)
{
	init_time++;
	if(init_time > 1)
	{
		std::cout << "can not init ParameterMgr twice!" << endl;
		return;
	}

	grid_r = -0.2;


	initDataMgrParameter();
	initDrawerParameter();		
	initGlareaParameter();
	initWLopParameter();
	initNormalSmootherParameter();
	initSkeletonParameter();
	initUpsamplingParameter();
}

ParameterMgr::~ParameterMgr(void)
{
}





void ParameterMgr::setGlobalParameter(QString paraName,Value& val)
{
	if(glarea.hasParameter(paraName))
		glarea.setValue(paraName, val);
	if(data.hasParameter(paraName))
		data.setValue(paraName, val);
	if(drawer.hasParameter(paraName))
		drawer.setValue(paraName, val);
	if(wLop.hasParameter(paraName))
		wLop.setValue(paraName, val);
	if(norSmooth.hasParameter(paraName))
		norSmooth.setValue(paraName, val);
	if (skeleton.hasParameter(paraName))
		skeleton.setValue(paraName, val);
	if (upsampling.hasParameter(paraName))
		upsampling.setValue(paraName, val);

	if (!glarea.hasParameter(paraName) && 
		  !data.hasParameter(paraName) &&
			!wLop.hasParameter(paraName) &&
			!norSmooth.hasParameter(paraName) &&
			!skeleton.hasParameter(paraName) &&
			!upsampling.hasParameter(paraName) &&
			!drawer.hasParameter(paraName))
	{
		std::cout << "out-of-dated para name: " << paraName.toStdString().c_str() << std::endl;
	}
}

void ParameterMgr::initDataMgrParameter()
{
	data.addParam(new RichDouble("Init Radius Para", 1.0));
	data.addParam(new RichDouble("Down Sample Num", 1000));
	data.addParam(new RichDouble("CGrid Radius", grid_r));
}


void ParameterMgr::initGlareaParameter()
{
	glarea.addParam(new RichString("Running Algorithm Name", "") );

	glarea.addParam(new RichBool("Light On or Off", false) );

	glarea.addParam(new RichBool("Show Normal", false) );
	
	glarea.addParam(new RichBool("Show Samples", true) );
	glarea.addParam(new RichBool("Show Samples Quad", false) );
	glarea.addParam(new RichBool("Show Samples Dot", true) );
	glarea.addParam(new RichBool("Show Samples Circle", false) );
	glarea.addParam(new RichBool("Show Samples Sphere", false) );

  glarea.addParam(new RichBool("Show Dual Samples", true) );
  glarea.addParam(new RichBool("Show Dual Samples Quad", false));
  glarea.addParam(new RichBool("Show Dual Samples Dot", true));
  glarea.addParam(new RichBool("Show Dual Samples Circle", false));
  glarea.addParam(new RichBool("Show Dual Samples Sphere", false));

  glarea.addParam(new RichBool("Show Dual Connection", false) );
  glarea.addParam(new RichBool("Pick Dual Point", false) );
	glarea.addParam(new RichBool("Show Skeltal Points", true));

	glarea.addParam(new RichBool("Show Original", false) );
	glarea.addParam(new RichBool("Show Original Quad", false) );
	glarea.addParam(new RichBool("Show Original Dot", true) );
	glarea.addParam(new RichBool("Show Original Circle", false) );
	glarea.addParam(new RichBool("Show Original Sphere", false) );


	glarea.addParam(new RichBool("Show Target Samples", true));
	glarea.addParam(new RichBool("Show Target Dual Samples", true));
	glarea.addParam(new RichBool("Show Correspondences", false));
	glarea.addParam(new RichBool("Show Eigen Directions", false));

	glarea.addParam(new RichBool("Show Bounding Box", false));
	glarea.addParam(new RichBool("Show Skeleton", false));

	glarea.addParam(new RichBool("Show Radius", true));
	glarea.addParam(new RichBool("Show All Radius", false));
	glarea.addParam(new RichBool("Show Radius Use Pick", true));
	glarea.addParam(new RichBool("Show Red Radius Line", true));
	glarea.addParam(new RichBool("Multiply Pick Point", false) );
	glarea.addParam(new RichBool("Random Erase", false));


  glarea.addParam(new RichDouble("Picked Index", 0));
	glarea.addParam(new RichBool("Show Cloest Dual Connection", true));

	glarea.addParam(new RichBool("GLarea Busying", false) );
  glarea.addParam(new RichBool("Algorithom Stop", false) );


	glarea.addParam(new RichPoint3f("Light Position", vcg::Point3f(-4.0, -4.0, -4.0)));
	glarea.addParam(new RichColor("Light Ambient Color", QColor(44, 44, 44)));
	//glarea.addParam(new RichColor("Light Diffuse Color", QColor(164, 241, 101)));
	glarea.addParam(new RichColor("Light Diffuse Color", QColor(200, 200, 200)));
	glarea.addParam(new RichColor("Light Specular Color", QColor(255, 255, 255)));

	//glarea.addParam(new RichPoint3f("Light Position", vcg::Point3f(4.0, 4.0, 4.0)));
	//glarea.addParam(new RichColor("Light Ambient Color", QColor(0.0, 0.0, 0.0)));
	//glarea.addParam(new RichColor("Light Diffuse Color", QColor(204, 204, 204)));
	//glarea.addParam(new RichColor("Light Specular Color", QColor(255, 255, 255)));

	glarea.addParam(new RichDouble("Snapshot Resolution", 2));
	glarea.addParam(new RichDouble("Snapshot Index", 1));
	glarea.addParam(new RichDouble("Radius Ball Transparency", 0.3));

	glarea.addParam(new RichBool("SnapShot Each Iteration", false));
	glarea.addParam(new RichBool("No Snap Radius", true));

	glarea.addParam(new RichBool("Need Snap Files", false));

  glarea.addParam(new RichBool("Fix ply mode", false));


	glarea.addParam(new RichDouble("Sample Confidence Color Scale", 0.5));
	glarea.addParam(new RichDouble("Point ISO Value Shift", -0.5));

	glarea.addParam(new RichDouble("Show Confidence Percentage", 0.8));
}


void ParameterMgr::initDrawerParameter()
{
	drawer.addParam(new RichBool("Doing Pick", false));
	drawer.addParam(new RichBool("Need Cull Points", false) );
	drawer.addParam(new RichBool("Use Pick Original", false));
	drawer.addParam(new RichBool("Use Pick Mode2", false) );
	drawer.addParam(new RichBool("Skeleton Light", true));
	drawer.addParam(new RichBool("Show Individual Color", false));

	drawer.addParam(new RichBool("Use Color From Normal", false));
  drawer.addParam(new RichBool("Use Differ Branch Color", false));
	drawer.addParam(new RichBool("Show Confidence Color", false));
	drawer.addParam(new RichBool("Show Segmentation Color", false));

	drawer.addParam(new RichBool("Show Feature Color", true));
  drawer.addParam(new RichBool("Draw Picked Point Neighbor", false));


	drawer.addParam(new RichDouble("Original Draw Width", 0.0010));
	drawer.addParam(new RichDouble("Sample Draw Width", 0.01));
	drawer.addParam(new RichDouble("Dual Sample Draw Width", 0.010));

	drawer.addParam(new RichDouble("Sample Dot Size", 6));
	drawer.addParam(new RichDouble("Dual Sample Dot Size", 6));
	drawer.addParam(new RichDouble("Original Dot Size", 1));
	drawer.addParam(new RichDouble("Normal Line Width", 2));
	drawer.addParam(new RichDouble("Normal Line Length", 0.06));

	drawer.addParam(new RichColor("Background Color", QColor(255, 255, 255) ));
	drawer.addParam(new RichColor("Normal Line Color", QColor(0, 0, 255) ));
	drawer.addParam(new RichColor("Sample Point Color", QColor(255, 0, 0) ));
	drawer.addParam(new RichColor("Original Point Color", QColor(22, 22, 22) ));
	drawer.addParam(new RichColor("Feature Color", QColor(0, 0, 255) ));
	drawer.addParam(new RichColor("Pick Point Color", QColor(128, 128, 0) ));
	drawer.addParam(new RichColor("Pick Point DNN Color", QColor(0, 0, 155) ));
	drawer.addParam(new RichColor("DLink Color", QColor(0, 0, 0)));
	drawer.addParam(new RichColor("Backface Color", QColor(175, 175, 175)));


	drawer.addParam(new RichColor("Skeleton Bone Color", QColor(200, 0, 0) ));
	drawer.addParam(new RichColor("Skeleton Node Color", QColor(50, 250, 50) ));
	drawer.addParam(new RichColor("Skeleton Branch Color", QColor(0, 0, 0)));
	drawer.addParam(new RichDouble("Skeleton Bone Width", 100)); // ./10000
	drawer.addParam(new RichDouble("Skeleton Node Size", 180)); // ./10000
	drawer.addParam(new RichDouble("Skeleton Branch Size", 30)); // abandoned
}


void ParameterMgr::initWLopParameter()
{
	wLop.addParam(new RichString("Algorithm Name", "WLOP") );
	wLop.addParam(new RichDouble("Num Of Iterate Time", 1));

	wLop.addParam(new RichDouble("CGrid Radius", grid_r));
	wLop.addParam(new RichDouble("H Gaussian Para", 4));
	wLop.addParam(new RichDouble("Repulsion Power", 1.0));
	wLop.addParam(new RichDouble("Average Power", 2.0));
	wLop.addParam(new RichBool("Need Compute Density", false));
	wLop.addParam(new RichBool("Need Compute PCA", false));
	wLop.addParam(new RichDouble("Repulsion Mu", 0.500));
	wLop.addParam(new RichDouble("Repulsion Mu2", 0.10));
	wLop.addParam(new RichDouble("Dual Radius", 0.33));

	wLop.addParam(new RichDouble("Local Angle Threshold", 45));

  wLop.addParam(new RichDouble("Dual Mu3", 0.35));
  wLop.addParam(new RichBool("Need Sample Average", false));
	wLop.addParam(new RichBool("Use Confidence", false));

	wLop.addParam(new RichBool("Run Anisotropic LOP", false));
  wLop.addParam(new RichBool("Run Dual WLOP", false));
  wLop.addParam(new RichBool("Run Step Forward", false));

  wLop.addParam(new RichBool("Run Skel WLOP", false));
  wLop.addParam(new RichBool("Run Dual Drag WLOP", false));
  wLop.addParam(new RichBool("Run Regularize Samples", false));
  wLop.addParam(new RichBool("Run Regularize Normals", false));

  wLop.addParam(new RichBool("Run Projection", false));
	wLop.addParam(new RichBool("Run Compute Confidence", false));
	wLop.addParam(new RichBool("Run Compute Distribution", false));
	wLop.addParam(new RichBool("Run Inner Clustering", false));
	wLop.addParam(new RichBool("Run Show Pick Distribution", false));
	wLop.addParam(new RichBool("Run Compute Correspondence", false));
	wLop.addParam(new RichBool("Run Detect Kite Points", false));
	wLop.addParam(new RichBool("Run Inner Points Classification", false));
	wLop.addParam(new RichBool("Run Compute Dual Index", false));


	wLop.addParam(new RichBool("WLOP test bool", false));


	wLop.addParam(new RichBool("Run Search Neighborhood", false));
	wLop.addParam(new RichBool("Run Smooth Neighborhood", false));
	wLop.addParam(new RichBool("Run Inner Points Regularization", false));

	wLop.addParam(new RichBool("Run Move Backward", false));
	wLop.addParam(new RichBool("Run Self WLOP", false));
	wLop.addParam(new RichBool("Run Normal Smoothing", false));

	wLop.addParam(new RichBool("Run Self PCA", false));
	wLop.addParam(new RichBool("Run Self Projection", false));

  wLop.addParam(new RichBool("Run 4PCS", false));

  wLop.addParam(new RichBool("Original Combine Sample", false));
  wLop.addParam(new RichBool("Need Averaging Movement", false));
	wLop.addParam(new RichBool("Need Similarity", false));


  wLop.addParam(new RichBool("Use Elliptical Original Neighbor", false));
  wLop.addParam(new RichBool("Use Adaptive Sample Neighbor", false));
  wLop.addParam(new RichBool("Use Adaptive Mu", false));
  wLop.addParam(new RichBool("Use Tangent Vector", false));
	
	wLop.addParam(new RichBool("Use Original Averaging KNN", false));
	wLop.addParam(new RichBool("Use Kite Points", false));

	wLop.addParam(new RichDouble("Sefl KNN", 15));


	wLop.addParam(new RichDouble("Eigen Neighborhood Para1", 0.5));
	wLop.addParam(new RichDouble("Eigen Neighborhood Para2", 3.0));

	wLop.addParam(new RichBool("Use Eigen Neighborhood", true));

	wLop.addParam(new RichBool("Use Ellipsoid Weight", false));
	wLop.addParam(new RichBool("Use Ellipsoid Repulsion", true));
	wLop.addParam(new RichBool("Use Separate Neighborhood", true));


  wLop.addParam(new RichBool("Run Compute Initial Sample Neighbor", false));
  wLop.addParam(new RichBool("Run Normal Reform", false));

	wLop.addParam(new RichDouble("Current Movement Error", 0.0));
	

	wLop.addParam(new RichDouble("Original Averaging KNN", 415));

	wLop.addParam(new RichDouble("sigmoid threshold", 0.4));

	wLop.addParam(new RichDouble("Density Confidence Threshold", 0.911));
	wLop.addParam(new RichDouble("Eigen Confidence Threshold", 0.97));
	wLop.addParam(new RichDouble("Mutual Distance Threshold", 0.1));


	wLop.addParam(new RichDouble("Progressive Min KNN", 20));
	wLop.addParam(new RichDouble("Progressive Max KNN", 2000));
	wLop.addParam(new RichBool("Run Progressive Neighborhood", false));

	wLop.addParam(new RichBool("Run Ellipsoid Fitting", false));
	wLop.addParam(new RichBool("Run MAT LOP", false));

	wLop.addParam(new RichBool("Run Compute Initial Neighborhood", false));

	wLop.addParam(new RichBool("Run Move Sample", false));
	wLop.addParam(new RichBool("Run Move Skel", false));

	wLop.addParam(new RichBool("Compute Eigen Directions", false));
	wLop.addParam(new RichBool("Compute Eigen Neighborhood", false));


	// generating inner points
	wLop.addParam(new RichDouble("Increasing Step Size", 0.01));
	wLop.addParam(new RichDouble("Local Neighbor Size For Inner Points", 0.040));
	wLop.addParam(new RichDouble("Local Neighbor Size For Surface Points", 0.1));
	wLop.addParam(new RichDouble("Inner Points Cooling Parameter", 0.5));


	//new parameters:
	wLop.addParam(new RichDouble("Average Closest Dist", 0.02301));
	wLop.addParam(new RichDouble("Search Dual Index Para", 2.1));
	wLop.addParam(new RichBool("Use Progressive Search Index", true));


	//void WLOP::runComputeConfidence()
	wLop.addParam(new RichBool("Run Estimate Average Dist Threshold", false));
	wLop.addParam(new RichBool("Use Average Dist Threshold", false));


	wLop.addParam(new RichDouble("Average Dist To Input Threshold", 0.11889));
	wLop.addParam(new RichDouble("Choose ADT Threshold Percentage", 0.50000));
	wLop.addParam(new RichDouble("Original Confidence KNN", 155));

	wLop.addParam(new RichDouble("Confidence Power", 3.0));

	//void WLOP::computeSampleSimilarityTerm(CMesh* samples)
	wLop.addParam(new RichDouble("Similarity Term Neighbor Para", 1.50000));
	wLop.addParam(new RichDouble("Similarity Length Outlier Threshold", 1.22));//no need?
	wLop.addParam(new RichDouble("KNN For Similarity", 15));//no need?


	//void WLOP::runRegularizeNormals()
	wLop.addParam(new RichDouble("Density Confidence Segment Threshold", 0.8));
	wLop.addParam(new RichDouble("Eigen Directional Threshold", 0.55000));
	wLop.addParam(new RichBool("Use Confidence To Combine Normal", false));


	//vector<Point3f> WLOP::computeNewSamplePositions(int& error_x)
	wLop.addParam(new RichDouble("Save Move Dist Along Normal Para", 0.55));
	wLop.addParam(new RichDouble("Big Repulsion Power", 5.0));
	//
	wLop.addParam(new RichDouble("Protect Small Tubular Para", 0.5));
	wLop.addParam(new RichDouble("Protect High Confidence Para", 0.90));

	wLop.addParam(new RichDouble("Data Outweigh Similarity Para", 1.0));

	wLop.addParam(new RichBool("Only Do Repuslion", false));
	wLop.addParam(new RichBool("Only Do Avergage", false));
	wLop.addParam(new RichBool("Use Confidence To Merge", false));

	//2015-2-2
	wLop.addParam(new RichBool("Dual Samples Represent Inner Points", false));
	wLop.addParam(new RichBool("Dual Samples Represent Skeltal Points", false));

	wLop.addParam(new RichBool("Run Tangential Motion", false));
	wLop.addParam(new RichBool("Run DLength Adjustment", false));

}

void ParameterMgr::initSkeletonParameter()
{
  /// 
  skeleton.addParam(new RichDouble("Repulsion Power", 1.0));
  skeleton.addParam(new RichDouble("Average Power", 2.0));

  /// 
  skeleton.addParam(new RichDouble("Num Of Iterate Time", 1));
	skeleton.addParam(new RichString("Algorithm Name", "Skeletonization") );
	

	skeleton.addParam(new RichDouble("CGrid Radius", grid_r));
	skeleton.addParam(new RichDouble("H Gaussian Para", 4));
	skeleton.addParam(new RichBool("Need Compute Density", true));
	
	
	skeleton.addParam(new RichDouble("Current Movement Error", 0.0));
	skeleton.addParam(new RichBool("Run Auto Wlop One Step", false));
	skeleton.addParam(new RichBool("Run Auto Wlop One Stage", false));
	skeleton.addParam(new RichBool("The Skeletonlization Process Should Stop", false));

	skeleton.addParam(new RichBool("Step1 Detect Skeleton Feature", false));
	skeleton.addParam(new RichBool("Step2 Run Search New Branchs", false));
	skeleton.addParam(new RichBool("Step3 Clean And Update Radius", false));
  skeleton.addParam(new RichBool("Run Skeletonlization", false));



	//init
	skeleton.addParam(new RichDouble("Max Iterate Time", 55));
	skeleton.addParam(new RichDouble("Stop And Grow Error", 0.0001));
	skeleton.addParam(new RichDouble("Initial Radius", -1.));
	skeleton.addParam(new RichDouble("Radius Update Speed", 0.5));

	//step0
	skeleton.addParam(new RichDouble("Repulsion Mu", 0.5));
	skeleton.addParam(new RichDouble("Repulsion Mu2", 0.15));
	skeleton.addParam(new RichDouble("Follow Sample Radius", 0.33));
	skeleton.addParam(new RichDouble("Follow Sample Max Angle", 80));// should add to UI
	skeleton.addParam(new RichDouble("Inactive And Keep Virtual Angle", 60)); // should add to UI
	skeleton.addParam(new RichDouble("Save Virtual Angle", 30)); // should add to UI

	skeleton.addParam(new RichDouble("Grow Accept Sigma", 0.8));// should add to UI
	skeleton.addParam(new RichDouble("Bad Virtual Angle", 101));// 2013-7-12

	//step1
	skeleton.addParam(new RichDouble("Combine Too Close Threshold", 0.01));
	skeleton.addParam(new RichDouble("Sigma KNN", 16));//this one is hard to determine, should be small for narrow region, but will lead to unnecessary small branches
	skeleton.addParam(new RichDouble("Eigen Feature Identification Threshold", 0.95));

	//step2
	skeleton.addParam(new RichDouble("Branches Search Angle", 25));
	skeleton.addParam(new RichDouble("Virtual Head Accecpt Angle", 25));
	skeleton.addParam(new RichDouble("Snake Search Max Dist Blue", 0.4));
	skeleton.addParam(new RichDouble("Accept Branch Size", 6)); // important, and hard to determine
	skeleton.addParam(new RichDouble("Branch Search Max Dist Yellow", 0.1));

	skeleton.addParam(new RichDouble("Branches Merge Max Dist", 0.08));
	skeleton.addParam(new RichDouble("Branch Search KNN", 12));
	skeleton.addParam(new RichDouble("Combine Similar Angle", 140));
	skeleton.addParam(new RichDouble("Grow Search Radius", 0.15));
	skeleton.addParam(new RichDouble("Add Accept Branch Size", 1));


	//step3
	skeleton.addParam(new RichDouble("Clean Near Branches Dist", 0.05));
	skeleton.addParam(new RichDouble("Fix Original Weight", 0.81));
	skeleton.addParam(new RichDouble("Curve Segment Length", 0.051));
	skeleton.addParam(new RichInt("Fix Original Mode", 4)); // 1 for noisy , 4 for clean

  skeleton.addParam(new RichBool("Run ALL Segment", false));
  skeleton.addParam(new RichBool("Need Segment Right Away", true));
  skeleton.addParam(new RichDouble("Max Stop Radius", 1.99));

	//strategy...
	skeleton.addParam(new RichBool("Use Nearby Combine Strategy", true));
	skeleton.addParam(new RichBool("Use Go Through Strategy", false));
	skeleton.addParam(new RichBool("Use Aggresive Growth Strategy", false));
	skeleton.addParam(new RichBool("Use Clean Points When Following Strategy", true));
	skeleton.addParam(new RichBool("Use All Connect Strategy", true));
	skeleton.addParam(new RichBool("Use Plus Perpendicular Dist Strategy", false));
	skeleton.addParam(new RichBool("Use Kill Too Close Strategy", false));
	skeleton.addParam(new RichBool("Use Compute Eigen Ignore Branch Strategy", true));
	skeleton.addParam(new RichBool("Use Virtual Group Merge Strategy", false));
	skeleton.addParam(new RichBool("Use Final Merge Strategy", true));
	skeleton.addParam(new RichBool("Use Search New Twice Strategy", false));
	skeleton.addParam(new RichBool("Inactive Overlap Strategy", false));
	skeleton.addParam(new RichBool("Move Overlap Strategy", false));
	skeleton.addParam(new RichBool("Use Virtual Near Body Stop Strategy", true));
	skeleton.addParam(new RichBool("Need To Keep Big Bug", false));
	skeleton.addParam(new RichDouble("Change Strategy Radius", 0.45));
	skeleton.addParam(new RichBool("Need Recentering", true));
}

void ParameterMgr::initNormalSmootherParameter()
{
	norSmooth.addParam(new RichString("Algorithm Name", "NormalSmooth") );

	norSmooth.addParam(new RichInt("PCA KNN", 9));
	norSmooth.addParam(new RichDouble("CGrid Radius", grid_r));
	norSmooth.addParam(new RichDouble("Sharpe Feature Bandwidth Sigma", 90));
	norSmooth.addParam(new RichBool("Run Anistropic PCA", false));
	norSmooth.addParam(new RichBool("Run Init Samples Using Normal", false));

	norSmooth.addParam(new RichInt("Number Of Iterate", 1));
	norSmooth.addParam(new RichInt("Number of KNN", 400));

	norSmooth.addParam(new RichDouble("PCA Threshold", 0.8));
}


void ParameterMgr::initUpsamplingParameter()
{
	upsampling.addParam(new RichString("Algorithm Name", "Upsampling") );

	upsampling.addParam(new RichDouble("CGrid Radius", 0.08) );
	upsampling.addParam(new RichInt("Number of Add Point", 5000) );
	upsampling.addParam(new RichDouble("Feature Sigma", 90));

	upsampling.addParam(new RichBool("Using Threshold Process", true) );
	upsampling.addParam(new RichDouble("Dist Threshold", 0.02));
	upsampling.addParam(new RichDouble("Edge Parameter", 0.0));
	upsampling.addParam(new RichDouble("Z Parameter", 0.1));

  upsampling.addParam(new RichBool("Use Constant Threshold", true) );
  upsampling.addParam(new RichBool("Run Predict Constant Threshold", false) );


	upsampling.addParam(new RichBool("Auto Recompute Radius For Dist", true) );
	upsampling.addParam(new RichDouble("Min Dist Rate", 2.0));

	upsampling.addParam(new RichDouble("New Point Avg Sigma", 15));
	upsampling.addParam(new RichBool("Use Avg Normal Method", false) );
	upsampling.addParam(new RichBool("Use Max Theta Psi Method", false));
	upsampling.addParam(new RichBool("Use Sigma Threshold Method", true));
	upsampling.addParam(new RichBool("Use No Psi Method", false));
	upsampling.addParam(new RichBool("Use Upsample For Closeby Surface", true));

	upsampling.addParam(new RichDouble("Upsample Radius", grid_r * 0.5) );
	upsampling.addParam(new RichBool("Use Proj New Term", false));
	upsampling.addParam(new RichBool("Run Projection", false));

  upsampling.addParam(new RichBool("Run Points Extrapolation", false));

  upsampling.addParam(new RichBool("Use Adaptive Upsampling", false));
  upsampling.addParam(new RichDouble("Density Threshold Dist", 0.05));


}


void ParameterMgr::outputAllParameters(std::ostream& out)
{
	outputParameters(out, glarea);
	outputParameters(out, data);
	outputParameters(out, drawer);
	outputParameters(out, wLop);
	outputParameters(out, norSmooth);
	outputParameters(out, skeleton);
	outputParameters(out, upsampling);
}

void ParameterMgr::outputParameters(std::ostream& out, RichParameterSet& para_set)
{
	QList<RichParameter*> paramList = para_set.paramList;
	QList<RichParameter*>::iterator iter;
	for (iter = paramList.begin(); iter != paramList.end(); ++iter)
	{
		RichParameter* para = *iter;
		QString name = para->name;
		if (para->val->isString())
		{
			continue;
		}
		//out << '"' << name.toStdString().c_str() << '"' << "|";
		out << name.toStdString().c_str() << "|";

		if (para->val->isBool())
		{
			out << "BoolValue|" << para->val->getBool() << std::endl;
		}
		else if (para->val->isDouble())
		{
			out << "DoubleValue|" << para->val->getDouble() << std::endl;
		}
		else if (para->val->isInt())
		{
			out << "IntValue|" << para->val->getInt() << std::endl;
		}
		else if (para->val->isColor())
		{
			out << "ColorValue|" << para->val->getColor().red() << ","
				                   << para->val->getColor().green() << "," 
				                   << para->val->getColor().blue()  << std::endl;
		}
		else if (para->val->isPoint3f())
		{
			out << "Point3fValue|" << para->val->getPoint3f().X() << ","
				                     << para->val->getPoint3f().Y()  << ","
				                     << para->val->getPoint3f().Z() << std::endl;
		}
		else
		{
			std::cout << "unexpected type!!! " << std::endl;
		}
		//out << std::endl;
		//out << para->val-> << endl;
	}
}


void ParameterMgr::inputAllParameters(std::istream& in)
{
	QString qstr;
	std::string str;
	while (!in.eof())
	{
		//in.getline(str, 200);
		std::getline(in, str, '\n');
		qstr = str.c_str();
		//std::cout << str.c_str() << std::endl << std::endl;

		QStringList qlist = qstr.split("|");
		QStringList::iterator iter;
		if (qlist.size() != 3)
		{
			std::cout << " para file wrong formate! " << std::endl;
			std::cout << str.c_str() << std::endl << std::endl;

			continue;
		}

		QString para_name = qlist[0];
		QString para_type = qlist[1];
		QString para_value = qlist[2];

		if (para_type == QString("BoolValue"))
		{
			bool b = para_value.toInt();
			setGlobalParameter(para_name, BoolValue(b));
		}
		else if (para_type == QString("DoubleValue"))
		{
			double tmp = para_value.toDouble();
			setGlobalParameter(para_name, DoubleValue(tmp));
		}
		else if (para_type == QString("IntValue"))
		{
			int tmp = para_value.toInt();
			setGlobalParameter(para_name, IntValue(tmp));
		}
		else if (para_type == QString("ColorValue"))
		{
			QStringList color_list = para_value.split(",");
			if (color_list.size() != 3)
			{
				std::cout << "wrong color!" << std::endl;
				continue;
			}
			int r, g, b;
			r = color_list[0].toInt();
			g = color_list[1].toInt();
			b = color_list[2].toInt();
			QColor color;
			color.setRgb(r, g, b);
			setGlobalParameter(para_name, ColorValue(color));
		}
		else if (para_type == QString("Point3fValue"))
		{
			QStringList float_list = para_value.split(",");
			if (float_list.size() != 3)
			{
				std::cout << "wrong point3f!" << std::endl;
				continue;
			}
			float x, y, z;
			x = float_list[0].toFloat();
			y = float_list[1].toFloat();
			z = float_list[2].toFloat();
			Point3f pt(x, y, z);
			
			setGlobalParameter(para_name, Point3fValue(pt));
		}
		else
		{
			std::cout << "type: " << para_type.toStdString().c_str() << std::endl;
		}

	}
}
