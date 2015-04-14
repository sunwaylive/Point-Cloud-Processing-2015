#include "WLOP.h"


WLOP::WLOP(RichParameterSet* _para)
{
	cout << "WLP constructed!!" << endl;
	para = _para;
	samples = NULL;
	original = NULL;
	dual_samples = NULL;
	skel_points = NULL;

	nTimeIterated = 0;
	error_x = 0.0;
}

WLOP::~WLOP(void)
{
	cout << "WLop destroy!! " << endl; 
}

void WLOP::clear()
{
	samples = NULL;
	original = NULL;
	dual_samples = NULL;
	skel_points = NULL;

}

void WLOP::setFirstIterate()
{
	nTimeIterated = 0;
}

void WLOP::setInput(DataMgr* pData)
{
	use_closest_dual = global_paraMgr.glarea.getBool("Show Cloest Dual Connection");
	use_kite_points = para->getBool("Use Kite Points");
	use_eigen_neighborhood = para->getBool("Use Eigen Neighborhood");
	use_ellipsoid_weight = para->getBool("Use Ellipsoid Weight");
	use_ellipsoid_repulsion = para->getBool("Use Ellipsoid Repulsion");

	repulsion_radius = para->getDouble("CGrid Radius");
	if (para->getBool("Use Separate Neighborhood"))
	{
		repulsion_radius = para->getDouble("CGrid Radius") * para->getDouble("Eigen Neighborhood Para1");
	}

	if(!pData->isSamplesEmpty() && !pData->isOriginalEmpty())
	{
		default_sphere = pData->default_sphere;
		CMesh* _samples = pData->getCurrentSamples();
		CMesh* _original = pData->getCurrentOriginal();

		target_samples = pData->getCurrentTargetSamples();
		target_dual_samples = pData->getCurrentTargetDualSamples();

		if (para->getBool("Run Skel WLOP"))
		{
			_samples = pData->getCurrentSkelPoints();
			_original = pData->getCurrentDualSamples();
		}

		if(_samples == NULL || _original == NULL)
		{
			cout << "ERROR: WLOP::setInput == NULL!!" << endl;
			return;
		}

		error_x = 0.0;
		samples = _samples;
		original = _original;

    dual_samples = pData->getCurrentDualSamples();
		skel_points = pData->getCurrentSkelPoints();

		if (!para->getBool("Run Skel WLOP") && para->getBool("Dual Samples Represent Skeltal Points"))
		{
			cout << "Dual Samples Represent Skeltal Points" << endl;
			dual_samples = pData->getCurrentSkelPoints();
		}
		else if (!para->getBool("Run Skel WLOP") && para->getBool("Dual Samples Represent Inner Points"))
		{
			cout << "Dual Samples Represent Inner Points" << endl;

			dual_samples = pData->getCurrentDualSamples();
		}
		else if (!para->getBool("Run Skel WLOP"))
		{
			cout << "Be Carefull! Dual Samples Representation is not clear!!" << endl;
			system("Pause");
		}

    if (para->getBool("Run Dual WLOP"))
    {
      samples = pData->getCurrentDualSamples();

			if (global_paraMgr.glarea.getBool("Show Samples"))
			{
				original = pData->getCurrentSamples();
			}
    }
		samples_density.assign(samples->vn, 1);
	}
	else
	{
		cout << "ERROR: WLOP::setInput: empty!!" << endl;
		return;
	}
}


void WLOP::initVertexes(bool clear_neighbor)
{
  bool use_current_neighbor = para->getBool("Use Adaptive Sample Neighbor");
	
  box.SetNull();
	CMesh::VertexIterator vi, vi_end;

	int i = 0;
	vi_end = samples->vert.end();
	for(vi = samples->vert.begin(); vi != vi_end; ++vi) 
	{
		vi->m_index = i++;
    if (!use_current_neighbor && clear_neighbor)
    {
      vi->neighbors.clear();
      vi->original_neighbors.clear();
    }
		
		if (vi->is_skel_ignore)
		{
			continue;
		}
		box.Add(vi->P());
	}
	samples->bbox = box;

	i = 0;
	vi_end = dual_samples->vert.end();
	for (vi = dual_samples->vert.begin(); vi != vi_end; ++vi)
	{
		vi->m_index = i++;
		if (!use_current_neighbor && clear_neighbor)
		{
			vi->neighbors.clear();
			vi->original_neighbors.clear();
		}

		if (vi->is_skel_ignore)
		{
			continue;
		}
		box.Add(vi->P());
	}
	dual_samples->bbox = box;



	vi_end = original->vert.end();
	i = 0;
	for(vi = original->vert.begin(); vi != vi_end; ++vi) 
	{
		vi->m_index = i++;
		box.Add(vi->P());
	}
	original->bbox = box;

	repulsion.assign(samples->vn, vcg::Point3f(0, 0, 0));
	repulsion_x_length.assign(samples->vn, 0.);
	repulsion_y_length.assign(samples->vn, 0.);
	repulsion_z_length.assign(samples->vn, 0.);

	//Eigen::Matrix3f iden = Eigen::Matrix3d::Identity(3,3);
	Matrix33f iden;
	iden.SetIdentity();

	repulsion_matA_set.assign(samples->vn, iden);

	average.assign(samples->vn, vcg::Point3f(0, 0, 0));

	repulsion_weight_sum.assign(samples->vn, 0);
	average_weight_sum.assign(samples->vn, 0);

  samples_average.assign(samples->vn, vcg::Point3f(0, 0, 0));
  samples_average_weight_sum.assign(samples->vn, 0);

	average_low_confidence.assign(samples->vn, vcg::Point3f(0, 0, 0));
	average_weight_sum_low_confidence.assign(samples->vn, 0);

	samples_similarity.assign(samples->vn, vcg::Point3f(0, 0, 0));
	//samples_similarity_weight_sum.assign(samples->vn, 0);

	if (para->getBool("Need Compute PCA"))
	{
		CVertex v;
		mesh_temp.assign(samples->vn, v);
	}
}


void WLOP::run()
{


  if (para->getBool("Run Projection"))
  {
    runProjection();
    return;
  }


	if (para->getBool("Run Anisotropic LOP"))
	{
		cout << "Run Anisotropic LOP" << endl;
	}

  if (para->getBool("Run Step Forward"))
  {
    cout << "Run Step Forward" << endl;
		smoothSkelDistance();
    return;
  }

  if (para->getBool("Run Compute Initial Sample Neighbor"))
  {
    //computeInitialSampleNeighbor();
		computeNearestNeighborDist();
    return;
  }



	if (para->getBool("Run MAT LOP"))
	{
		cout << "Run MAT LOP removed" << endl;
		//runMatLOP();
		return;
	}

	if (para->getBool("Run Detect Kite Points"))
	{
		cout << "Run Detect Kite Points removed" << endl;
		//runDetectKitePoitns();
		return;
	}

  if (para->getBool("Run Dual Drag WLOP"))
  {
		cout << "Run Dual Drag WLOP removed" << endl;
    return;
  }


  if (para->getBool("Run Regularize Samples"))
  {
    runRegularizeSamples();
    return;
  }

  if (para->getBool("Run Regularize Normals"))
  {
    runRegularizeNormals(samples, skel_points);
    return;
  }

	if (para->getBool("Run Compute Confidence"))
	{
		runComputeConfidence();
		//runComputeHoleConfidence();
		return;
	}

	if (para->getBool("Run Compute Distribution"))
	{
		runComputeDistribution();
		return;
	}

	if (para->getBool("Run Inner Clustering"))
	{
		//runComputeInnerClusering();
		return;
	}



	if (para->getBool("Run Show Pick Distribution"))
	{
		//runShowPickDistribution();
		return;
	}

	if (para->getBool("Run Compute Correspondence"))
	{
		runComputeCorrespondence();
		return;
	}

	if (para->getBool("Run Progressive Neighborhood"))
	{
		//runProgressiveNeighborhood();
		return;
	}

	if (para->getBool("Run Inner Points Classification"))
	{
		//innerpointsClassification();
		return;
	}

	if (para->getBool("Run Ellipsoid Fitting"))
	{
		runEllipsoidFitting();
		return;
	}


	if (para->getBool("Run Search Neighborhood"))
	{
		runSearchNeighborhood();
		return;
	}



	if (para->getBool("Run Smooth Neighborhood"))
	{
		runComputeAverageDistThreshold();
		//runSmoothNeighborhood();
		return;
	}

	if (para->getBool("Run Inner Points Regularization"))
	{
		runInnerPointsRegularization();
		return;
	}


	if (para->getBool("Run Move Backward"))
	{
		runMoveBackward();
		return;
	}

	if (para->getBool("Run Normal Smoothing"))
	{
		runNormalSmoothing();
		return;
	}

	
	if (para->getBool("Run Self WLOP"))
	{
		runSelfWLOP();

		nTimeIterated++;
		cout << "Iterated: " << nTimeIterated << endl;
		return;
	}

	if (para->getBool("Run Self PCA"))
	{
		runSelfPCA();
		return;
	}

	if (para->getBool("Run Self Projection"))
	{
		runSelfProjection();
		return;
	}

	if (para->getBool("Run Compute Initial Neighborhood"))
	{
		computeInitialNeighborSize();
		return;
	}

	if (para->getBool("Run Move Sample"))
	{
		runMoveSample();
		return;
	}

	if (para->getBool("Run Move Skel"))
	{
		runMoeveSkel();
		return;
	}


	if (para->getBool("Compute Eigen Directions"))
	{
		if (para->getBool("Run Skel WLOP"))
		{
			runComputeEigenDirections(samples, original);
		}
		else
		{
			runComputeEigenDirections(dual_samples, samples);
		}
		return;
	}

	if (para->getBool("Compute Eigen Neighborhood"))
	{
		runComputeEigenNeighborhood(dual_samples, samples);
		return;
	}

	if (para->getBool("Run Skel WLOP"))
	{
		runSkelWlop();
		//runComputeEigenNeighborhood(samples, original);
		return;
	}


	if (para->getBool("Run Compute Dual Index"))
	{
		computeDualIndex(samples, dual_samples);
		//runComputeEigenNeighborhood(samples, original);
		return;
	}

	if (para->getBool("Run Tangential Motion"))
	{
		runTangentialMotion();
		return;
	}

 	if (para->getBool("Run DLength Adjustment"))
 	{
 		runDlengthAdjustment();
 		return;
 	}

	if (para->getBool("Run Estimate Average Dist Threshold"))
	{
		//runComputeAverageDistThreshold();
		//runComputeConfidence();
    runEstimateParameters();
		return;
	}

  if (para->getBool("Run 4PCS"))
  {
    run4PCS();
    return;
  }


	//int nTimes = para->getDouble("Num Of Iterate Time");
	for(int i = 0; i < 1; i++)
	{ 
		cout << "run iterate" << endl;
		iterate();
		
		nTimeIterated ++;
		cout << "Iterated: " << nTimeIterated << endl;
	}
	cout << "**************iterate Number: " << nTimeIterated << endl;
}

void WLOP::computeAverageAddSampleTerm(CMesh* samples, CMesh* original)
{
  double average_power = para->getDouble("Average Power");
  bool need_density = para->getBool("Need Compute Density");
  double radius = para->getDouble("CGrid Radius"); 

  double radius2 = radius * radius;
  double iradius16 = -para->getDouble("H Gaussian Para")/radius2;

	bool use_tangent = para->getBool("Use Tangent Vector");

  for(int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];

    for (int j = 0; j < v.original_neighbors.size(); j++)
    {
      CVertex& t = original->vert[v.original_neighbors[j]];

      Point3f diff = v.P() - t.P();
			double proj_dist = diff * v.N();
			double proj_dist2 = proj_dist * proj_dist;

      double dist2  = diff.SquaredNorm();

      double w = 1;
      if (para->getBool("Run Anisotropic LOP"))
      {
        double len = sqrt(dist2);
        if(len <= 0.001 * radius) len = radius*0.001;
        double hn = diff * v.N();
        double phi = exp(hn * hn * iradius16);
        w = phi / pow(len, 2 - average_power);
      }
      else if (average_power < 2)
      {
        double len = sqrt(dist2);
        if(len <= 0.001 * radius) len = radius*0.001;
        w = exp(dist2 * iradius16) / pow(len, 2 - average_power);
      }
      else
      {
        w = exp(dist2 * iradius16);
      }

      if (need_density)
      {
        w *= original_density[t.m_index];
      }

			if (use_tangent)
			{
				Point3f proj_point = v.P() + v.N() * proj_dist;
				average[i] += proj_point * w;
			}
			else
			{
				average[i] += t.P() * w;
			}
      average_weight_sum[i] += w;  
    }

    for (int k = 0; k < v.neighbors.size(); k++)
    {
      CVertex& t = samples->vert[v.neighbors[k]];

      Point3f diff = v.P() - t.P();
			double proj_dist = diff * v.N();
			double proj_dist2 = proj_dist * proj_dist;

      double dist2  = diff.SquaredNorm();

      double w = 1;
      if (para->getBool("Run Anisotropic LOP"))
      {
        double len = sqrt(dist2);
        if(len <= 0.001 * radius) len = radius*0.001;
        double hn = diff * v.N();
        double phi = exp(hn * hn * iradius16);
        w = phi / pow(len, 2 - average_power);
      }
      else if (average_power < 2)
      {
        double len = sqrt(dist2);
        if(len <= 0.001 * radius) len = radius*0.001;
        w = exp(dist2 * iradius16) / pow(len, 2 - average_power);
      }
      else
      {
        w = exp(dist2 * iradius16);
      }

      if (need_density)
      {
        w *= 1.0 / (samples_density[t.m_index] * samples_density[t.m_index]);
      }


			if (use_tangent)
			{
				Point3f proj_point = v.P() + v.N() * proj_dist;
				average[i] += proj_point * w;
			}
			else
			{
				average[i] += t.P() * w;
			}

      average_weight_sum[i] += w;  
    }
  }
}

void WLOP::computeAverageTerm(CMesh* samples, CMesh* original)
{
	double average_power = para->getDouble("Average Power");
	bool need_density = para->getBool("Need Compute Density");
	bool use_self_wlop = para->getBool("Run Self WLOP");

//	bool use_elliptical_neighbor = para->getBool("Use Elliptical Original Neighbor");
	bool use_elliptical_neighbor = false;

	bool use_nearest_neighbor = para->getBool("Use Adaptive Sample Neighbor");

	double radius = para->getDouble("CGrid Radius");

	double fix_original_weight = global_paraMgr.skeleton.getDouble("Fix Original Weight");
	//bool use_fixed_original = (para->getBool("Run Skel WLOP") && global_paraMgr.glarea.getBool("Show Samples"));
	//bool use_fixed_original = para->getBool("Run Skel WLOP");

	//bool need_normal_weight = para->getBool("Need Compute Density");
	bool need_normal_weight = para->getBool("Use Elliptical Original Neighbor");
	//bool need_normal_weight = false;
//	bool need_anti_normal = para->getBool("Use Elliptical Original Neighbor");
	bool use_ellipsoid_weight = para->getBool("Use Ellipsoid Weight");


	double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
	double sigma_threshold = pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2);

	double radius2 = radius * radius;
	double iradius16 = -para->getDouble("H Gaussian Para") / radius2;
	double iradius16_proj = -16/ radius2;


	double close_threshold = radius2 / 16;
	bool use_tangent = para->getBool("Use Tangent Vector");


	bool L2 = para->getBool("Need Sample Average");
	bool use_confidence = para->getBool("Use Confidence");

	bool run_anisotropic_lop = para->getBool("Run Anisotropic LOP");
	bool run_skel_lop = para->getBool("Run Skel WLOP");

	if (para->getBool("Run Skel WLOP"))
	{
		use_elliptical_neighbor = use_nearest_neighbor = use_tangent = L2 = use_confidence = need_normal_weight = false;
	}

  bool only_adjust_dlength = para->getBool("Run DLength Adjustment");

	cout << "dual" << endl;
	cout << "Original Size:" << samples->vert[0].original_neighbors.size() << endl;
  if (need_normal_weight)
  {
    cout << "using normal weight!!!!!!!";
  }
  //need_normal_weight = true;


	for (int i = 0; i < samples->vert.size(); i++)
	{

		CVertex& v = samples->vert[i];

		if (i < 5)
		{
			cout << "original size:" << v.original_neighbors.size() << endl;
		}

		for (int j = 0; j < v.original_neighbors.size(); j++)
		{
			CVertex& t = original->vert[v.original_neighbors[j]];

			Point3f diff = t.P() - v.P();
			double proj_dist = diff * v.N();
			//double proj_dist2 = proj_dist * proj_dist;

			
			double dist2 = diff.SquaredNorm();

			double w = 1;

// 			if (use_nearest_neighbor)
// 			{
// 				float new_radius = radius + v.nearest_neighbor_dist;
// 				float new_radius2 = new_radius * new_radius;
// 				iradius16 = -4 / radius2;
// 			}

			if (run_anisotropic_lop)
			{
				double len = sqrt(dist2);
				if (len <= 0.001 * radius) len = radius*0.001;
				double hn = diff * v.N();
				double phi = exp(hn * hn * iradius16);
				w = phi / pow(len, 2 - average_power);
			}
// 			else if (use_ellipsoid_weight && run_skel_lop)
// 			{
// 				if (i < 3 && j < 3)
// 				{
// 					cout << "use_ellipsoid_weight && run_skel_lop" << endl;
// 				}
// 				double ellipsoid_dist0 = (v.eigen_vector0 * diff) / v.eigen_value0;
// 				double ellipsoid_dist1 = (v.eigen_vector1 * diff) / v.eigen_value1;
// 				double ellipsoid_dist2 = (v.eigen_vector2 * diff) / v.eigen_value2;
// 				double ellipsoid_dist_square = ellipsoid_dist0 * ellipsoid_dist0 +
// 					ellipsoid_dist1 * ellipsoid_dist1 + ellipsoid_dist2 * ellipsoid_dist2;
// 
// 				w = exp(dist2 * (-1. / radius2));
// 				//w = exp(ellipsoid_dist_square * iradius16);
// 
// 			}
			else if (average_power < 2)
			{

				double len = sqrt(dist2);
				if (len <= 0.001 * radius) len = radius*0.001;
				w = exp(dist2 * iradius16) / pow(len, 2 - average_power);
			}
			else
			{
				w = exp(dist2 * iradius16);
			}

 			if (need_normal_weight)
 			{

 
 				double normal_diff = exp(-pow(1 - v.N() * t.N(), 2) / sigma_threshold);
 				w *= normal_diff;

        if (i < 15)
        {
          cout << "need_normal_weight " << normal_diff << endl;
        }
 			}

			if (need_density && !original_density.empty())
			{
//  				if (i < 2 && j < 5)
//  				{
// 					cout << "need_original_density " << original_density[t.m_index]  << endl;
//  				}

				w *= original_density[t.m_index];
			}

// 			if (L2 /*|| use_elliptical_neighbor*/ /*|| use_self_wlop*/)
// 			{
// 				if (i < 2)
// 				{
// 					cout << "L2L2L2L2L2" << endl;
// 				}
// 				w = 1.0;
// 			}

// 			if (need_anti_normal)
// 			{
// 				w = 1.0;
// 			}

// 			if (use_fixed_original)
// 			{
// 				if (i < 2)
// 				{
// 					cout << "use_fixed_originaluse_fixed_originaluse_fixed_original" << endl;
// 				}
// 				if (dual_samples->vert[t.m_index].is_fixed_sample)
// 				{
// 					w *= fix_original_weight;
// 				}
// 			}

// 			if (use_tangent)
// 			{
// 				Point3f proj_point = v.P() + v.N() * proj_dist;
// 				average[i] += proj_point * w;
// 			}
// 			else
// 			{
// 				average[i] += t.P() * w;
// 			}
// 			average_weight_sum[i] += w;

// 			if (use_confidence && use_tangent && v.eigen_confidence > 0.0)
// 			{
// 				Point3f proj_point = v.P() + v.N() * proj_dist;
// 				average_low_confidence[i] += proj_point * w;
// 
// 				average[i] += t.P() * w;
// 			}
// 			else

      if (use_tangent /*&& !only_adjust_dlength*/ /*&& !use_confidence*/)
			{
				Point3f proj_point = v.P() + v.N() * proj_dist;
				average[i] += proj_point * w;
			}
			else 
			{
				average[i] += t.P() * w;
			}
			average_weight_sum[i] += w;

		}
	}
}

void WLOP::computeDLinkRepulsionTerm(CMesh* samples)
{
	double radius = para->getDouble("CGrid Radius");
	double radius2 = radius * radius;
	double iradius16 = -para->getDouble("H Gaussian Para") / radius2;
	bool use_tangent = para->getBool("Use Tangent Vector");

	double repulsion_power = para->getDouble("Repulsion Power");
	bool need_density = para->getBool("Need Compute Density");

	double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
	double sigma_threshold = pow(max(1e-8, cos(sigma / 180.0*3.1415926)), 2);

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		CVertex& dual_v = dual_samples->vert[i];

		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& t = samples->vert[v.neighbors[j]];
			
			Point3f diff = v.P() - t.P();

			if (use_tangent)
			{
				diff = GlobalFun::getTangentVector(diff, v.N());
			}

			Point3f d_link = v.P() - dual_v.P();
			//double link_weight = d_link * diff;
			//link_weight = d_link.Normalize() * diff.Normalize();
			//link_weight = 1.0;
			//Point3f repul = diff *(diff * d_link);
			//double direction_diff = exp(pow(d_link.Normalize() * diff.Normalize(), 2) / sigma_threshold);
			//double direction_diff = 1-abs(d_link.Normalize() * diff.Normalize());

			double link_weight = abs(d_link * diff);

			if (i < 5)
			{
				//cout << "dir diff: " << direction_diff << endl;
			}

			double dist2 = diff.SquaredNorm();
			double len = sqrt(dist2);
			if (len <= 0.001 * radius) len = radius*0.001;

			double theta = exp(dist2*iradius16);
			double theta2 = 2 * theta;
			//double weight = theta2 * pow(1.0 / len, repulsion_power);
			double weight = theta2 * link_weight;

			repulsion[i] += diff * weight;
			repulsion_weight_sum[i] += weight;
		}
	}
}


void WLOP::computeRepulsionTerm(CMesh* samples)
{
	double repulsion_power = para->getDouble("Repulsion Power");
	bool need_density = para->getBool("Need Compute Density");
	double radius = para->getDouble("CGrid Radius"); 

	double radius2 = radius * radius;
	double iradius16 = -para->getDouble("H Gaussian Para")/radius2;
  bool use_tangent = para->getBool("Use Tangent Vector");
	bool run_skel_wlop = para->getBool("Run Skel WLOP");
	bool run_tangential_movement = para->getBool("Run Tangential Motion");


	if (use_ellipsoid_repulsion && (run_skel_wlop || para->getBool("Run Move Sample")))
	{
		repulsion_power = para->getDouble("Big Repulsion Power");
	}

	if (run_tangential_movement)
	{
		repulsion_power = para->getDouble("Big Repulsion Power");
	}


	if (para->getBool("Run Skel WLOP"))
	{
		use_tangent = false;
	}
	
	bool run_skel_lop = para->getBool("Run Skel WLOP");
	cout << endl<< endl<< "Sample Neighbor Size:" << samples->vert[0].neighbors.size() << endl<< endl;
	for(int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];

		Point3f eigenX = v.eigen_vector0 * v.eigen_value0;
		Point3f eigenY = v.eigen_vector1 * v.eigen_value1;
		Point3f eigenZ = v.eigen_vector2 * v.eigen_value2;

		//Eigen::Matrix3d eigen_MAT;
		Matrix33f eigen_MAT;
		eigen_MAT.SetRow(0, eigenX);
		eigen_MAT.SetRow(1, eigenY);
		eigen_MAT.SetRow(2, eigenZ);

		Matrix33f eigen_Trans;
		eigen_Trans = eigen_MAT.Transpose();

		if (i < 4 && use_ellipsoid_weight && run_skel_lop)
		{
			cout << "eigen matric:!! " << eigen_Trans[0][0] <<
				"  " << eigen_Trans[0][1] << "  " << eigen_Trans[1][1] << "  " << eigen_Trans[2][2] << endl;
		}

		Matrix33f eigen_multi = eigen_Trans * eigen_MAT;
		Point3f col0 = eigen_multi.GetColumn(0);
		Point3f col1 = eigen_multi.GetColumn(1);
		Point3f col2 = eigen_multi.GetColumn(2);

		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& t = samples->vert[v.neighbors[j]];
			Point3f diff = v.P() - t.P();

			//Point3f diff2 = v.P() - t.P();;
      if (use_tangent || run_tangential_movement)
      {
        diff = GlobalFun::getTangentVector(diff, v.N());
      }

			double dist2 = diff.SquaredNorm();
			double len = sqrt(dist2);
			if(len <= 0.001 * radius) len = radius*0.001;

			double w = exp(dist2*iradius16);

			if (use_ellipsoid_weight && run_skel_lop)
			{
				if (i < 3)
				{
					cout << "use_ellipsoid_weight && run_skel_lop" << endl;
				}
				double ellipsoid_dist0 = (v.eigen_vector0 * diff) / v.eigen_value0;
				double ellipsoid_dist1 = (v.eigen_vector1 * diff) / v.eigen_value1;
				double ellipsoid_dist2 = (v.eigen_vector2 * diff) / v.eigen_value2;
				double ellipsoid_dist_square = ellipsoid_dist0 * ellipsoid_dist0 +
					                             ellipsoid_dist1 * ellipsoid_dist1 + 
																			 ellipsoid_dist2 * ellipsoid_dist2;

				//w = exp(dist2 * (-1. / radius2));
				w = exp(ellipsoid_dist_square * iradius16);
			}

			double rep = w * pow(1.0 / len, repulsion_power);
			//double rep = w;

			if (need_density)//2015
			{
// 				if (i<2 && j<3)
// 				{
// 					cout << "sample density: " << samples_density[t.m_index] << endl;
// 				}
				rep *= samples_density[t.m_index];
			}

			if (use_ellipsoid_repulsion && run_skel_wlop)
			{
				repulsion_weight_sum[i] += rep;

				//cout << "i " << i << endl;
				repulsion_x_length[i] += diff * v.eigen_vector0 * rep;
				repulsion_y_length[i] += diff * v.eigen_vector1 * rep;
				repulsion_z_length[i] += diff * v.eigen_vector2 * rep;
			}
			else
			{
				repulsion[i] += diff * rep;
				repulsion_weight_sum[i] += rep;
			}
		}
	}
}

void WLOP::computeSampleAverageTerm(CMesh* samples)
{
  double repulsion_power = para->getDouble("Repulsion Power");
  bool need_density = para->getBool("Need Compute Density");
  double radius = para->getDouble("CGrid Radius"); 

  double radius2 = radius * radius;
  double iradius16 = -para->getDouble("H Gaussian Para")/radius2;

  cout << endl<< endl<< "Sample Neighbor Size:" << samples->vert[0].neighbors.size() << endl<< endl;
  for(int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];
    for (int j = 0; j < v.neighbors.size(); j++)
    {
      CVertex& t = samples->vert[v.neighbors[j]];
      Point3f diff = v.P() - t.P();

      double dist2  = diff.SquaredNorm();

      double w = exp(dist2*iradius16);

      if (need_density)
      {
        w *= samples_density[t.m_index];
      }

      samples_average[i] += t.P() * w;  
      samples_average_weight_sum[i] += w;
    }
  }
}



// fail??
// void WLOP::computeSampleSimilarityTerm(CMesh* samples)
// {
// 	GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
// 	GlobalFun::computeEigenWithTheta(dual_samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));
// 
// 	GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 1, false, "WlopParaDlg::computeSampleSimilarityTerm()");
// 
// 	double radius = para->getDouble("CGrid Radius") *2.0;
// 	double radius2 = radius * radius;
// 	double iradius16 = -4.0 / radius2;
// 	double iradius16_perpend = -4.0 / radius2;
// 
// 	double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
// 	double sigma_threshold = pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2);
// 
// 	bool use_confidence = para->getBool("Use Confidence");
// 	if (use_confidence)
// 	{
// 		cout << "use confidence" << endl;
// 	}
// 
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		CVertex& v = samples->vert[i];
// 
// 		int neighbor_idx = v.neighbors[0];
// 
// 		CVertex& dual_v = dual_samples->vert[neighbor_idx];
// 
// 		Point3f diff = v.P() - dual_v.P();
// 		
// 		double proj_dist = diff * v.N();		
// 		//v.skel_radius = abs(proj_dist);
// 
// 		Point3f proj_p = dual_v.P() + dual_v.eigen_vector0 * proj_dist;
//     v.skel_radius = GlobalFun::computeEulerDist(v.P(), proj_p);
// 		
// 		v.dual_index = neighbor_idx;
// 
// 		samples->bbox.Add(v.P());
// 	}
// 
// 	//GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, 50, false, "WlopParaDlg::computeSampleSimilarityTerm()");
// 	GlobalFun::computeBallNeighbors(samples, NULL, radius, samples->bbox);
// 
// 
// 	vector<double> new_radiuses;
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		new_radiuses.push_back(samples->vert[i].skel_radius);
// 	}
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		CVertex& v = samples->vert[i];
// 		CVertex& dual_v = dual_samples->vert[v.dual_index];
// 		Point3f v_outward_direction = (v.P() - dual_v.P()).Normalize();
// 
// 		if (use_kite_points && !v.is_boundary)
// 		{
// 			continue;
// 		}
// 
// 		double sum_radius = 0;
// 		double sum_weight = 0;
// 		double weight = 1;
// 
// 		if (v.neighbors.size() < 3)
// 		{
// 			continue;
// 		}
// 
// 		for (int j = 0; j < v.neighbors.size(); j++)
// 		{
// 			CVertex& t = samples->vert[v.neighbors[j]];
// 			CVertex& dual_t = dual_samples->vert[t.dual_index];
// 
// 			Point3f t_outward_direction = (t.P() - dual_t.P()).Normalize();
// 
// 			float dist2 = (v.P() - t.P()).SquaredNorm();
// 
// 			float dist_diff = exp(dist2 * iradius16);
// 			//double direction_diff = exp(-pow(1 - pow(v_outward_direction * t_outward_direction, 2), 2) / sigma_threshold);
// 
// 			//double direction_diff = exp(-pow(1 - v_outward_direction * t_outward_direction, 2) / sigma_threshold);
// 			double direction_diff = exp(-pow(1 - v.N() * t.N(), 2) / sigma_threshold);
// 			weight = direction_diff * dist_diff;
// 
// 
// 			if (use_confidence)
// 			{
// 				weight *= (t.eigen_confidence * t.eigen_confidence);
// 			}
// 
// 			sum_radius += t.skel_radius * weight;
// 			sum_weight += weight;
// 		}
// 
// 		double new_radius = v.skel_radius;
// 
// 		if (!v.neighbors.empty())
// 		{
// 			new_radiuses[i] = sum_radius / sum_weight;
// 		}
// 	}
// 
// 	GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 1, false, "WlopParaDlg::computeSampleSimilarityTerm()");
// 
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		CVertex& v = samples->vert[i];
// 		v.recompute_m_render();
// 
// 		int neighbor_idx = v.neighbors[0];
// 		CVertex dual_v = dual_samples->vert[neighbor_idx];
// 		double movement = v.skel_radius - new_radiuses[i];
// 
// 		Point3f diff = v.P() - dual_v.P();
// 
// // 		if (diff * v.N() < 0)
// // 		{
// // 			samples_similarity[i] = v.P() - v.N() * movement;
// // 		}
// // 		else
// // 		{
// // 			samples_similarity[i] = v.P() + v.N() * movement;
// // 		}
// 
//  		double proj_dist = diff * dual_v.eigen_vector0;
//  		Point3f proj_p = dual_v.P() + dual_v.eigen_vector0 * proj_dist;
//  
//  		//Point3f direction = (v.P() - dual_v.P()).Normalize();
//  		Point3f direction = (v.P() - proj_p).Normalize();
//  
//  		//samples_similarity[i] = dual_v.P() + direction * new_radiuses[i];
//  		samples_similarity[i] = proj_p + direction * new_radiuses[i];
//  
//  		//samples_similarity[i] = v.P();
// 	}
// }


// 21-11-2014
//void WLOP::computeSampleSimilarityTerm(CMesh* samples)
//{
//	GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
//	GlobalFun::computeEigenWithTheta(dual_samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));
//
//	GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 1, false, "WlopParaDlg::computeSampleSimilarityTerm()");
//
//	double radius = para->getDouble("CGrid Radius") * 2.0;
//	double radius2 = radius * radius;
//	double iradius16 = -4.0 / radius2;
//	double iradius16_perpend = -4.0 / radius2;
//
//	double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
//	double sigma_threshold = pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2);
//
//	bool use_confidence = para->getBool("Use Confidence");
//	if (use_confidence)
//	{
//		cout << "use confidence" << endl;
//	}
//
//	for (int i = 0; i < samples->vert.size(); i++)
//	{
//		CVertex& v = samples->vert[i];
//
//		int neighbor_idx = v.neighbors[0];
//
//		CVertex& dual_v = dual_samples->vert[neighbor_idx];
//
//		Point3f diff = v.P() - dual_v.P();
//		double proj_dist = diff * dual_v.eigen_vector0;
//		Point3f proj_p = dual_v.P() + dual_v.eigen_vector0 * proj_dist;
//
//		v.skel_radius = GlobalFun::computeEulerDist(v.P(), proj_p);
//		v.dual_index = neighbor_idx;
//
//		samples->bbox.Add(v.P());
//	}
//
//	//GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, 50, false, "WlopParaDlg::computeSampleSimilarityTerm()");
//	GlobalFun::computeBallNeighbors(samples, NULL, radius, samples->bbox);
//
//
//	vector<double> new_radiuses;
//	for (int i = 0; i < samples->vert.size(); i++)
//	{
//		new_radiuses.push_back(samples->vert[i].skel_radius);
//	}
//	for (int i = 0; i < samples->vert.size(); i++)
//	{
//		CVertex& v = samples->vert[i];
//		CVertex& dual_v = dual_samples->vert[v.dual_index];
//		Point3f v_outward_direction = (v.P() - dual_v.P()).Normalize();
//
//		if (use_kite_points && !v.is_boundary)
//		{
//			continue;
//		}
//
//		double sum_radius = 0;
//		double sum_weight = 0;
//		double weight = 1;
//
//		if (v.neighbors.size() < 3)
//		{
//			continue;
//		}
//
//		for (int j = 0; j < v.neighbors.size(); j++)
//		{
//			CVertex& t = samples->vert[v.neighbors[j]];
//			CVertex& dual_t = dual_samples->vert[t.dual_index];
//
//			Point3f t_outward_direction = (t.P() - dual_t.P()).Normalize();
//
//			float dist2 = (v.P() - t.P()).SquaredNorm();
//
//			float dist_diff = exp(dist2 * iradius16);
//			//double direction_diff = exp(-pow(1 - pow(v_outward_direction * t_outward_direction, 2), 2) / sigma_threshold);
//			
//			double direction_diff = exp(-pow(1 - v_outward_direction * t_outward_direction, 2) / sigma_threshold);
//			weight = direction_diff * dist_diff;
//
//
//			if (use_confidence)
//			{
//				weight *= (t.eigen_confidence * t.eigen_confidence);
//			}
//
//			sum_radius += t.skel_radius * weight;
//			sum_weight += weight;
//		}
//
//		double new_radius = v.skel_radius;
//
//		if (!v.neighbors.empty())
//		{
//			new_radiuses[i] = sum_radius / sum_weight;
//		}
//	}
//
//	GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 1, false, "WlopParaDlg::computeSampleSimilarityTerm()");
//
//	for (int i = 0; i < samples->vert.size(); i++)
//	{
//		CVertex& v = samples->vert[i];
//		v.recompute_m_render();
//
//		int neighbor_idx = v.neighbors[0];
//		CVertex dual_v = dual_samples->vert[neighbor_idx];
//
//		Point3f diff = v.P() - dual_v.P();
//		double proj_dist = diff * dual_v.eigen_vector0;
//		Point3f proj_p = dual_v.P() + dual_v.eigen_vector0 * proj_dist;
//
//		//Point3f direction = (v.P() - dual_v.P()).Normalize();
//		Point3f direction = (v.P() - proj_p).Normalize();
//
//		//samples_similarity[i] = dual_v.P() + direction * new_radiuses[i];
//		samples_similarity[i] = proj_p + direction * new_radiuses[i];
//
//		//samples_similarity[i] = v.P();
//	}
//}


void WLOP::computeSampleSimilarityTerm(CMesh* samples)
{
	double simi_neighbor_radius_para = para->getDouble("Similarity Term Neighbor Para");
	double length_threshold = para->getDouble("Similarity Length Outlier Threshold");

	computeDualIndex(samples, dual_samples);

	double radius = para->getDouble("CGrid Radius") * simi_neighbor_radius_para;
	double radius2 = radius * radius;
	double iradius16 = -4.0 / radius2;
	double iradius16_perpend = -4.0 / radius2;

	double length_threshold_dist = para->getDouble("CGrid Radius") /** length_threshold*/;
	double length_threshold_dist2 = length_threshold_dist * length_threshold_dist;
	double iradius16_length = -4.0 / length_threshold_dist2;

	double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
	double sigma_threshold = pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2);

	int similarity_KNN = para->getDouble("KNN For Similarity");

	bool use_confidence = para->getBool("Use Confidence");
	if (use_confidence)
	{
		cout << "use confidence in similarity term" << endl;
	}

	// compute current lengthes
	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		CVertex& dual_v = dual_samples->vert[v.dual_index];

		Point3f diff = v.P() - dual_v.P();
		//double proj_dist = abs(diff * v.N());
		double proj_dist = sqrt(diff.SquaredNorm());

		v.skel_radius = proj_dist;
		samples->bbox.Add(v.P());
	}

	// compute neighbor from dual points
	if (similarity_KNN > 2)
	{
		GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, similarity_KNN, false, "Similarity KNN");
		GlobalFun::computeAnnNeigbhors(dual_samples->vert, dual_samples->vert, similarity_KNN, false, "Similarity KNN");
	}
	else
	{
		GlobalFun::computeBallNeighbors(samples, NULL, radius, samples->bbox);
		GlobalFun::computeBallNeighbors(dual_samples, NULL, radius, samples->bbox);
	}

	double length_threshold_percentage = para->getDouble("Similarity Length Outlier Threshold");


	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		CVertex& dual_v = dual_samples->vert[i];

		set<int> neighbors_set;
		for (int j = 0; j < v.neighbors.size(); j++)
		{
			neighbors_set.insert(v.neighbors[j]);
		}
		for (int j = 0; j < dual_v.neighbors.size(); j++)
		{
			neighbors_set.insert(dual_v.neighbors[j]);
		}
		set<int>::iterator iter;
		v.neighbors.clear();
		for (iter = neighbors_set.begin(); iter != neighbors_set.end(); ++iter)
		{
			v.neighbors.push_back(*iter);
		}
		//samples->vert[i].neighbors = dual_samples->vert[i].neighbors;
	}

	vector<double> new_radiuses;
	for (int i = 0; i < samples->vert.size(); i++)
	{
		new_radiuses.push_back(samples->vert[i].skel_radius);
	}

	// compute averaging lengthes
	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		CVertex& dual_v = dual_samples->vert[v.dual_index];
		Point3f v_outward_direction = v.N();

// 		if (use_kite_points && !v.is_boundary)
// 		{
// 			continue;
// 		}

		double sum_radius = 0;
		double sum_weight = 0;
		double weight = 1;

		if (v.neighbors.size() < 2)
		{
			cout << "similarity term too small neighborhood" << endl;
			continue;
		}

		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& t = samples->vert[v.neighbors[j]];
			CVertex& dual_t = dual_samples->vert[t.dual_index];

			Point3f t_outward_direction = t.N();

			//double dist2 = (v.P() - t.P()).SquaredNorm();
      double dist2 = (dual_v.P() - dual_t.P()).SquaredNorm();

			double dist_diff = exp(dist2 * iradius16);

			//double sign_dist = t.skel_radius - v.skel_radius;
			//double length_dist = abs(sign_dist);
//  			if (sign_dist > length_threshold_dist)
//  			{
//  				continue;
//  			}
// 
//   			double length_dist = v.skel_radius - t.skel_radius;
//   			double length_dist2 = length_dist *length_dist;
// 				double length_diff = exp(length_dist2 * iradius16_length);
			//double direction_diff = exp(-pow(1 - pow(v_outward_direction * t_outward_direction, 2), 2) / sigma_threshold);

			//double direction_diff = exp(-pow(1 - v_outward_direction * t_outward_direction, 2) / sigma_threshold);
      double angle_diff = abs(v_outward_direction * t_outward_direction);
      double direction_diff = exp(-pow(1 - v_outward_direction * t_outward_direction, 2) / sigma_threshold);


// 			if (similarity_KNN > 2)
// 			{
// 				direction_diff = 1.0;
// 			}
			
			weight = direction_diff * dist_diff /** length_diff*/;


			if (use_confidence && t.eigen_confidence > 1e-5)
			{
				weight *= (t.eigen_confidence * t.eigen_confidence);
			}

			sum_radius += t.skel_radius * weight;
			sum_weight += weight;
		}

		if (!v.neighbors.empty() && sum_weight > 1e-5)
		{
			new_radiuses[i] = sum_radius / sum_weight;
		}
	}

	computeDualIndex(samples, dual_samples);
	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		v.recompute_m_render();

		CVertex dual_v = dual_samples->vert[v.dual_index];

// 		Point3f backward_v = v.P() - v.N() * v.skel_radius;
// 		Point3f forward_v = backward_v + v.N() * new_radiuses[i];

		Point3f direction = (v.P() - dual_v.P()).Normalize();
		Point3f forward_v = dual_v.P() + direction * new_radiuses[i];
		samples_similarity[i] = forward_v;

// 		Point3f diff = v.P() - dual_v.P();
// 		double proj_dist = diff * dual_v.eigen_vector0;
// 		Point3f proj_p = dual_v.P() + dual_v.eigen_vector0 * proj_dist;
// 
// 		//Point3f direction = (v.P() - dual_v.P()).Normalize();
// 		Point3f direction = (v.P() - proj_p).Normalize();
// 
// 		//samples_similarity[i] = dual_v.P() + direction * new_radiuses[i];
// 		samples_similarity[i] = proj_p + direction * new_radiuses[i];
		
		//samples_similarity[i] = v.P();

		//sample
	}
}


void WLOP::computeDensity(bool isOriginal, double radius, CMesh* samples, CMesh* original)
{
	CMesh* mesh;
	if (isOriginal)
	{
		mesh = original;
	}
	else
	{
// 		if (para->getBool("Run Skel WLOP"))
// 		{
// 			cout << "density for Skel WLOP Sample" << endl;
// 			mesh = dual_samples;
// 		}
		mesh = samples;
	}

	double radius2 = radius * radius;
	double iradius16 = -para->getDouble("H Gaussian Para") / radius2;

	for(int i = 0; i < mesh->vert.size(); i++)
	{
		CVertex& v = mesh->vert[i];

		if (isOriginal)
		{
			original_density[i] = 1.;
		}
		else
		{
			samples_density[i] = 1.;
		}

		vector<int>* neighbors = &v.neighbors;

		for (int j = 0; j < neighbors->size(); j++)
		{
			CVertex& t = mesh->vert[(*neighbors)[j]];
			double dist2  = (v.P() - t.P()).SquaredNorm();
			double den = exp(dist2*iradius16);

			if (isOriginal)
			{
				original_density[i] += den;
// 
// 				if (i<2 && j <2)
// 				{
// 					cout << "density::::" << den << endl;
// 				}
			}
			else
			{
				samples_density[i] += den;
			}
		}
	}

	for(int i = 0; i < mesh->vert.size(); i++)
	{
		if (isOriginal)
		{
			CVertex& v = mesh->vert[i];
			original_density[i] = 1. / original_density[i];
		}
		else
		{
			samples_density[i] = sqrt(samples_density[i]);
		}
	}

}


void WLOP::runDlengthAdjustment()
{
  use_adaptive_mu = para->getBool("Use Adaptive Mu");
  is_sample_close_to_original.assign(samples->vert.size(), false);
  bool use_tangent = para->getBool("Use Tangent Vector");

  Timer time;

  initVertexes(true);

// 
//   time.start("Sample Sample Neighbor Tree");
//   GlobalFun::computeBallNeighbors(samples, NULL,
//     para->getDouble("CGrid Radius"), samples->bbox);
//   time.end();


//   if (nTimeIterated == 0)
//   {
//     if (para->getBool("Need Compute Density"))
//     {
//       double local_density_para = 1.0;
//       time.start("Original Original Neighbor Tree");
//       GlobalFun::computeBallNeighbors(original, NULL,
//         para->getDouble("CGrid Radius") * local_density_para, original->bbox);
//       time.end();
// 
//       time.start("Compute Original Density");
//       original_density.assign(original->vn, 0);
// 
//       computeDensity(true, para->getDouble("CGrid Radius") * local_density_para, samples, original);
//       time.end();
//     }
//   }

//   if (para->getBool("Need Compute Density"))
//   {
//     time.start("Compute Density For Sample");
//     computeDensity(false, para->getDouble("CGrid Radius"), samples, original);
//     time.end();
//   }


  time.start("Sample Original Neighbor Tree!!!");
  GlobalFun::computeBallNeighbors(samples, original,
    para->getDouble("CGrid Radius"), box);
  time.end();
  
  time.start("Compute Average Term");
  computeAverageTerm(samples, original);
  time.end();

  for (int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];

    if (global_paraMgr.upsampling.getBool("Use Adaptive Upsampling") && !v.is_fixed_sample)
    {
      continue;
    }

    Point3f avg_point = v.P();
    if (average_weight_sum[i] > 1e-20)
    {
      avg_point = average[i] / average_weight_sum[i];
    }

    v.P() = avg_point;
  }
}

double WLOP::iterate()
{
  use_adaptive_mu = para->getBool("Use Adaptive Mu");
  is_sample_close_to_original.assign(samples->vert.size(), false);
  bool use_tangent = para->getBool("Use Tangent Vector");

	if (para->getBool("Use Adaptive Sample Neighbor"))
	{
		computeNearestNeighborDist();
	}

	Timer time;

	initVertexes(true);


	time.start("Sample Sample Neighbor Tree");
	GlobalFun::computeBallNeighbors(samples, NULL,
		para->getDouble("CGrid Radius"), samples->bbox);
	time.end();


	if (nTimeIterated == 0) 
	{
		if (para->getBool("Need Compute Density"))
		{
			double local_density_para = 1.0;
			time.start("Original Original Neighbor Tree");
			GlobalFun::computeBallNeighbors(original, NULL, 
				para->getDouble("CGrid Radius") * local_density_para, original->bbox);
			time.end();

			time.start("Compute Original Density");
			original_density.assign(original->vn, 0);

			computeDensity(true, para->getDouble("CGrid Radius") * local_density_para, samples, original);
			time.end();
		}
	}

	if (para->getBool("Need Compute Density"))
	{
		time.start("Compute Density For Sample");
		computeDensity(false, para->getDouble("CGrid Radius"), samples, original);
		time.end();
	}

	if (para->getBool("Use Original Averaging KNN"))
	{
		int knn = para->getDouble("Original Averaging KNN");
		time.start("Sample Original Neighbor Tree!!!");
		GlobalFun::computeAnnNeigbhors(original->vert, samples->vert, knn, false, "averaging KNN");
		for (int i = 0; i < samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];
			v.original_neighbors = v.neighbors;
		}
		time.end();

		time.start("Sample Sample Neighbor Tree");
		GlobalFun::computeBallNeighbors(samples, NULL,
			para->getDouble("CGrid Radius"), samples->bbox);
		time.end();
	}
	else
	{
		time.start("Sample Original Neighbor Tree!!!");
		GlobalFun::computeBallNeighbors(samples, original,
			para->getDouble("CGrid Radius"), box);
		time.end();
	}

   time.start("Compute Average Term");
 	computeAverageTerm(samples, original);
 	time.end();


 	time.start("Compute Repulsion Term");
 	computeRepulsionTerm(samples);
 	time.end();


  if (para->getBool("Need Similarity"))
  {
    time.start("Compute Similarity Term");
    computeSampleSimilarityTerm(samples);
    time.end();
 	}

	
  vector<Point3f> new_sample_positions;
  vector<float> move_proj_vec;

  double radius = para->getDouble("CGrid Radius");
  double radius2 = radius * radius;
  double iradius16 = -para->getDouble("H Gaussian Para")/radius2;

	bool use_confidence = para->getBool("Use Confidence");
	if (use_confidence)
	{
		GlobalFun::normalizeConfidence(samples->vert, 0);
	}

	int error_x = 0;
  if (para->getBool("Need Averaging Movement"))
  {
		vector<Point3f> new_sample_positions = computeNewSamplePositions(error_x);
		vector<float> move_proj_vec;
		vector<Point3f> tangent_moves;
		move_proj_vec.assign(samples->vert.size(), 0.);
		tangent_moves.assign(samples->vert.size(), Point3f(0., 0., 0.));

		for (int i = 0; i < samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];
			Point3f temp_p = new_sample_positions[i];

			Point3f move_vector = temp_p - v.P();
			float move_proj = move_vector * v.N();

			move_proj_vec[i] = move_proj;
			tangent_moves[i] = move_vector - v.N() * move_proj;
		}
 
		for (int i = 0; i < samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];

			float sum_move_proj = 0.0;
			float sum_w = 0.0;

			if (v.neighbors.empty())
			{
				continue;
			}

			for (int j = 0; j < v.neighbors.size(); j++)
			{
				int neighbor_idx = v.neighbors[j];
				CVertex& t = samples->vert[neighbor_idx];
				float neighbor_move_proj = move_proj_vec[neighbor_idx];

				float dist2 = GlobalFun::computeEulerDistSquare(v.P(), t.P());
				float w = exp(dist2 * iradius16);

				sum_move_proj += w * neighbor_move_proj;
				sum_w += w;
			}

			float avg_move_proj = sum_move_proj / sum_w;

			if (sum_w > 0.0)
			{
				//v.P() += (v.N() * avg_move_proj) +tangent_moves[i];
				v.P() = new_sample_positions[i];

			}

		}
  }
  else
  {
		vector<Point3f> new_positions = computeNewSamplePositions(error_x);

		for (int i = 0; i < new_positions.size(); i++)
		{
			samples->vert[i].P() = new_positions[i];

			//samples->vert[i].N() = repulsion[i] * (0.5 / repulsion_weight_sum[i]);;
		}
  }

	para->setValue("Current Movement Error", DoubleValue(error_x));
	cout << "****finished compute WLOP error:	" << error_x << endl;

	if (para->getBool("Need Compute PCA"))
	{
		time.start("Recompute PCA");
		recomputePCA_Normal();
		time.end();
	}

	return error_x;
}


vector<Point3f> WLOP::computeNewSamplePositions(int& error_x)
{
	double mu = para->getDouble("Repulsion Mu");
	bool need_similarity = para->getBool("Need Similarity");
	bool use_tangent = para->getBool("Use Tangent Vector");
	bool use_confidence = para->getBool("Use Confidence");
	bool only_do_repulsion = para->getBool("Only Do Repuslion");
	bool only_do_average = para->getBool("Only Do Repuslion");

	bool use_self_wlop = para->getBool("Run Self WLOP");

	double average_dist = para->getDouble("Average Closest Dist");

	double radius = para->getDouble("CGrid Radius");

	double save_move_threshold_along_normal = para->getDouble("CGrid Radius") * para->getDouble("Save Move Dist Along Normal Para");

	double radius_threshold = radius * 0.6;

	bool use_confidence_to_combine = para->getBool("Use Confidence To Merge");

	double protect_small_tubular_radius_para = para->getDouble("Protect Small Tubular Para");
	double protect_high_confidence_para = para->getDouble("Protect High Confidence Para");
	double data_outweigh_simimlarity_para = para->getDouble("Data Outweigh Similarity Para");

  bool only_adjust_dlength = para->getBool("Run DLength Adjustment");

	Point3f c;

	vector<Point3f> new_pos(samples->vert.size());
	for (int i = 0; i < new_pos.size(); i++)
	{
		new_pos[i] = samples->vert[i].P();
	}

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		CVertex& dual_v = dual_samples->vert[i];
		Point3f diff = v.P() - dual_v.P();
		double dlink_length = sqrt(diff.SquaredNorm());

		c = v.P();

		//save_move_threshold_along_normal = dlink_length * para->getDouble("Save Move Dist Along Normal Para");

// 		if (v.is_fixed_sample)
// 		{
// 			new_pos[i] = v.P();
// 			continue;
// 		}

		if (use_tangent)
		{

			if (use_self_wlop)
			{
				if (v.is_fixed_sample || !(repulsion_weight_sum[i] > 1e-20 && mu > 0))
				{
					new_pos[i] = v.P();
				}
				else
				{
					new_pos[i] = v.P() + repulsion[i] * (mu / repulsion_weight_sum[i]);
				}
			}
			else if (only_do_repulsion && repulsion_weight_sum[i] > 1e-20 && mu > 0)
			{
				new_pos[i] = v.P() + repulsion[i] * (mu / repulsion_weight_sum[i]);
			}
			else if (only_do_average && average_weight_sum[i] > 1e-20)
			{
				new_pos[i] = average[i] / average_weight_sum[i];
			}
			else if (need_similarity)
			{
        if (use_confidence /*&& v.eigen_confidence < 0.85*/)
				{
					//using this
          Point3f avg_point = v.P();
          if (average_weight_sum[i] > 1e-20)
          {
            avg_point = average[i] / average_weight_sum[i];
          }
					Point3f sim_point = samples_similarity[i];

					double dist_avg = GlobalFun::computeEulerDist(v.P(), avg_point);
					if (dist_avg > save_move_threshold_along_normal)
					{
						avg_point = v.P();
						//v.is_skel_virtual = true;
						//continue;
					}
					else
					{
						//v.is_skel_virtual = false;
					}

					double dist_sim = GlobalFun::computeEulerDist(v.P(), sim_point);
					if (dist_sim > save_move_threshold_along_normal /*&& (sim_point - v.P())*v.N() > 0*/)
					{
						sim_point = v.P();
						//v.is_skel_branch = true;
						//continue;
					}
					else
					{
						//v.is_skel_branch = false;
					}


					double dist = GlobalFun::computeEulerDist(avg_point, sim_point);


//  					else if (v.eigen_confidence < 0.4)
//  					{
//  						new_pos[i] = sim_point;
//  					}
					v.is_skel_virtual = false;
					v.is_boundary = false;
					v.is_skel_branch = false;
					v.is_fixed_original = false;

					if (v.eigen_confidence > protect_high_confidence_para)
					{
						new_pos[i] = avg_point;
						v.is_skel_branch = true; // blue

					}
					else if (dual_v.eigen_confidence > 0.95 && dlink_length < (average_dist*protect_small_tubular_radius_para))
					{

						new_pos[i] = sim_point;
						v.is_skel_virtual = true; // gray

					}
					else if (dist < radius_threshold)
					{
            
						if (use_confidence_to_combine && v.eigen_confidence > 0)
						{
              cout << "rare points" << endl;

							//new_pos[i] = avg_point * 0.5 + sim_point * 0.5;
							if (v.eigen_confidence > 0.99)
							{
								new_pos[i] = avg_point;
							}
							else if (v.eigen_confidence < 0.1)
							{
								new_pos[i] = sim_point;
							}
							else
							{
								new_pos[i] = avg_point * v.eigen_confidence + sim_point * (1 - v.eigen_confidence);

								v.is_boundary = true; //orange
							}
						}
						else
						{
							if (dist_avg < dist_sim * data_outweigh_simimlarity_para)
							{
								new_pos[i] = avg_point;

								v.is_boundary = true; //orange
								v.is_fixed_original = true; // green // black

							}
							else
							{
								new_pos[i] = sim_point;
								v.is_fixed_original = true; // green
							}
							//new_pos[i] = avg_point * 0.5 + sim_point * 0.5;
						}
					}
					else
					{
						if (dist_avg < dist_sim * data_outweigh_simimlarity_para)
						{
							new_pos[i] = avg_point;

							v.is_boundary = true; //orange
							v.is_fixed_original = true; // green // black

						}
						else
						{
							new_pos[i] = sim_point;

							v.is_fixed_original = true; // green 

						}
					}
				}
				else if (average_weight_sum[i] > 1e-20)
				{
					if (dual_v.eigen_confidence > 0.95 && dlink_length < (average_dist*3.0))
					{
						Point3f sim_point = samples_similarity[i];
						double dist_sim = GlobalFun::computeEulerDist(v.P(), sim_point);
						if (dist_sim > save_move_threshold_along_normal && (sim_point - v.P())*v.N() > 0)
						{
							sim_point = v.P();
						}

						new_pos[i] = sim_point;
						v.is_skel_branch = true;
						//cout << "tubular" << endl;
					}
					else
					{
						new_pos[i] = average[i] / average_weight_sum[i];
					}
				}
			}
			else if (average_weight_sum[i] > 1e-20)
			{
				new_pos[i] = average[i] / average_weight_sum[i];
			}

      if (!only_adjust_dlength)
      {
        if (repulsion_weight_sum[i] > 1e-20 && mu > 0)
        {
          new_pos[i] += repulsion[i] * (mu / repulsion_weight_sum[i]);
        }
        else
        {
          cout << "maybe no neighbor for repulsion" << endl;
          v.is_skel_virtual = true;
        }

        if (repulsion_weight_sum[i] > 1e-20)
        {
          Point3f diff = v.P() - c;
          double move_error = sqrt(diff.SquaredNorm());
          error_x += move_error;
        }
      }

			//}
		}
		else
		{

			if (mu > 1.0)
			{
				new_pos[i] = v.P();
			}
			else if (average_weight_sum[i] > 1e-20)
			{
				new_pos[i] = average[i] / average_weight_sum[i];
			}

			if (repulsion_weight_sum[i] > 1e-20 && mu > 0)
			{
				new_pos[i] += repulsion[i] * (mu / repulsion_weight_sum[i]);
			}



			if (repulsion_weight_sum[i] > 1e-20)
			{
				Point3f diff = v.P() - c;
				double move_error = sqrt(diff.SquaredNorm());
				error_x += move_error;
			}
		}
	}
	error_x = error_x / samples->vn;

	return new_pos;
}



void WLOP::recomputePCA_Normal()
{
  CMesh temp_mesh;
  //if (para->getBool("Use Adaptive Sample Neighbor"))
  if (false)
  {
    for(int i = 0; i < samples->vn; i++)
    { 
      temp_mesh.vert.push_back(samples->vert[i]);
    }
    temp_mesh.vn = samples->vert.size();

    GlobalFun::computeUndirectedNormal(&temp_mesh, para->getDouble("CGrid Radius"));
 
    for(int i = 0; i < samples->vn; i++)
    {
      Point3f& new_normal = temp_mesh.vert[i].N();
      //Point3f& new_normal = temp_mesh.vert[i].eigen_vector1;

      CVertex& v = samples->vert[i];
      if (v.N() * new_normal > 0)
      {
        v.N() = new_normal;
      }
      else
      {
        v.N() = -new_normal;
      }
      v.recompute_m_render();
    }
  }
  else
  {
		
    for(int i = 0; i < samples->vn; i++)
    {
      mesh_temp[i].P() = samples->vert[i].P();
    }

    int knn = global_paraMgr.norSmooth.getInt("PCA KNN");
    vcg::NormalExtrapolation<vector<CVertex> >::ExtrapolateNormals(mesh_temp.begin(), mesh_temp.end(), knn, -1);

    for(int i = 0; i < samples->vn; i++)
    {
      Point3f& new_normal = mesh_temp[i].N();
      CVertex& v = samples->vert[i];
      if (v.N() * new_normal > 0)
      {
        v.N() = new_normal;
      }
      else
      {
        v.N() = -new_normal;
      }
      v.recompute_m_render();
    }
  }

}

void WLOP::stepForward()
{
  float pace_size = 0.1 * para->getDouble("CGrid Radius");
  for (int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];
    v.P() += v.N() * pace_size;
  }
  //recomputePCA_Normal();
}


void WLOP::runEstimateParameters()
{
  GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, 1, false, "WlopParaDlg::runEstimateParameters()()");

  double sum_dist = 0;
  for (int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];
    CVertex& t = samples->vert[v.neighbors[0]];

    sum_dist += GlobalFun::computeEulerDist(v.P(), t.P());
  }
  double avg_dist = sum_dist / samples->vert.size();

  global_paraMgr.wLop.setValue("Increasing Step Size", DoubleValue(avg_dist * 0.05));
  global_paraMgr.wLop.setValue("Local Neighbor Size For Inner Points", DoubleValue(avg_dist * 1.5));
  global_paraMgr.wLop.setValue("Local Neighbor Size For Surface Points", DoubleValue(avg_dist * 4.0));


  cout << "Increasing Step Size" << avg_dist * 0.005 << endl;
  cout << "Local Neighbor Size For Inner Points" << avg_dist * 1.5 << endl;
  cout << "Increasing Step Size" << avg_dist * 5.0 << endl;

}

void WLOP::runComputeAverageDistThreshold()
{
	double original_knn = para->getDouble("Original Averaging KNN");
	GlobalFun::computeAnnNeigbhors(original->vert, samples->vert, original_knn, false, "WlopParaDlg::runRegularizeNormals()");

	vector<double> avg_dists(samples->vert.size());
	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];

		double sum_dist = 0;
		double sum_weight = 0;
		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& t = original->vert[v.neighbors[j]];
			sum_dist += GlobalFun::computeEulerDist(v.P(), t.P());
			sum_weight += 1.0;
		}

		//v.eigen_confidence = sum_dist / sum_weight;
		v.nearest_neighbor_dist = sum_dist / sum_weight;
		avg_dists[i] = v.nearest_neighbor_dist;
	}

	std::sort(avg_dists.begin(), avg_dists.end());
	double percentage = para->getDouble("Choose ADT Threshold Percentage");
	int index = avg_dists.size() * percentage;
	double threshold = avg_dists[index];

	para->setValue("Average Dist To Input Threshold", DoubleValue(threshold));

	cout << "dist threshold!  " << para->getDouble("Average Dist To Input Threshold") << endl;
}

void WLOP::runComputeConfidence()
{
 	for (int i = 0; i < samples->vert.size(); i++)
 	{
 		CVertex& v = samples->vert[i];
 		v.eigen_confidence = 0.0;
 	}
// 	return;
	cout << "comput confidence" << endl;

	double original_knn = para->getDouble("Original Confidence KNN");

	GlobalFun::computeAnnNeigbhors(original->vert, samples->vert, original_knn, false, "WlopParaDlg::runComputeConfidence()");
	
	bool use_dist_threshold = para->getBool("Use Average Dist Threshold");
	//bool use_dist_threshold = use_ellipsoid_weight;
	double theshold = para->getDouble("Average Dist To Input Threshold");
	double sigma = 45;
	double sigma_threshold = pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2);

	double radius = para->getDouble("CGrid Radius");	
	double radius2 = radius * radius;
	double iradius16 = -4.0 / radius2;


	//ofstream outfile("average_neighbor_dist.txt");
	int cnt_hole = 0;
	double min_dist = 10000;
	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];

		double sum_dist = 0;
		double sum_weight = 0;
		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& t = original->vert[v.neighbors[j]];
			sum_dist += GlobalFun::computeEulerDist(v.P(), t.P())/* * normal_diff*/;
			sum_weight += 1.0;
		}

		if (sum_weight > 0)
		{
			v.eigen_confidence = sum_dist / sum_weight;
			if (min_dist > v.eigen_confidence)
			{
				min_dist = v.eigen_confidence;
			}
		}
		else
		{
			cout << "NO neighbor???" << endl;
			v.eigen_confidence = 1;
		}
		v.nearest_neighbor_dist = v.eigen_confidence;
		//outfile << v.nearest_neighbor_dist << endl;
	}
	cout << "how many hole points in compute confidence?: " << cnt_hole << endl;



	if (use_dist_threshold)
	{
		for (int i = 0; i < samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];

			if (v.eigen_confidence < theshold)
			{
				v.eigen_confidence = min_dist;
			}
			else
			{
				cnt_hole++;
			}
		}
	}

	GlobalFun::normalizeConfidence(samples->vert, 0);
	double confidenc_power = para->getDouble("Confidence Power");

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		//v.neighbors.clear();
		//v.eigen_confidence = 1 - v.eigen_confidence;

		if (use_dist_threshold)
		{
			v.eigen_confidence = 1 - v.eigen_confidence;
			//v.eigen_confidence = pow(1 - v.eigen_confidence, 2);
		}
		else
		{
			v.eigen_confidence = pow(1 - v.eigen_confidence, confidenc_power);
		}
		
	}
	GlobalFun::normalizeConfidence(samples->vert, 0);

 	if (use_dist_threshold)
 	{
 		GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius") * 2.0, samples->bbox);
 		GlobalFun::smoothConfidences(samples, para->getDouble("CGrid Radius"));
 		GlobalFun::normalizeConfidence(samples->vert, 0);
 	}

	//ofstream outfile2("average_neighbor_dist2.txt");
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		CVertex& v = samples->vert[i];
// 		//outfile2 << v.eigen_confidence << endl;
// 	}


	cout << "finshed compute confidence#######" << endl;

}


// this one is useful
//void WLOP::runComputeConfidence()
//{
////  	//double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
////  	int knn = global_paraMgr.norSmooth.getInt("PCA KNN");
////  	double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
////  
////  	double sigma_threshold = pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2);
////  
////  	GlobalFun::computeAnnNeigbhors(original->vert, samples->vert, knn, false, "WlopParaDlg::runRegularizeNormals()");
////  
////  	for (int i = 0; i < samples->vert.size(); i++)
////  	{
////  		CVertex& v = samples->vert[i];
////  
////  		v.eigen_confidence = 0.0;
////  		for (int j = 0; j < v.neighbors.size(); j++)
////  		{
////  			int neighbor_idx = v.neighbors[0];
////  			CVertex& t = original->vert[neighbor_idx];
////  			double normal_diff = exp(-pow(1 - v.N() * t.N(), 2) / sigma_threshold);
////  
////  			v.eigen_confidence += normal_diff;
////  		}
////  		v.eigen_confidence /= v.neighbors.size();
////  	}
////  
////  	//GlobalFun::normalizeConfidence(samples->vert, 0.0);
////  
////  	return;
//
//
//
//// 	Timer time;
//// 	time.start("Samples Initial");
//// 	GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
//// 	GlobalFun::computeEigenWithTheta(dual_samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));
//// 	time.end();
//// 
//// 	GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 1, false, "WlopParaDlg::runRegularizeNormals()");
//// 
//// 	for (int i = 0; i < samples->vert.size(); i++)
//// 	{
//// 		CVertex& v = samples->vert[i];
//// 
//// 		int neighbor_idx = v.neighbors[0];
//// 
//// 		CVertex& dual_v = dual_samples->vert[neighbor_idx];
//// 		v.eigen_confidence = dual_v.eigen_confidence;
//// 	}
//// 
//// 	GlobalFun::normalizeConfidence(samples->vert, 0);
//// 
//// 	return;
//
//
//
//	// 	cout << "compute confidence" << endl;
//	// 	cout << samples->vert.size() << "	" << original->vert.size() << endl;
//	// 	return;
//	double radius = para->getDouble("CGrid Radius");
//	
//	double radius2 = radius * radius;
//	double iradius16 = -4.0 / radius2;
//
//	double original_knn = para->getDouble("Original Confidence KNN");
//
//	GlobalFun::computeAnnNeigbhors(original->vert, samples->vert, original_knn, false, "WlopParaDlg::runRegularizeNormals()");
//
//	//ofstream outfile("nearest_neighbor_dist.txt");
//
//	for (int i = 0; i < samples->vert.size(); i++)
//	{
//		CVertex& v = samples->vert[i];
//
//		double sum_dist = 0;
//		for (int j = 0; j < v.neighbors.size(); j++)
//		{
//			CVertex& t = original->vert[v.neighbors[j]];
//			sum_dist += GlobalFun::computeEulerDist(v.P(), t.P());
//		}
//
//		v.eigen_confidence = sum_dist / v.neighbors.size();
//
//		//outfile << v.nearest_neighbor_dist << endl;
//	}
//
//
//	GlobalFun::normalizeConfidence(samples->vert, 0);
//
//	for (int i = 0; i < samples->vert.size(); i++)
//	{
//		CVertex& v = samples->vert[i];
//		v.neighbors.clear();
//		//v.eigen_confidence = 1 - v.eigen_confidence;
//
//		v.eigen_confidence = pow(1-v.eigen_confidence, 2);
//	}
//	GlobalFun::normalizeConfidence(samples->vert, 0);
//
// 
////    	GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
////    	vector<double> sum_confidences(samples->vn, 1e-15);
////  		vector<double> sum_weights(samples->vn, 1e-15);
////   
////    	for (int i = 0; i < samples->vert.size(); i++)
////    	{
////    		CVertex& v = samples->vert[i];
////    		
////    		double sum_confidence = 0;
////    		double weight_sum = 0;
////    		for (int j = 0; j < v.neighbors.size(); j++)
////    		{
////    			CVertex& t = samples->vert[v.neighbors[j]];
////    			double dist2 = GlobalFun::computeEulerDistSquare(v.P(), t.P());
////    
////    			double dist_diff = exp(dist2 * iradius16);
////    			double w = dist_diff ;
////    
////    			sum_confidence += w * t.eigen_confidence;
////    			weight_sum += w;
////    		}
////    
////    		sum_confidences[i] = sum_confidence;
////    		sum_weights[i] = weight_sum;
////    	}
////   
////    	for (int i = 0; i < samples->vert.size(); i++)
////    	{
////    		CVertex& v = samples->vert[i];
////    
////  			//v.eigen_confidence = 0.0;
////    		if (sum_weights[i] > 1e-8)
////    		{
////    			v.eigen_confidence = sum_confidences[i] / sum_weights[i];
////    		}
////    	}
//
//
//	// 	double threshold = para->getDouble("sigmoid threshold");
//	// 	//double threshold_shift = 0.5 - threshold;
//	// 	for (int i = 0; i < samples->vert.size(); i++)
//	// 	{
//	// 		CVertex& v = samples->vert[i];
//	// 		double cofidence = (v.eigen_confidence - threshold) * 10;
//	// 		v.eigen_confidence = 1.0 / (1.0 + exp(-cofidence));
//	// 		//v.eigen_confidence = 1.0 / (1.0 + exp(-v.eigen_confidence));
//	// 
//	// 	}
//
//// 	ofstream out_file("confidence.txt");
//// 	for (int i = 0; i < samples->vert.size(); i++)
//// 	{
//// 		out_file << samples->vert[i].eigen_confidence << endl;
//// 	}
//// 	out_file.close();
//// 
//// 
// 	cout << "finshed compute confidence#######" << endl;
//
//	// 	ofstream out_file2("dual_vector_length.txt");
//	// 	for (int i = 0; i < samples->vert.size(); i++)
//	// 	{
//	// 		CVertex& v = samples->vert[i];
//	// 		CVertex& dual_v = dual_samples->vert[i];
//	// 
//	// 		double length = GlobalFun::computeEulerDist(v.P(), dual_v.P());
//	// 		out_file2 << length << endl;
//	// 	}
//	// 	out_file2.close();
//
//}













void WLOP::smoothSkelDistance()
{
	initVertexes(true);

	computeSampleSimilarityTerm(samples);
	for (int i = 0; i < samples_similarity.size(); i++)
	{
		CVertex& v = samples->vert[i];
		v.P() = samples_similarity[i];
	}

}


void WLOP::computeInitialSampleNeighbor()
{
  
  Timer time;
  time.start("Sample Sample Neighbor Tree");
  GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
  time.end();

}

void WLOP::computeNearestNeighborDist()
{
	GlobalFun::computeAnnNeigbhors(original->vert, samples->vert, 1, false, "WlopParaDlg::runRegularizeNormals()");

	//ofstream outfile("nearest_neighbor_dist.txt");

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		CVertex& t = original->vert[v.neighbors[0]];

		v.nearest_neighbor_dist = GlobalFun::computeEulerDist(v.P(), t.P());
	}

	//outfile.close();
}



void WLOP::runSkelWlop()
{
  cout << "runSkelWlop" << endl;

  Timer time;

  initVertexes(true);

	if (use_eigen_neighborhood)
	{
		cout << "use_eigen_neighborhood" << endl;
		time.start("Samples runComputeEigenNeighborhood");
		runComputeEigenNeighborhood(samples, original);

		if (nTimeIterated == 0)
		{
			GlobalFun::computeBallNeighbors(original, NULL,
				para->getDouble("CGrid Radius"), original->bbox);

			//cout << endl  << "compute original_density original_density original_density" << endl;
			original_density.assign(original->vn, 0);
			if (para->getBool("Need Compute Density"))
			{
				cout << "compute density!!" << endl;
				computeDensity(true, para->getDouble("CGrid Radius") /**0.5*/, samples, original);
			}
		}
		time.end();
	}
	else
	{
		time.start("Samples Initial");
		GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
		GlobalFun::computeEigenWithTheta(samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));
		time.end();

		if (nTimeIterated == 0)
		{
			time.start("Original Initial");
			GlobalFun::computeBallNeighbors(original, NULL,
				para->getDouble("CGrid Radius"), original->bbox);

			original_density.assign(original->vn, 0);
			if (para->getBool("Need Compute Density"))
			{
				computeDensity(true, para->getDouble("CGrid Radius"), samples, original);
			}
			time.end();
		}

		time.start("Sample Original neighbor");
		GlobalFun::computeBallNeighbors(samples, original,
			para->getDouble("CGrid Radius"), box);
		time.end();
	}
	

	time.start("computeAverageTerm");
	computeAverageTerm(samples, original);
	time.end();


	if (!use_eigen_neighborhood)
	{
		GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
		GlobalFun::computeEigenWithTheta(samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));
	}
	else
	{
		time.start("Samples runComputeEigenNeighborhood");
		runComputeEigenNeighborhood(samples, original);
		time.end();
	}

	if (para->getBool("Need Compute Density"))
	{
		
		time.start("Compute Density For Sample");
		computeDensity(false, para->getDouble("CGrid Radius") *0.5, samples, original);
		time.end();
	}

	time.start("computeRepulsionTerm");
	computeRepulsionTerm(samples);
	time.end();


  double min_sigma = GlobalFun::getDoubleMAXIMUM();
  double max_sigma = -1;
  for (int i = 0; i < samples->vn; i++)
  {
    CVertex& v = samples->vert[i];
    if (v.eigen_confidence < min_sigma)
    {
      min_sigma = v.eigen_confidence;
    }
    if (v.eigen_confidence > max_sigma)
    {
      max_sigma = v.eigen_confidence;
    }
  }

//  double mu_max = global_paraMgr.skeleton.getDouble("Repulsion Mu");
// 	double mu_min = global_paraMgr.skeleton.getDouble("Repulsion Mu2");
//   double mu_length = abs(mu_max - mu_min);
//  double sigma_length = abs(max_sigma - min_sigma);
	double mu = para->getDouble("Repulsion Mu");
  Point3f c;
  int moving_num = 0;
  double max_error = 0;

  for(int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];
    if (v.is_fixed_sample || v.is_skel_ignore)
    {
      continue;
    }
    c = v.P();


		if (mu > 1.0)
		{
			if (i <5)
			{
				cout << "no using first term" << endl;
			}
		}
    else if (average_weight_sum[i] > 1e-20)
    {
			if (use_kite_points)
			{
				if (i <3)
				{
					cout << "use_kite_points use_kite_points" << endl;
				}
				CVertex& sample_point = original->vert[i];
				Point3f avg_point = average[i] / average_weight_sum[i];
				double proj_dist = abs((avg_point - sample_point.P()) * sample_point.N());
				Point3f new_point = sample_point.P() - sample_point.N() * proj_dist;
				v.P() = new_point;
			}
			else
			{
				v.P() = average[i] / average_weight_sum[i];
			}
    }

		if (repulsion_weight_sum[i] > 1e-20 && mu > 0)
		{
 			if (use_ellipsoid_repulsion )
 			{
				repulsion_x_length[i] = repulsion_x_length[i] * (mu / repulsion_weight_sum[i]);
				repulsion_y_length[i] = repulsion_y_length[i] * (mu / repulsion_weight_sum[i]);
				repulsion_z_length[i] = repulsion_z_length[i] * (mu / repulsion_weight_sum[i]);

				Point3f final_dir = Point3f(0, 0, 0);
				double threshold = 0.3;
				if (v.eigen_value0 > threshold)
				{
					final_dir  += v.eigen_vector0 * v.eigen_value0 * repulsion_x_length[i];
				}
				if (v.eigen_value1 > threshold)
				{
					final_dir += v.eigen_vector1 * v.eigen_value1 * repulsion_y_length[i];
				}
				if (v.eigen_value2 > threshold)
				{
					final_dir += v.eigen_vector2 * v.eigen_value2 * repulsion_z_length[i];
				}
			
				if (i < 5)
				{
					cout << "eigen values final: " << v.eigen_value0 << "	" << v.eigen_value1 << "	" << v.eigen_value2 << "	" << endl;
				}

				v.P() += final_dir;
				//v.N() = final_dir.Normalize();
 			}
 			else
			{
				v.P() += repulsion[i] * (mu / repulsion_weight_sum[i]);
				//v.N() = repulsion[i].Normalize();
			}
		}

    if (average_weight_sum[i] > 1e-20 && repulsion_weight_sum[i] > 1e-20 )
    {
      moving_num++;
      Point3f diff = v.P() - c; 
      double move_error = sqrt(diff.SquaredNorm());

      error_x += move_error; 
    }
  }
  error_x = error_x / moving_num;

  para->setValue("Current Movement Error", DoubleValue(error_x));
  cout << "****finished compute Skeletonization error:	" << error_x << endl;
  return ;
}





void WLOP::runRegularizeSamples()
{
  cout << "runRegularizeSamples" << endl;

	if (global_paraMgr.glarea.getBool("Show Dual Samples"))
	{
		samples = dual_samples;
	}

  GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
  GlobalFun::computeEigenWithTheta(samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));

  int branch_KNN = global_paraMgr.skeleton.getDouble("Sigma KNN");
   //int branch_KNN = 20;
   GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, branch_KNN, false, "void Skeletonization::searchNewBranches()");
  
  double feature_sigma = global_paraMgr.skeleton.getDouble("Eigen Feature Identification Threshold");
  cout << "feature sigma:" << feature_sigma << endl;

  vector<Point3f> new_sample_set;
  new_sample_set.assign(samples->vert.size(), Point3f(0, 0, 0));
  for (int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];
    Point3f direction = v.eigen_vector0.Normalize();

    if (v.eigen_confidence < feature_sigma)
    {
      continue;;
    }

    Point3f front_nearest_p = v.P();
    Point3f back_nearest_p = v.P();
    double front_nearest_dist = 1000000.;
    double back_nearest_dist = -1000000.;

    for (int j = 0; j < v.neighbors.size(); j++)
    {
      int neighbor_idx = v.neighbors[j];
      CVertex& t = samples->vert[neighbor_idx];

      Point3f diff = t.P() - v.P();
      double proj_dist = diff * direction;

      if (proj_dist > 0)
      {
        if (proj_dist < front_nearest_dist)
        {
          front_nearest_p = t.P();
          front_nearest_dist = proj_dist;
        }
      }
      else
      {
        if (proj_dist > back_nearest_dist)
        {
          back_nearest_p = t.P();
          back_nearest_dist = proj_dist;
        }
      }
    }

    if (front_nearest_dist > 100000 || back_nearest_dist < -100000)
    {
      cout << "end pointssss" << endl;
      continue;
    }

    v.P() = (front_nearest_p + back_nearest_p) / 2.0;
  }
}



void WLOP::runRegularizeNormals(CMesh* samples, CMesh* dual_samples)
{
	computeDualIndex(samples, dual_samples);

 	GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
 	GlobalFun::computeEigenWithTheta(dual_samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));

	bool use_confidence = para->getBool("Use Confidence");
	bool use_confidence_to_combine = para->getBool("Use Confidence To Combine Normal");


	double confidence_threshold = para->getDouble("Density Confidence Segment Threshold");
	double eigen_directional_threshold = para->getDouble("Eigen Directional Threshold");
	int hole_points = 0;
	int tubular_points = 0;

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];

		hole_points++;

		CVertex& dual_v = dual_samples->vert[v.dual_index];
		Point3f diff = v.P() - dual_v.P();

		Point3f dir = diff;


    //这里应该用另外的判断法则：假设骨架的洞都被填上
		if (use_confidence && v.eigen_confidence > confidence_threshold)
		{
			if (dir.Normalize() * v.N() < 0)
			{
				v.N() *= -1;
			}

			continue;
		}

		if (use_confidence && dual_v.eigen_confidence < eigen_directional_threshold)
		{
			if (dir.Normalize() * v.N() < 0)
			{
				v.N() *= -1;
			}

			continue;
		}

		dir = diff.Normalize();

		if (use_confidence_to_combine)
		{
			if (dir*v.N() < 0 && v.eigen_confidence <= 1.1)
			{
				v.N() = ((dir * (1 - v.eigen_confidence) + (-v.N()) * v.eigen_confidence)).Normalize();
			}
			else
			{
				v.N() = ((dir * (1-v.eigen_confidence) + v.N() * v.eigen_confidence)).Normalize();
			}
		}
		else
		{
			if (dir*v.N() < 0)
			{
				v.N() = ((dir + (-v.N())) / 2.0).Normalize();
			}
			else
			{
				v.N() = ((dir + v.N()) / 2.0).Normalize();
			}
		}

	}

	cout << "!!how many hole points in runRegularizeNormals? " << hole_points << endl;
	cout << "!!how many tubular points in runRegularizeNormals? " << tubular_points << endl;

	return;
}





void WLOP::runComputeCorrespondence()
{
	vector<SphereSlots> source_distributions = computeDistributions(samples, dual_samples);
	vector<SphereSlots> target_distributions = computeDistributions(target_samples, target_dual_samples);

	Timer time;
	time.start("correspondence");
	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];

		SphereSlots& src_slots = source_distributions[i];
		double min_dist = GlobalFun::getDoubleMAXIMUM();
		int min_index = -1;
		for (int j = 0; j < target_distributions.size(); j++)
		{
			SphereSlots& target_slots = target_distributions[j];
			double dist = src_slots.computeDistance(target_slots);
			//double dist = src_slots.computeL1Distance(target_slots);

			if (dist < min_dist)
			{
				min_dist = dist;
				min_index = j;
			}
		}

		if (min_index >= 0)
		{
			v.target_index = min_index;
			v.eigen_confidence = min_dist;
		}
	}
	time.end();

	GlobalFun::normalizeConfidence(samples->vert, 0);
	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		v.eigen_confidence = 1 - v.eigen_confidence;
	}
}

vector<SphereSlots> WLOP::computeDistributions(CMesh* samples, CMesh* dual_samples)
{
	cout << "runComputeDistribution" << endl;

	int knn = global_paraMgr.norSmooth.getInt("PCA KNN");
	cout << "knn:	" << knn << endl;
	GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
	//GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
	//GlobalFun::computeAnnNeigbhors(dual_samples->vert, dual_samples->vert, knn, false, "void runComputeDistribution()");

	if (use_closest_dual)
	{
		GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 1, false, "WlopParaDlg::runRegularizeNormals()");
	}
	//return;

	double radius = para->getDouble("CGrid Radius");
	double radius2 = radius * radius;
	double iradius16 = -para->getDouble("H Gaussian Para") / radius2;

	vector<SphereSlots> sphere_distributions(dual_samples->vert.size(), default_sphere);
	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex dual_v = dual_samples->vert[i];
		CVertex& v = samples->vert[i];

		if (use_closest_dual)
		{
			int neighbor_idx = v.neighbors[0];
			dual_v = dual_samples->vert[neighbor_idx];
		}
		//Point3f direction = (v.P() - dual_v.P()).Normalize();
		Point3f direction = v.P() - dual_v.P();

		SphereSlots& sphere_slots = sphere_distributions[i];
		//sphere_slots.rearrangeSlots(direction.Normalize());
		
		updateSphereSlots(sphere_slots, direction, 1.0);

		for (int j = 0; j < dual_v.neighbors.size(); j++)
		{
			int dual_neighbor_idx = dual_v.neighbors[j];
			CVertex dual_t = dual_samples->vert[dual_neighbor_idx];
			CVertex& dual_t_v = samples->vert[dual_neighbor_idx];

			double dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_t_v.P());
			double dist_diff = exp(dist2 * iradius16);

			if (use_closest_dual)
			{
				int neighbor_idx = dual_t_v.neighbors[0];
				dual_t = dual_samples->vert[neighbor_idx];
			}

			Point3f direction = (dual_t_v.P() - dual_t.P()).Normalize();
			updateSphereSlots(sphere_slots, direction, dist_diff);
		}

		//v.eigen_confidence = getSphereSlotsConfidence(sphere_slots);
	}

	//GlobalFun::normalizeConfidence(samples->vert, 0.0);

	return sphere_distributions;
}


void WLOP::runComputeDistribution()
{
	if (global_paraMgr.glarea.getBool("Show Skeleton"))
	{
		bool use_cloest = global_paraMgr.glarea.getBool("Show Cloest Dual Connection");

		for (int i = 0; i < samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];
			CVertex dual_v = dual_samples->vert[i];

			if (use_closest_dual)
			{
				int neighbor_idx = v.neighbors[0];
				dual_v = dual_samples->vert[neighbor_idx];
			}
			double dist = GlobalFun::computeEulerDist(v.P(), dual_v.P());

			v.eigen_confidence = dist;
		}

		GlobalFun::normalizeConfidence(samples->vert, 0.0);
		return;
	}


	cout << "runComputeDistribution" << endl;

	int knn = global_paraMgr.norSmooth.getInt("PCA KNN");
	cout << "knn:	" << knn <<endl;
	GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
	//GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
	//GlobalFun::computeAnnNeigbhors(dual_samples->vert, dual_samples->vert, knn, false, "void runComputeDistribution()");

	if (use_closest_dual)
	{
		GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 1, false, "");
	}
	//return;

	double radius = para->getDouble("CGrid Radius");
	double radius2 = radius * radius;
	double iradius16 = -para->getDouble("H Gaussian Para") / radius2;

	vector<SphereSlots> sphere_distributions(dual_samples->vert.size(), default_sphere);
	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex dual_v = dual_samples->vert[i];
		CVertex& v = samples->vert[i];

		if (use_closest_dual)
		{
			int neighbor_idx = v.neighbors[0];
			dual_v = dual_samples->vert[neighbor_idx];
		}
		//Point3f direction = (v.P() - dual_v.P()).Normalize();
		Point3f direction = v.P() - dual_v.P();

		SphereSlots& sphere_slots = sphere_distributions[i];
		updateSphereSlots(sphere_slots, direction, 1.0);

		for (int j = 0; j < dual_v.neighbors.size(); j++)
		{
			int dual_neighbor_idx = dual_v.neighbors[j];
			CVertex dual_t = dual_samples->vert[dual_neighbor_idx];
			CVertex& dual_t_v = samples->vert[dual_neighbor_idx];

			double dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_t_v.P());
			double dist_diff = exp(dist2 * iradius16);

			if (use_closest_dual)
			{
				int neighbor_idx = dual_t_v.neighbors[0];
				dual_t = dual_samples->vert[neighbor_idx];
			}

			Point3f direction = (dual_t_v.P() - dual_t.P()).Normalize();
			updateSphereSlots(sphere_slots, direction, dist_diff);
		}

		v.eigen_confidence = getSphereSlotsConfidence(sphere_slots);
	}

	GlobalFun::normalizeConfidence(samples->vert, 0.0);

}

void WLOP::updateSphereSlots(SphereSlots& sphere_slots, Point3f dir, double dist_weight)
{
	double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
	double sigma_threshold = pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2);

	double length = sqrt(dir.SquaredNorm());
	dir = dir.Normalize();
	//bool use_length_confidence = para->getBool("Use Confidence");
	bool use_length_confidence = true;

	for (int i = 0; i < sphere_slots.size(); i++)
	{
		SphereSlot& sphere_slot = *sphere_slots.getSlot(i);
		Point3f slot_direction = sphere_slot.slot_direction;

		double direction_diff = exp(-pow(1 - slot_direction * dir, 2) / sigma_threshold);
		
		if (use_length_confidence)
		{
			sphere_slot.slot_value += direction_diff * length * dist_weight;
		}
		else
		{
			sphere_slot.slot_value += direction_diff * dist_weight;
		}
	}
}

double WLOP::getSphereSlotsConfidence(SphereSlots& sphere_slots)
{
	double sum_confidence = 0.0;

	bool use_length_confidence = para->getBool("Use Confidence");

	for (int i = 0; i < sphere_slots.size(); i++)
	{
		SphereSlot& sphere_slot = *sphere_slots.getSlot(i);

		double confidence = sphere_slot.slot_value < 2.0 ? sphere_slot.slot_value : 1.0;
		
		if (use_length_confidence)
		{
			confidence = sphere_slot.slot_value;
		}

		//double confidence = sphere_slot.slot_value;
		sum_confidence += confidence;
	}
	return sum_confidence;
}






void WLOP::runProjection()
{
	double radius = para->getDouble("CGrid Radius");
	double percentage = para->getDouble("Repulsion Mu") * 0.1;
	GlobalFun::addOutliers(samples, percentage, radius);
}




void WLOP::runMoveSample()
{
	if (true)
	{
// 		GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
// 		GlobalFun::computeEigenWithTheta(dual_samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));

		initVertexes(true);

		GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
		GlobalFun::computeEigenWithTheta(samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));

		Timer time;
		time.start("Compute Repulsion Term");
		computeRepulsionTerm(samples);
		time.end();

		double mu = 0.4;
		for (int i = 0; i < dual_samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];
			CVertex& dual_v = dual_samples->vert[i];

			Point3f diff = v.P() - dual_v.P();

			double proj_dist = diff * dual_v.N();
			v.P() = dual_v.P() + dual_v.N() * proj_dist;

			if (repulsion_weight_sum[i] > 1e-10)
			{
				v.P() += repulsion[i] * (mu / repulsion_weight_sum[i]);
			}
		}
	}
	else
	{



	}


}

void WLOP::runMoeveSkel()
{
	GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
	GlobalFun::computeEigenWithTheta(dual_samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));

	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		CVertex& dual_v = dual_samples->vert[i];

		Point3f diff = v.P() - dual_v.P();
		//double proj_dist = diff * dual_v.eigen_vector0;
		//v.P() -= dual_v.eigen_vector0 * proj_dist;


		// 		double proj_dist = diff * dual_v.N();
		// 		v.P() = dual_v.P() + dual_v.N() * proj_dist;

		double proj_dist = diff * dual_v.eigen_vector0;
		//dual_v.N() = dual_v.eigen_vector0;
		dual_v.P() += dual_v.eigen_vector0 * proj_dist;

		Point3f new_dir = (v.P() - dual_v.P()).Normalize();
		dual_v.N() = new_dir;
	}

}

void WLOP::runEllipsoidFitting()
{
	GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
	GlobalFun::computeEigenWithTheta(dual_samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));

	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		CVertex& dual_v = dual_samples->vert[i];

		Point3f diff = v.P() - dual_v.P();

		double proj_dist = diff * dual_v.eigen_vector0;
		//dual_v.N() = dual_v.eigen_vector0;
		dual_v.P() += dual_v.eigen_vector0 * proj_dist;
	}

}



void WLOP::runInnerPointsRegularization()
{
	cout << "runInnerPointsRegularization" << endl;

	double max_skel_radius = 0.0;

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];

		if (v.skel_radius > max_skel_radius)
		{
			max_skel_radius = v.skel_radius;
		}
	}
	cout << "Max skel radius:  " << max_skel_radius << endl;

	Timer time;
	time.start("Sample Original neighbor");
	GlobalFun::computeBallNeighbors(dual_samples, samples,
		                              max_skel_radius * 0.5, box);
	time.end();

	vector<double> farthest_proj_dist(dual_samples->vert.size(), 0.0);
	vector<Point3f> data_term(dual_samples->vert.size());
	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];
		CVertex& v = samples->vert[i]; 
		double threshold = v.skel_radius * 0.5;
		double skel_radius2 = threshold * threshold;

		for (int j = 0; j < dual_v.original_neighbors.size(); j++)
		{
			int idx = dual_v.original_neighbors[j];
			CVertex& t = samples->vert[idx];

			double dist2 = GlobalFun::computeEulerDistSquare(v.P(), t.P());
			if (dist2 > skel_radius2)
			{
				continue;
			}

			//double proj_dist = abs((dual_v.P() - t.P()) * dual_v.N());
			double proj_dist = abs((v.P() - t.P()) * dual_v.N());

			if (proj_dist > farthest_proj_dist[i])
			{
				farthest_proj_dist[i] = proj_dist;
			}
		}
	}

	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];
		CVertex& v = samples->vert[i];

		double move_dist = farthest_proj_dist[i] * 0.5;
		data_term[i] = v.P() - v.N() * move_dist;
	}

	// compute repulsion term
	vector<Point3f> regular_term(dual_samples->vert.size());

	repulsion.assign(samples->vn, vcg::Point3f(0, 0, 0));
	repulsion_weight_sum.assign(samples->vn, 0);

	//bool need_density = para->getBool("Need Compute Density");
	bool need_density = false;

	double radius = para->getDouble("CGrid Radius");
	double radius2 = radius * radius;
	double iradius16 = -para->getDouble("H Gaussian Para") / radius2;
	double repulsion_power = para->getDouble("Repulsion Power");


	//Timer time;
	time.start("Sample Original neighbor");
	GlobalFun::computeBallNeighbors(dual_samples, NULL,	radius, box);
	time.end();

	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];

		if (dual_v.neighbors.empty())
		{
			continue;
		}

		for (int j = 0; j < dual_v.neighbors.size(); j++)
		{
			CVertex& dual_t = dual_samples->vert[dual_v.neighbors[j]];
			Point3f diff = dual_v.P() - dual_t.P();

			double dist2 = diff.SquaredNorm();
			double len = sqrt(dist2);
			if (len <= 0.001 * radius) len = radius*0.001;

			double w = exp(dist2*iradius16);
			double rep = w * pow(1.0 / len, repulsion_power);

// 			if (need_density)
// 			{
// 				rep *= samples_density[dual_t.m_index];
// 			}

			repulsion[i] += diff * rep;
			repulsion_weight_sum[i] += rep;
		}

		if (repulsion_weight_sum[i] > 1e-10)
		{
			regular_term[i] = repulsion[i] / repulsion_weight_sum[i];
		}
	}

	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];

		dual_v.P() = data_term[i] + regular_term[i] * 0.25;
	}
}

void WLOP::runSearchNeighborhood()
{
	cout << "runSearchNeighborhood" << endl;

}

void WLOP::runSmoothNeighborhood()
{
	//double radius2 = para->getDouble("CGrid Radius");
	GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
	GlobalFun::smoothConfidences(samples, para->getDouble("CGrid Radius"));
	GlobalFun::normalizeConfidence(samples->vert, 0.0);

	return;

	cout << "runSmoothNeighborhood" << endl;
	double local_radius = para->getDouble("CGrid Radius");

	//double local_radius = 0.06;
	
	double radius = local_radius;
	double radius2 = radius * radius;
	double iradius16 = -4 / radius2;

	Timer time;
	time.start("Sample Original neighbor");
	GlobalFun::computeBallNeighbors(samples, NULL,
		para->getDouble("CGrid Radius") , box);
	time.end();

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		v.eigen_confidence = v.skel_radius;
	}

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];

		if (v.eigen_confidence > 1e-10)
		{
			continue;
		}

		cout << "new 2eigen_confidence: ";

		if (v.neighbors.empty())
		{
			continue;
		}

		
		vector<int>* neighbors = &v.neighbors;

		if (v.neighbors.empty())
		{
			continue;
		}

		v.eigen_confidence = 0.0;

		double sum_confidence = 0;
		double weight_sum = 0;

		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& t = samples->vert[(*neighbors)[j]];

			if (t.skel_radius < 1e-10)
			{
				continue;
			}
			double dist2 = (v.P() - t.P()).SquaredNorm();
			double w = exp(dist2*iradius16);

			sum_confidence += w * t.eigen_confidence;
			weight_sum += w;
		}

		if (weight_sum >= 1e-10)
		{
			v.eigen_confidence += sum_confidence / weight_sum;
		}

		cout << v.eigen_confidence << endl;
	}

	//normalizeConfidence(mesh->vert, 0.0);
// 	vector<double> confidences(samples->vert.size());
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		confidences[i] = samples->vert[i].eigen_confidence;
// 	}




	GlobalFun::smoothConfidences(samples, local_radius);

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		v.skel_radius = v.eigen_confidence;
	}
}


void WLOP::computeConstNeighborhoodUsingKNN(int knn)
{
	GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, knn, false, "WlopParaDlg::runMoveBackward()");
	cout << "end knn" << endl;
}

void WLOP::computeConstNeighborhoodUsingRadius(double radius)
{
	GlobalFun::computeBallNeighbors(samples, NULL, radius, box);
}

void WLOP::runMoveBackward()
{
	double curr_radius = global_paraMgr.wLop.getDouble("CGrid Radius");

	double stop_angle_threshold = para->getDouble("Local Angle Threshold");
	double stop_neighbor_size = para->getDouble("Local Neighbor Size For Inner Points");
	double step_size = para->getDouble("Increasing Step Size");
	double cooling_parameter = para->getDouble("Inner Points Cooling Parameter");
	Timer time;

	initVertexes(true);

  time.start("Sample@@@@dual_samples neighbor");
  GlobalFun::computeBallNeighbors(dual_samples, samples,
    stop_neighbor_size, box);
  time.end();


	time.start("Sample0000dual_samples neighbor");
	GlobalFun::computeBallNeighbors(dual_samples, NULL,
		stop_neighbor_size, box);
	time.end();


	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];

		if (dual_v.is_fixed_sample)
		{
			continue;
		}

		bool keep_going = true;
		for (int j = 0; j < dual_v.neighbors.size(); j++)
		{
			CVertex& t = dual_samples->vert[dual_v.neighbors[j]];

			double angle = GlobalFun::computeRealAngleOfTwoVertor(dual_v.N(), t.N());

			if (angle > stop_angle_threshold )
			{
				keep_going = false;
				break;
			}
		}

		if (!keep_going)
		{
			dual_v.is_fixed_sample = true;
      dual_v.moving_speed = 0.0;

			if (dual_v.skel_radius < 1e-15)
			{
				dual_v.skel_radius = curr_radius;
			}
		}
	}


	// compute move step
	//int knn = para->getDouble("Sefl KNN");
	//GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, knn, false, "WlopParaDlg::runMoveBackward()");
	//cout << "end knn" << endl;
	double local_radius = para->getDouble("Local Neighbor Size For Surface Points");
	computeConstNeighborhoodUsingRadius(local_radius);
  double radius2 = local_radius * local_radius;
  double iradius16 = -4 / radius2;

	vector<double> real_step_size(samples->vert.size(), 0.0);
  vector<double> new_speed(samples->vert.size(), 0.0);

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];
		CVertex& v = samples->vert[i];

		if (v.neighbors.empty())
		{
			continue;
		}
    
    double sum_weight = 0.;
    double sum_speed = 0.;
    double weight = 0;

		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& dual_t = dual_samples->vert[v.neighbors[j]];
      CVertex& t = samples->vert[v.neighbors[j]];

      Point3f diff = t.P() - v.P();
      double dist2 = diff.SquaredNorm();
      weight = exp(dist2 * iradius16);

      sum_speed += dual_t.moving_speed * weight;
      sum_weight += weight;
		}

    if (weight > 0)
    {
      new_speed[i] = sum_speed / sum_weight;
    }
	}

	for (int i = 0; i < samples->vert.size(); i++)
	{
    CVertex& dual_v = dual_samples->vert[i];
    dual_v.moving_speed = new_speed[i];

		real_step_size[i] = step_size * dual_v.moving_speed;
	}

// 	// bi-laplician begin
// 	vector<double> real_step_size2(samples->vert.size(), 0.0);
// 
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		CVertex& dual_v = dual_samples->vert[i];
// 		CVertex& v = samples->vert[i];
// 
// 		v.m_index = i;
// 		dual_v.m_index = i;
// 
// 		for (int j = 0; j < v.neighbors.size(); j++)
// 		{
// 			CVertex& dual_t = dual_samples->vert[v.neighbors[j]];
// 
// 			real_step_size2[i] += real_step_size[i];
// 		}
// 
// 		real_step_size2[i] /= v.neighbors.size();
// 	}
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		real_step_size[i] = real_step_size2[i];
// 	}
// 	// bi-laplician end

	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];

		if (dual_v.is_fixed_sample)
		{
			continue;
		}

		dual_v.P() -= dual_v.N() * real_step_size[i];
	}

	return;


 	time.start("Sample0000dual_samples neighbor");
 	GlobalFun::computeBallNeighbors(dual_samples, NULL,
 		stop_neighbor_size, box);
 	time.end();

// 	cout << "!!!time.start(Sample0000dual_samples neighbor)" << endl;
// 	return;
// 	cout << "@@time.start(Sample0000dual_samples neighbor)" << endl;


	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];
    CVertex& v = samples->vert[i];

		if (dual_v.is_fixed_sample || dual_v.neighbors.empty())
		{
			continue;
		}

		if ((real_step_size[i] / step_size) < cooling_parameter)
		{
			dual_v.is_fixed_sample = true;
			//dual_v.is_boundary = true;
			continue;
		}

		bool keep_going = true;
		for (int j = 0; j < dual_v.neighbors.size(); j++)
		{
			CVertex& dual_t = dual_samples->vert[dual_v.neighbors[j]];
      CVertex& t = samples->vert[dual_v.neighbors[j]];

      double dist_v_t = GlobalFun::computeEulerDist(v.P(), t.P());
      if (dist_v_t < stop_neighbor_size * 2.0)
      {
        //cout << "no local static can stop" << endl;
        continue;
      }

      double angle = GlobalFun::computeRealAngleOfTwoVertor(dual_v.N(), dual_t.N());

			if (angle > stop_angle_threshold /*|| dual_v.N() * t.N() < 0*/)
			{
        cout << "angel: " << angle << endl;
				keep_going = false;
				break;
			}
		}

		if (keep_going)
		{
			dual_v.is_fixed_sample = false;
		}
		else
		{
			dual_v.is_fixed_sample = true;


			if (dual_v.skel_radius < 1e-15)
			{
				dual_v.skel_radius = curr_radius;
			}
		}
	}

	//GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, knn, false, "WlopParaDlg::runMoveBackward()");
	computeConstNeighborhoodUsingRadius(local_radius);
}


void WLOP::runSelfWLOP()
{
	use_adaptive_mu = para->getBool("Use Adaptive Mu");
	is_sample_close_to_original.assign(samples->vert.size(), false);
	bool use_tangent = para->getBool("Use Tangent Vector");

	Timer time;
	original = samples;

	initVertexes(true);


	time.start("Sample Sample Neighbor Tree");
	int knn = para->getDouble("Sefl KNN");

 	GlobalFun::computeBallNeighbors(dual_samples, NULL,
		para->getDouble("CGrid Radius"), dual_samples->bbox);

//	GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, knn, false, "WlopParaDlg::runNormalSmoothing()");

	time.end();

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		CVertex& dual_v = dual_samples->vert[i];

		dual_v.original_neighbors = dual_v.neighbors;
	}

	time.start("Compute Average Term");
	computeAverageTerm(dual_samples, samples);
	time.end();


	time.start("Compute Repulsion Term");
	computeRepulsionTerm(dual_samples);
	time.end();

	vector<Point3f> new_sample_positions;
	vector<float> move_proj_vec;

	double radius = para->getDouble("CGrid Radius");
	double radius2 = radius * radius;
	double iradius16 = -para->getDouble("H Gaussian Para") / radius2;

	bool use_confidence = para->getBool("Use Confidence");
	if (use_confidence)
	{
		GlobalFun::normalizeConfidence(samples->vert, 0);
	}

	int error_x = 0;

	samples = dual_samples;
	vector<Point3f> new_positions = computeNewSamplePositions(error_x);
	for (int i = 0; i < new_positions.size(); i++)
	{
		samples->vert[i].P() = new_positions[i];
	}
	

	para->setValue("Current Movement Error", DoubleValue(error_x));
	cout << "****finished compute WLOP error:	" << error_x << endl;

	if (para->getBool("Need Compute PCA"))
	{

	}

	return ;
}


void WLOP::runNormalSmoothing()
{
	cout << "Normal smoothing" << endl;

	double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
	double radius = para->getDouble("CGrid Radius");

	double radius2 = radius * radius;
	double iradius16 = -4 / radius2;
// 	int knn = para->getDouble("Sefl KNN");
// 
// 	//GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
// 	//GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
// 	GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, knn, false, "WlopParaDlg::runNormalSmoothing()");
	double local_radius = para->getDouble("Local Neighbor Size For Surface Points");
	computeConstNeighborhoodUsingRadius(local_radius);

	vector<Point3f> normal_sum;
	vector<double> normal_weight_sum;
	vector<double> proj_dist;

	normal_sum.assign(dual_samples->vert.size(), Point3f(0., 0., 0.));
	normal_weight_sum.assign(dual_samples->vert.size(), 0);
	//proj_dist.assign(dual_samples->vert.size(), 0);


	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];
		CVertex& v = samples->vert[i];

   		if (dual_v.is_fixed_sample)
   		{
   			continue;
   		}

		for (int j = 0; j < v.neighbors.size(); j++)
		{
			int neighbor_idx = v.neighbors[j];
			if (neighbor_idx < 0 || neighbor_idx >= dual_samples->vert.size())
			{
				cout << "bad index " << neighbor_idx << endl;
				continue;
			}
			CVertex& dual_t = dual_samples->vert[neighbor_idx];
			CVertex& t = samples->vert[neighbor_idx];

			//Point3f diff = dual_v.P() - dual_t.P();
			Point3f diff = v.P() - t.P();

			//double project_dist = dual_t.N() *(diff);
			double project_dist = t.N() *(diff);

			double dist2 = diff.SquaredNorm();

			double weight;

// 			Point3f vm(dual_v.N());
// 			Point3f tm(dual_t.N());
			Point3f vm(v.N());
			Point3f tm(t.N());
			Point3f d = vm - tm;
			double psi = exp(-pow(1 - vm*tm, 2) / pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2));
			double theta = exp(dist2*iradius16);
			weight = psi * theta;
			weight = max(weight, 1e-10);

			normal_weight_sum[i] += weight;
			normal_sum[i] += tm * weight;
			//proj_dist[i] += project_dist * weight;

			//cout << "1" << endl;

		}
	}

	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];

   	if (dual_v.is_fixed_sample)
   	{
   		continue;
   	}

		if (normal_weight_sum[i] > 1e-6)
		{
			//cout << "2" << endl;
			dual_v.N() = normal_sum[i] / normal_weight_sum[i];
			dual_v.N().Normalize();
			//dual_v.P() -= dual_v.N() * (proj_dist[i] / normal_weight_sum[i]);
		}
	}


	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		dual_samples->vert[i].recompute_m_render();
	}
}


void WLOP::runSelfPCA()
{
// 	int knn = para->getDouble("Sefl KNN");
// 	GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, knn, false, "WlopParaDlg::runNormalSmoothing()");
	double radius = para->getDouble("Local Neighbor Size For Surface Points");
	computeConstNeighborhoodUsingRadius(radius);


	vector<Point3f> normals_save(dual_samples->vert.size());
	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		normals_save[i] = dual_samples->vert[i].N();
		dual_samples->vert[i].neighbors = samples->vert[i].neighbors;

		if (use_eigen_neighborhood)
		{
			dual_samples->vert[i].is_fixed_sample = true;
		}
		
	}

	
	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];
		CVertex& v = samples->vert[i];

		std::vector<Point3f> ptVec;

		for (int j = 0; j < v.neighbors.size(); j++)
		{
			int neighbor_idx = v.neighbors[j];
			if (neighbor_idx < 0 || neighbor_idx >= dual_samples->vert.size())
			{
				cout << "bad index " << neighbor_idx << endl;
				continue;
			}
			CVertex& dual_t = dual_samples->vert[neighbor_idx];

			ptVec.push_back(dual_t.P());
		}

		Plane3f plane;
		FitPlaneToPointSet(ptVec, plane);
		dual_v.N() = plane.Direction();
	}

	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		if (normals_save[i] * dual_samples->vert[i].N() < 0)
		{
			dual_samples->vert[i].N() *= -1;
		}
	}

	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		dual_samples->vert[i].recompute_m_render();
	}

}


void WLOP::runSelfProjection()
{
  cout << "sefl projection freeze" << endl;

  double static_radius = para->getDouble("Local Neighbor Size For Surface Points");
  GlobalFun::computeBallNeighbors(samples, NULL, static_radius, samples->bbox);

  for (int k = 0; k < 5; k++)
  {
    int unsettle_number = 0;

    for (int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];
      CVertex& dual_v = dual_samples->vert[i];

      if (!dual_v.is_fixed_sample)
      {
        unsettle_number++;
      }

      v.skel_radius = GlobalFun::computeEulerDist(v.P(), dual_v.P());
    }

    if (unsettle_number < 1)
    {
      break;
    }

    cout << "unsettle_number: " << unsettle_number << endl;

    for (int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];
      CVertex& dual_v = dual_samples->vert[i];

      if (dual_v.is_fixed_sample)
      {
        continue;
      }

      int fixed_neighbor_number = 0;
      double sum_length = 0;

      for (int j = 0; j < v.neighbors.size(); j++)
      {
        CVertex& t = samples->vert[v.neighbors[j]];
        CVertex& dual_t = dual_samples->vert[v.neighbors[j]];
        if (!dual_t.is_fixed_sample)
        {
          continue;
        }

        sum_length += t.skel_radius;
        fixed_neighbor_number++;
      }

      if (fixed_neighbor_number < 1)
      {
        continue;
      }

      double length = sum_length / fixed_neighbor_number;

      dual_v.P() = v.P() - v.N() * length;
      dual_v.is_fixed_sample = true;
    }


  }

}

// void WLOP::runSelfProjection()
// {
//   cout << "sefl projection" << endl;
// 
//   double static_radius = para->getDouble("Local Neighbor Size For Surface Points");
//   GlobalFun::computeBallNeighbors(samples, NULL, static_radius, samples->bbox);
// 
//   for (int k = 0; k < 5; k++)
//   {
//     int unsettle_number = 0;
// 
//     for (int i = 0; i < samples->vert.size(); i++)
//     {
//       CVertex& v = samples->vert[i];
//       CVertex& dual_v = dual_samples->vert[i];
// 
//       if (!dual_v.is_fixed_sample)
//       {
//         unsettle_number++;
//       }
// 
//       v.skel_radius = GlobalFun::computeEulerDist(v.P(), dual_v.P());
//     }
// 
//     if (unsettle_number < 1)
//     {
//       break;
//     }
// 
//     cout << "unsettle_number: " << unsettle_number << endl;
// 
//     for (int i = 0; i < samples->vert.size(); i++)
//     {
//       CVertex& v = samples->vert[i];
//       CVertex& dual_v = dual_samples->vert[i];
// 
//       if (dual_v.is_fixed_sample)
//       {
//         continue;
//       }
// 
//       int fixed_neighbor_number = 0;
//       double sum_length = 0;
// 
//       for (int j = 0; j < v.neighbors.size(); j++)
//       {
//         CVertex& t = samples->vert[v.neighbors[j]];
//         CVertex& dual_t = dual_samples->vert[v.neighbors[j]];
//         if (!dual_t.is_fixed_sample)
//         {
//           continue;
//         }
// 
//         sum_length += t.skel_radius;
//         fixed_neighbor_number++;
//       }
// 
//       if (fixed_neighbor_number < 1)
//       {
//         continue;
//       }
// 
//       double length = sum_length / fixed_neighbor_number;
// 
//       dual_v.P() = v.P() - v.N() * length;
//       dual_v.is_fixed_sample = true;
//     }
// 
// 
//   }
// 
// }


// void WLOP::runSelfProjection()
// {
// 	cout << "Normal smoothing" << endl;
// 
// 	double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
// 	double radius = para->getDouble("CGrid Radius");
// 
// 	double radius2 = radius * radius;
// 	double iradius16 = -4 / radius2;
// 	//int knn = para->getDouble("Sefl KNN");
// 
// 	//GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
// 	//GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
// 	//GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, knn, false, "WlopParaDlg::runNormalSmoothing()");
// 	double local_radius = para->getDouble("Local Neighbor Size For Surface Points");
// 	computeConstNeighborhoodUsingRadius(local_radius);
// 
// 	vector<Point3f> normal_sum;
// 	vector<double> normal_weight_sum;
// 	vector<double> proj_dist;
// 
// 	normal_sum.assign(dual_samples->vert.size(), Point3f(0., 0., 0.));
// 	normal_weight_sum.assign(dual_samples->vert.size(), 0);
// 	proj_dist.assign(dual_samples->vert.size(), 0);
// 
// 
// 	for (int i = 0; i < dual_samples->vert.size(); i++)
// 	{
// 		CVertex& dual_v = dual_samples->vert[i];
// 		CVertex& v = samples->vert[i];
// 
// 		if (dual_v.is_fixed_sample)
// 		{
// 			continue;
// 		}
// 
// 		for (int j = 0; j < v.neighbors.size(); j++)
// 		{
// 			int neighbor_idx = v.neighbors[j];
// 			if (neighbor_idx < 0 || neighbor_idx >= dual_samples->vert.size())
// 			{
// 				cout << "bad index " << neighbor_idx << endl;
// 				continue;
// 			}
// 			CVertex& dual_t = dual_samples->vert[neighbor_idx];
// 			CVertex& t = samples->vert[neighbor_idx];
// 
//  			Point3f diff = dual_v.P() - dual_t.P();
//  			double project_dist = dual_t.N() *(diff);
// // 			Point3f diff = v.P() - t.P();
// // 			double project_dist = t.N() *(diff);
// 			double dist2 = diff.SquaredNorm();
// 
// 			double weight;
// 
//  			Point3f vm(dual_v.N());
//  			Point3f tm(dual_t.N());
// // 			Point3f vm(v.N());
// // 			Point3f tm(t.N());
// 			Point3f d = vm - tm;
// 			double psi = exp(-pow(1 - vm*tm, 2) / pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2));
// 			double theta = exp(dist2*iradius16);
// 			weight = psi * theta;
// 			weight = max(weight, 1e-10);
// 
// 			normal_weight_sum[i] += weight;
// 			normal_sum[i] += tm * weight;
// 			proj_dist[i] += project_dist * weight;
// 		}
// 	}
// 
// 	for (int i = 0; i < dual_samples->vert.size(); i++)
// 	{
// 		CVertex& dual_v = dual_samples->vert[i];
// 
// 		if (dual_v.is_fixed_sample)
// 		{
// 			continue;
// 		}
// 
// 		if (normal_weight_sum[i] > 1e-6)
// 		{
// 			//cout << "2" << endl;
// 			dual_v.N() = normal_sum[i] / normal_weight_sum[i];
// 			dual_v.N().Normalize();
// 			dual_v.P() -= dual_v.N() * (proj_dist[i] / normal_weight_sum[i]);
// 		}
// 	}
// 
// 
// 	for (int i = 0; i < dual_samples->vert.size(); i++)
// 	{
// 		dual_samples->vert[i].recompute_m_render();
// 	}
// }


void WLOP::compute_neighbor_weights(vector<CVertex>& samples,
	                                  vector<CVertex>& target,
	                                  vector< vector<double>>& neighbor_weights,
	                                  double radius,
	                                  double sigma,
																		WeightType type)
{

	vector< vector<int>> neighbors_indexes;
	for (int i = 0; i < samples.size(); i++)
	{
		CVertex& v = samples[i];

		vector<int> neighbor_indexes = v.neighbors;
		neighbors_indexes.push_back(neighbor_indexes);
	}

	compute_neighbor_weights(samples, target, neighbors_indexes, neighbor_weights, radius, sigma);
}

void WLOP::compute_neighbor_weights(vector<CVertex>& samples,
	                                  vector<CVertex>& target,
	                                  vector< vector<int>>& neighbors_indexes,
																		vector< vector<double>>& neighbors_weights,
	                                  double radius,
	                                  double sigma,
																		WeightType type)
{

	double radius2 = radius * radius;
	double iradius16 = -4 / radius2;

	neighbors_weights.clear();

	for (int i = 0; i < samples.size(); i++)
	{
		CVertex& v = samples[i];

		vector<double> neighbor_weights(neighbors_indexes[i].size());
		for (int j = 0; j < neighbors_indexes[i].size(); j++)
		{
			int neighbor_index = neighbors_indexes[i][j];
			CVertex& t = target[j];

			double dist2 = GlobalFun::computeEulerDistSquare(v.P(), t.P());
			double dist_weight = exp(dist2*iradius16);
			double normal_weight = exp(-pow(1 - v.N()*t.N(), 2) / pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2));

			if (type == Bilateral)
			{
				neighbor_weights[j] = dist_weight * normal_weight;
			}
			else if (type == DistanceDiff)
			{
				neighbor_weights[j] = dist_weight;
			}
			else if (type == NormalDiff)
			{
				neighbor_weights[j] = dist_weight;
			}
			else
			{
				neighbor_weights[j] = 1.0;
			}
		}
		neighbors_weights.push_back(neighbor_weights);
	}
}

void WLOP::runComputeEigenNeighborhood(CMesh* dual_samples, CMesh* samples)
{
	runComputeEigenDirections(dual_samples, samples);


	bool use_separate_neighborhood = para->getBool("Use Separate Neighborhood");

	cout << "WLOP:: runComputeEigenNeighborhood" << endl;
	//vector<Point3f> axis(3);
	vector<double> lengthes(3);

	vector< vector<double>> ranges;
	double max_length = 0;

	double eigin_para2 = global_paraMgr.wLop.getDouble("Eigen Neighborhood Para2");
	double eigin_para1 = global_paraMgr.wLop.getDouble("Eigen Neighborhood Para1");
	double radius = para->getDouble("CGrid Radius");

	for (int i = 0; i < dual_samples->vert.size(); i++)
	{
		CVertex& dual_v = dual_samples->vert[i];

		Point3f x_diff = dual_v.eigen_vector0 * radius*eigin_para1 + dual_v.eigen_vector0 * dual_v.eigen_value0 * radius*eigin_para2;
		Point3f y_diff = dual_v.eigen_vector1 * radius*eigin_para1 + dual_v.eigen_vector1 * dual_v.eigen_value1 * radius*eigin_para2;
		Point3f z_diff = dual_v.eigen_vector2 * radius*eigin_para1 + dual_v.eigen_vector2 * dual_v.eigen_value2 * radius*eigin_para2;
		vector<double> lengthes(3);

		lengthes[0] = x_diff.Norm();
		lengthes[1] = y_diff.Norm();
		lengthes[2] = z_diff.Norm();

		if (lengthes[0] > max_length)
		{
			max_length = lengthes[0];
		}
		ranges.push_back(lengthes);
	}

	if (use_separate_neighborhood)
	{
		GlobalFun::computeBallNeighbors(dual_samples, samples, radius, dual_samples->bbox);
	}
	else
	{
		GlobalFun::computeBallNeighbors(dual_samples, samples, max_length, dual_samples->bbox);
		
		for (int i = 0; i < dual_samples->vert.size(); i++)
		{
			CVertex& dual_v = dual_samples->vert[i];
			vector<int> filtered_neighbor;

			vector<double> range = ranges[i];
			for (int j = 0; j < dual_v.original_neighbors.size(); j++)
			{
				int index = dual_v.original_neighbors[j];

				CVertex& t = samples->vert[index];

				Point3f diff = t.P() - dual_v.P();

				double proj_x = abs(diff*dual_v.eigen_vector0);
				double proj_y = abs(diff*dual_v.eigen_vector1);
				double proj_z = abs(diff*dual_v.eigen_vector2);

				if (proj_x < range[0] && proj_y < range[1] && proj_z < range[2])
				{
					filtered_neighbor.push_back(index);
				}

			}

			dual_v.original_neighbors = filtered_neighbor;
		}
	}

	//filter out
	if (1) //!para->getBool("WLOP test bool")
	{
		GlobalFun::computeBallNeighbors(dual_samples, NULL, max_length, dual_samples->bbox);
		for (int i = 0; i < dual_samples->vert.size(); i++)
		{
			CVertex& dual_v = dual_samples->vert[i];
			vector<int> filtered_neighbor;

			vector<double> range = ranges[i];
			for (int j = 0; j < dual_v.neighbors.size(); j++)
			{
				int index = dual_v.neighbors[j];

				CVertex& t = dual_samples->vert[index];

				Point3f diff = t.P() - dual_v.P();

				double proj_x = abs(diff*dual_v.eigen_vector0);
				double proj_y = abs(diff*dual_v.eigen_vector1);
				double proj_z = abs(diff*dual_v.eigen_vector2);

				if (proj_x < range[0] && proj_y < range[1] && proj_z < range[2])
				{
					filtered_neighbor.push_back(index);
				}

			}
			
			dual_v.neighbors = filtered_neighbor;
		}
	}
}

void WLOP::runComputeEigenDirections(CMesh* dual_samples, CMesh* samples)
{
 	GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
 	GlobalFun::computeEigenWithTheta(dual_samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));
	GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);

	//if (!para->getBool("WLOP test bool"))
	if (1)//2015 no time for this
	{
		for (int i = 0; i < dual_samples->vert.size(); i++)
		{
			CVertex& dual_v = dual_samples->vert[i];
			CVertex& v = samples->vert[i];

			vector<Point3f> dirs;
			Point3f dir = (v.P() - dual_v.P()).Normalize();
			dirs.push_back(dir);

			for (int j = 0; j < dual_v.neighbors.size(); j++)
			{
				int index = dual_v.neighbors[j];

				CVertex& dual_t = dual_samples->vert[index];
				CVertex& t = samples->vert[index];

// 				if (j < 3)
// 				{
// 					cout << "v t  " << t.P()[0] << "	" << dual_t.P()[0] << endl;
// 				}

				dir = (t.P() - dual_t.P()).Normalize();
				dirs.push_back(dir);
			}

			double sum_proj_x = 0.0;
			double sum_proj_y = 0.0;
			double sum_proj_z = 0.0;

			for (int j = 0; j < dirs.size(); j++)
			{
				Point3f dir = dirs[j];
				sum_proj_x += abs(dir * dual_v.eigen_vector0);
				sum_proj_y += abs(dir * dual_v.eigen_vector1);
				sum_proj_z += abs(dir * dual_v.eigen_vector2);

// 				if (j < 3)
// 				{
// 					cout << "innter  " << dir[0] << "	" << dual_v.eigen_vector0[0] << endl;
// 				}
			}

			sum_proj_x /= dirs.size();
			sum_proj_y /= dirs.size();
			sum_proj_z /= dirs.size();


// 			if (i < 5)
// 			{
// 				cout << "dir size: " << dirs.size() << endl;
// 				cout << "old eigen: " << dual_v.eigen_value0 << "  " << dual_v.eigen_value1 << "  " << dual_v.eigen_value2 << "  " << endl;
// 				cout << "sum_proj_: " << sum_proj_z << "  " << sum_proj_z << "  " << sum_proj_z << "  " << endl;;
// 			}

			double sum_sum_proj = sum_proj_x + sum_proj_y + sum_proj_z;
			dual_v.eigen_value0 = (sum_sum_proj - sum_proj_x) / sum_sum_proj;
			dual_v.eigen_value1 = (sum_sum_proj - sum_proj_y) / sum_sum_proj;
			dual_v.eigen_value2 = (sum_sum_proj - sum_proj_z) / sum_sum_proj;

			double sum_eigen = dual_v.eigen_value0 + dual_v.eigen_value1 + dual_v.eigen_value2;
			dual_v.eigen_value0 /= sum_eigen;
			dual_v.eigen_value1 /= sum_eigen;
			dual_v.eigen_value2 /= sum_eigen;

			dual_v.eigen_value0 = pow(dual_v.eigen_value0, 8);
			dual_v.eigen_value1 = pow(dual_v.eigen_value1, 8);
			dual_v.eigen_value2 = pow(dual_v.eigen_value2, 8);

			sum_eigen = dual_v.eigen_value0 + dual_v.eigen_value1 + dual_v.eigen_value2;
			dual_v.eigen_value0 /= sum_eigen;
			dual_v.eigen_value1 /= sum_eigen;
			dual_v.eigen_value2 /= sum_eigen;

			if (i < 5)
			{
// 				cout << "new sum_proj: " << sum_proj_x << "  "
// 					    << sum_proj_y<< "  "
// 							<< sum_proj_z << "  sum: "
// 							<< sum_sum_proj << "  " << endl;;
				//cout << "new eigen: " << dual_v.eigen_value0 << "  " << dual_v.eigen_value1 << "  " << dual_v.eigen_value2 << "  " << endl;;
			}
		}
	}
}

void WLOP::computeInitialNeighborSize()
{
	GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, 1, false, "WlopParaDlg::computeInitialNeighborSize()");

	double sum_dist = 0;
	double sum_n = 0;
	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		if (v.neighbors.empty())
		{
			continue;
		}
		sum_n++;
		CVertex& t = samples->vert[v.neighbors[0]];

		double dist = GlobalFun::computeEulerDist(v.P(), t.P());
		sum_dist += dist;
	}
	double average_dist = sum_dist / sum_n;

	global_paraMgr.setGlobalParameter("CGrid Radius", DoubleValue(average_dist * 3.0));
	global_paraMgr.upsampling.setValue("Dist Threshold", DoubleValue(average_dist * average_dist));
	
	global_paraMgr.wLop.setValue("Local Neighbor Size For Inner Points", DoubleValue(average_dist * 3.0));
	global_paraMgr.wLop.setValue("Local Neighbor Size For Surface Points", DoubleValue(average_dist * 5.0));
	global_paraMgr.wLop.setValue("Increasing Step Size", DoubleValue(average_dist * 0.5));
	global_paraMgr.wLop.setValue("Average Closest Dist", DoubleValue(average_dist));


	GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
}


void WLOP::computeDualIndex(CMesh* samples, CMesh* dual_samples)
{

  if (!dual_samples || dual_samples->vert.empty())
  {
    return;
  }

  bool use_progressive_search = para->getBool("Use Progressive Search Index");
  double search_dual_index_para = para->getDouble("Search Dual Index Para");

  if (!global_paraMgr.glarea.getBool("Show Cloest Dual Connection"))
  {
    cout << "use fixed index" << endl;

    for (int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];
      CVertex& dual_v = dual_samples->vert[i];

      v.dual_index = i;
    }

    return;
  }


  cout << "computeDualIndex" << endl;
  Timer timer;
  timer.start("computeDualIndex");

  GlobalFun::computeBallNeighbors(dual_samples, NULL, search_dual_index_para * para->getDouble("CGrid Radius"), dual_samples->bbox);
  bool use_cloest = global_paraMgr.glarea.getBool("Show Cloest Dual Connection");

  for (int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];
    //CVertex& dual_v = dual_samples->vert[i];
    CVertex& dual_v = dual_samples->vert[v.dual_index];

    double min_dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_v.P());
    //double min_dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_v.P(), v.N());



    // 		int dual_idx = i;
    // 		int current_idx = i;
    int dual_idx = v.dual_index;
    int current_idx = v.dual_index;


    for (int j = 0; j < dual_v.neighbors.size(); j++)
    {
      int index = dual_v.neighbors[j];
      CVertex& dual_t = dual_samples->vert[index];

      double dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_t.P());
      //double dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_v.P(), v.N());
      if (dist2 < min_dist2)
      {
        min_dist2 = dist2;
        dual_idx = index;
      }
    }
    


    int max_iterate = 0;
    while ((use_progressive_search && dual_idx != current_idx) || max_iterate++ < 15)
    {
      current_idx = dual_idx;

      CVertex dual_v2 = dual_samples->vert[current_idx];

      if (dual_v2.neighbors.empty())
      {
        break;
      }

      double min_dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_v2.P());
      //double min_dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_v.P(), v.N());

      for (int j = 0; j < dual_v2.neighbors.size(); j++)
      {
        int index = dual_v2.neighbors[j];
        CVertex& dual_t2 = dual_samples->vert[index];

        double dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_t2.P());
        //double dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_v.P(), v.N());

        if (dist2 < min_dist2)
        {
          min_dist2 = dist2;
          dual_idx = index;
        }
      }

    }

    v.dual_index = dual_idx;
  }

  timer.end();

  cout << "$$$$$$$$$$$$$  done new index $$$$$$$$$$$" << endl;

}

// void WLOP::computeDualIndex(CMesh* samples, CMesh* dual_samples)
// {
// 	if (!dual_samples || dual_samples->vert.empty())
// 	{
// 		return;
// 	}
// 
// 	bool use_progressive_search = para->getBool("Use Progressive Search Index");
// 	double search_dual_index_para = para->getDouble("Search Dual Index Para");
// 
// 	if (!global_paraMgr.glarea.getBool("Show Cloest Dual Connection"))
// 	{
// 		cout << "use fixed index" << endl;
// 
// 		for (int i = 0; i < samples->vert.size(); i++)
// 		{
// 			CVertex& v = samples->vert[i];
// 			CVertex& dual_v = dual_samples->vert[i];
// 
// 			v.dual_index = i;
// 		}
// 
// 		return;
// 	}
// 
// 
// 	cout << "computeDualIndex" << endl;
// 	Timer timer;
// 	timer.start("computeDualIndex");
// 
// 	GlobalFun::computeBallNeighbors(dual_samples, NULL, search_dual_index_para * para->getDouble("CGrid Radius"), dual_samples->bbox);
// 	bool use_cloest = global_paraMgr.glarea.getBool("Show Cloest Dual Connection");
// 
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		CVertex& v = samples->vert[i];
// 		//CVertex& dual_v = dual_samples->vert[i];
//     CVertex& dual_v = dual_samples->vert[v.dual_index];
// 
// 		double min_dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_v.P());
// // 		int dual_idx = i;
// // 		int current_idx = i;
//     int dual_idx = v.dual_index;
//     int current_idx = v.dual_index;
// 
// 		if (use_cloest)
// 		{
// 			for (int j = 0; j < dual_v.neighbors.size(); j++)
// 			{
// 				int index = dual_v.neighbors[j];
// 				CVertex& dual_t = dual_samples->vert[index];
// 
// 				double dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_t.P());
// 				if (dist2 < min_dist2)
// 				{
// 					min_dist2 = dist2;
// 					dual_idx = index;
// 				}
// 			}
// 		}
// 
// 
// 		int max_iterate = 0;
// 		while ( (use_progressive_search && dual_idx != current_idx) || max_iterate++ < 15 )
// 		{
// 			current_idx = dual_idx;
// 
// 			CVertex dual_v2 = dual_samples->vert[current_idx];
// 
// 			if (dual_v2.neighbors.empty())
// 			{
// 				break;
// 			}
// 
// 			double min_dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_v2.P());
// 			for (int j = 0; j < dual_v2.neighbors.size(); j++)
// 			{
// 				int index = dual_v2.neighbors[j];
// 				CVertex& dual_t2 = dual_samples->vert[index];
// 
// 				double dist2 = GlobalFun::computeEulerDistSquare(v.P(), dual_t2.P());
// 				if (dist2 < min_dist2)
// 				{
// 					min_dist2 = dist2;
// 					dual_idx = index;
// 				}
// 			}
// 
// 		}
// 
// 		v.dual_index = dual_idx;
// 	}
// 
// 	timer.end();
// }


void WLOP::runTangentialMotion()
{
	initVertexes(true);

	GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
	GlobalFun::computeEigenWithTheta(samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));

	Timer time;
	time.start("Compute Repulsion Term");
	computeRepulsionTerm(samples);
	time.end();

	//double mu = 0.4 * 1.0;
	double mu = para->getDouble("Repulsion Mu");

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		CVertex dual_v = dual_samples->vert[i];

		Point3f diff = v.P() - dual_v.P();

		double proj_dist = diff * dual_v.N();
		v.P() = dual_v.P() + dual_v.N() * proj_dist;

 		if (repulsion_weight_sum[i] > 1e-10)
 		{
 			v.P() += repulsion[i] * (mu / repulsion_weight_sum[i]);
 		}

 		diff = v.P() - dual_v.P();
 		v.N() = ((diff.Normalize() + v.N()) / 2.0).Normalize();
		v.recompute_m_render();
	}
}



void WLOP::run4PCS()
{
  cout << "run4PCS" << endl;

  vector<Point3D> set1, set2;
  Match4PCSOptions options;

  vector<cv::Point2f> tex_coords1, tex_coords2;
  vector<cv::Point3f> normals1, normals2;
  vector<tripple> tris1, tris2;
  vector<std::string> mtls1, mtls2;

  IOManager iomananger;

  // First input.
  //std::string input1 = "input0.obj";
  std::string input1 = "box1.obj";

  // Second input.
  //std::string input2 = "input1.obj";
  std::string input2 = "box2.obj";

  // Read the inputs.
  if (!iomananger.ReadObject(input1.c_str(), set1, tex_coords1, normals1, tris1,
    mtls1)) {
    perror("Can't read input set1");
    exit(-1);
  }

  if (!iomananger.ReadObject(input2.c_str(), set2, tex_coords2, normals2, tris2,
    mtls2)) {
    perror("Can't read input set2");
    exit(-1);
  }



  // Output. The transformed second input.
  std::string output = "output.obj";

  // Delta (see the paper).
  double delta = 5.0;

  // Estimated overlap (see the paper).
  double overlap = 0.2;

  // Threshold of the computed overlap for termination. 1.0 means don't terminate
  // before the end.
  double thr = 1.0;

  // Maximum norm of RGB values between corresponded points. 1e9 means don't use.
  double max_color = 150;

  // Number of sampled points in both files. The 4PCS allows a very aggressive
  // sampling.
  int n_points = 200;

  // Maximum angle (degrees) between corresponded normals.
  double norm_diff = 90.0;

  // Maximum allowed computation time.
  int max_time_seconds = 10;

  bool use_super4pcs = false;

  // Set parameters.
  cv::Mat mat = cv::Mat::eye(4, 4, CV_64F);
  options.overlap_estimation = overlap;
  options.sample_size = n_points;
  options.max_normal_difference = norm_diff;
  options.max_color_distance = max_color;
  options.max_time_seconds = max_time_seconds;
  options.delta = delta;

  // Match and return the score (estimated overlap or the LCP).  
  float score = 0;

  if (use_super4pcs){
//     MatchSuper4PCS matcher(options);
//     cout << "Use Super4PCS" << endl;
//     score = matcher.ComputeTransformation(set1, &set2, &mat);
  }
  else{
    Match4PCS matcher(options);
    cout << "Use old 4PCS" << endl;
    score = matcher.ComputeTransformation(set1, &set2, &mat);
  }
  cout << "Score: " << score << endl;
  cerr << score << endl;

  printf(
    "\n\n%25.3f %25.3f %25.3f %25.3f\n%25.3f %25.3f %25.3f %25.3f\n%25.3f "
    "%25.3f %25.3f %25.3f\n%25.3f %25.3f %25.3f %25.3f\n\n",
    mat.at<double>(0, 0), mat.at<double>(0, 1), mat.at<double>(0, 2),
    mat.at<double>(0, 3), mat.at<double>(1, 0), mat.at<double>(1, 1),
    mat.at<double>(1, 2), mat.at<double>(1, 3), mat.at<double>(2, 0),
    mat.at<double>(2, 1), mat.at<double>(2, 2), mat.at<double>(2, 3),
    mat.at<double>(3, 0), mat.at<double>(3, 1), mat.at<double>(3, 2),
    mat.at<double>(3, 3));

  iomananger.WriteObject((char *)output.c_str(),
    set2,
    tex_coords2,
    normals2,
    tris2,
    mtls2);

}