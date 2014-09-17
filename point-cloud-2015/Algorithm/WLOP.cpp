#include "WLOP.h"


WLOP::WLOP(RichParameterSet* _para)
{
	cout << "WLP constructed!!" << endl;
	para = _para;
	samples = NULL;
	original = NULL;
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
}

void WLOP::setFirstIterate()
{
	nTimeIterated = 0;
}

void WLOP::setInput(DataMgr* pData)
{
	if(!pData->isSamplesEmpty() && !pData->isOriginalEmpty())
	{
		CMesh* _samples = pData->getCurrentSamples();
		CMesh* _original = pData->getCurrentOriginal();

    if (para->getBool("Run Skel WLOP"))
    {
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

    if (para->getBool("Run Dual WLOP"))
    {
      samples = pData->getCurrentDualSamples();
      original = pData->getCurrentSamples();
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


	vi_end = original->vert.end();
	i = 0;
	for(vi = original->vert.begin(); vi != vi_end; ++vi) 
	{
		vi->m_index = i++;
		box.Add(vi->P());
	}
	original->bbox = box;


	repulsion.assign(samples->vn, vcg::Point3f(0, 0, 0));
	average.assign(samples->vn, vcg::Point3f(0, 0, 0));

	repulsion_weight_sum.assign(samples->vn, 0);
	average_weight_sum.assign(samples->vn, 0);

  samples_average.assign(samples->vn, vcg::Point3f(0, 0, 0));
  samples_average_weight_sum.assign(samples->vn, 0);

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
    stepForward();
    return;
  }

  if (para->getBool("Run Compute Initial Sample Neighbor"))
  {
    computeInitialSampleNeighbor();
    return;
  }

  if (para->getBool("Run Skel WLOP"))
  {
    runSkelWlop();
    return;
  }

  if (para->getBool("Run Dual Drag WLOP"))
  {
    pioneer_points_id.clear();
    pioneer_points_position.clear();
    pioneer_points_origininal.clear();

    computeJointNeighborhood();
    samples = dual_samples;
    if (para->getBool("Need Sample Average"))
    {
      runDragWlop();
    }
    return;
  }

  if (para->getBool("Run Regularize Samples"))
  {
    runRegularizeSamples();
    return;
  }

  if (para->getBool("Run Regularize Normals"))
  {
    runRegularizeNormals();
    return;
  }

	//int nTimes = para->getDouble("Num Of Iterate Time");
	for(int i = 0; i < 1; i++)
	{ 
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

  for(int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];

    for (int j = 0; j < v.original_neighbors.size(); j++)
    {
      CVertex& t = original->vert[v.original_neighbors[j]];

      Point3f diff = v.P() - t.P();
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

      average[i] += t.P() * w;  
      average_weight_sum[i] += w;  
    }

    for (int k = 0; k < v.neighbors.size(); k++)
    {
      CVertex& t = samples->vert[v.neighbors[k]];

      Point3f diff = v.P() - t.P();
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

      average[i] += t.P() * w;  
      average_weight_sum[i] += w;  
    }
  }
}

void WLOP::computeAverageTerm(CMesh* samples, CMesh* original)
{
	double average_power = para->getDouble("Average Power");
	bool need_density = para->getBool("Need Compute Density");
  bool use_elliptical_neighbor = para->getBool("Use Elliptical Original Neighbor");
	double radius = para->getDouble("CGrid Radius"); 

	double radius2 = radius * radius;
	double iradius16 = -para->getDouble("H Gaussian Para")/radius2;
  
  double close_threshold = radius2 / 16;
  bool use_tangent = para->getBool("Use Tangent Vector");


  bool L2 = para->getBool("Need Sample Average");
  
  cout << "dual" << endl;
	cout << "Original Size:" << samples->vert[0].original_neighbors.size() << endl;
  for(int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];

		for (int j = 0; j < v.original_neighbors.size(); j++)
		{
			CVertex& t = original->vert[v.original_neighbors[j]];
			
			Point3f diff = t.P() - v.P();
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

      if (use_elliptical_neighbor)
      {
        w *= exp(proj_dist2 * iradius16);
      }

			if (need_density)
			{
				w *= original_density[t.m_index];
			}

      if (L2)
      {
        w = 1.0;
      }

      if (use_tangent)
      {
//         if (i < 10)
//         {
//           cout << "tangent" << endl;
//         }
        Point3f proj_point = v.P() + v.N() * proj_dist;
        average[i] += proj_point * w;  

      }
      else
      {
        average[i] += t.P() * w;  
      }
			average_weight_sum[i] += w;  


      //if (use_adaptive_mu && !is_sample_close_to_original[v.m_index])
      //{
      //  if (dist2 < close_threshold)
      //  {
      //    is_sample_close_to_original[v.m_index] = true;
      //  }
      //}
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

	cout << endl<< endl<< "Sample Neighbor Size:" << samples->vert[0].neighbors.size() << endl<< endl;
	for(int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& t = samples->vert[v.neighbors[j]];
			Point3f diff = v.P() - t.P();

      if (use_tangent)
      {
        diff = GlobalFun::getTangentVector(diff, v.N());
      }

			double dist2  = diff.SquaredNorm();
			double len = sqrt(dist2);
			if(len <= 0.001 * radius) len = radius*0.001;

			double w = exp(dist2*iradius16);
			double rep = w * pow(1.0 / len, repulsion_power);

			if (need_density)
			{
				rep *= samples_density[t.m_index];
			}

			repulsion[i] += diff * rep;  
			repulsion_weight_sum[i] += rep;
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
//       double len = sqrt(dist2);
//       if(len <= 0.001 * radius) len = radius*0.001;

      double w = exp(dist2*iradius16);
      //double rep = w * pow(1.0 / len, repulsion_power);

      if (need_density)
      {
        w *= samples_density[t.m_index];
      }

      //w = 1.0;

      samples_average[i] += t.P() * w;  

//       if (i < 5)
//       {
//         cout << t.P().X() << "  " << t.P().Y() << " " << t.P().Z() << endl;
//       }
      samples_average_weight_sum[i] += w;
    }
  }
}

void WLOP::computeDensity(bool isOriginal, double radius)
{
	CMesh* mesh;
	if (isOriginal)
	{
		mesh = original;
	}
	else
	{
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

void WLOP::updateAdaptiveNeighbor()
{
  double radius = para->getDouble("CGrid Radius"); 
  double radius2 = radius * radius;

  for (int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];

    vector<int> new_neighbors;
    for (int j = 0; j < v.neighbors.size(); j++)
    {
      CVertex& t = samples->vert[v.neighbors[j]];

      float dist2 = GlobalFun::computeEulerDistSquare(v.P(), t.P());
      //Point3f diff = t.P() - v.P();
      //diff.normalized();

      if (dist2 < radius2 && (t.N() * v.N()) > 0)
      {
        new_neighbors.push_back(v.neighbors[j]);
      }
    }

    v.neighbors.clear();
    for (int j = 0; j < new_neighbors.size(); j++)
    {
      v.neighbors.push_back(new_neighbors[j]);
    }
  }
}

double WLOP::iterate()
{
  use_adaptive_mu = para->getBool("Use Adaptive Mu");
  is_sample_close_to_original.assign(samples->vert.size(), false);
  bool use_tangent = para->getBool("Use Tangent Vector");

	Timer time;

	initVertexes(true);

	time.start("Sample Original Neighbor Tree!!!");
	GlobalFun::computeBallNeighbors(samples, original, 
		para->getDouble("CGrid Radius"), box);
	time.end();

  if (!para->getBool("Use Adaptive Sample Neighbor"))
  {
    time.start("Sample Sample Neighbor Tree");
    GlobalFun::computeBallNeighbors(samples, NULL, 
      para->getDouble("CGrid Radius"), samples->bbox);
    time.end();
  }
  else
  {
    updateAdaptiveNeighbor();
  }

	
	if (nTimeIterated == 0) 
	{
		if (para->getBool("Need Compute Density"))
		{
			double local_density_para = 0.95;
			time.start("Original Original Neighbor Tree");
			GlobalFun::computeBallNeighbors(original, NULL, 
				para->getDouble("CGrid Radius") * local_density_para, original->bbox);
			time.end();

			time.start("Compute Original Density");
			original_density.assign(original->vn, 0);

			computeDensity(true, para->getDouble("CGrid Radius") * local_density_para);
			time.end();
		}
	}

	if (para->getBool("Need Compute Density"))
	{
		time.start("Compute Density For Sample");
		computeDensity(false, para->getDouble("CGrid Radius"));
		time.end();
	}

	time.start("Sample Original Neighbor Tree!!!");
	GlobalFun::computeBallNeighbors(samples, original, 
		para->getDouble("CGrid Radius"), box);
	time.end();

  time.start("Compute Average Term");
  if (para->getBool("Original Combine Sample"))
  {
    computeAverageAddSampleTerm(samples, original);
  }
  else
  {
    computeAverageTerm(samples, original);
  }
	time.end();

	time.start("Compute Repulsion Term");
	computeRepulsionTerm(samples);
	time.end();

  bool need_sample_average_term = false;
  if (para->getBool("Need Sample Average"))
  {
    time.start("Compute Repulsion Term");
    computeSampleAverageTerm(samples);
    time.end();
    need_sample_average_term = true;
  }

	double mu = para->getDouble("Repulsion Mu");
  double mu3 = para->getDouble("Sample Average Mu3");

  double current_mu = para->getDouble("Repulsion Mu");

	Point3f c;

  vector<Point3f> new_sample_positions;
  vector<float> move_proj_vec;

  double radius = para->getDouble("CGrid Radius");
  double radius2 = radius * radius;
  double iradius16 = -para->getDouble("H Gaussian Para")/radius2;

  if (use_adaptive_mu)
  {
    for (int i = 0; i < samples->vert.size(); i++)
    {
      if (is_sample_close_to_original[i])
      {
        samples->vert[i].is_fixed_sample = true;
      }
      else
      {
        samples->vert[i].is_fixed_sample = false;
      }
    }
  }

  if (para->getBool("Need Averaging Movement"))
  {
    new_sample_positions.assign(samples->vert.size(), Point3f(0.,0.,0.));
    move_proj_vec.assign(samples->vert.size(), 0.);

    for (int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];
      Point3f temp_p = v.P();

      if (average_weight_sum[i] > 1e-20)
      {
        temp_p = average[i] / average_weight_sum[i];
      }

      if (use_adaptive_mu)
      {
        if (is_sample_close_to_original[i])
        {
          current_mu = mu;
        }
        else
        {
          current_mu = mu3;
        }
      }

      if (repulsion_weight_sum[i] > 1e-20 && current_mu > 0)
      {
        temp_p += repulsion[i] * (current_mu / repulsion_weight_sum[i]);
      }

      Point3f move_vector = temp_p - v.P();
      float move_proj = move_vector * v.N();
      move_proj_vec[i] = move_proj;
    }

    for (int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];

      float sum_move_proj = 0.0;
      float sum_w = 0.0;

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
        v.P() += v.N() * avg_move_proj;
      }

    }
  }
  else
  {
    for(int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];
      c = v.P();

      if (average_weight_sum[i] > 1e-20)
      {
        v.P() = average[i] / average_weight_sum[i];
      }

      if (repulsion_weight_sum[i] > 1e-20 && mu > 0)
      {
        v.P() +=  repulsion[i] * (mu / repulsion_weight_sum[i]);
      }

      if (need_sample_average_term && samples_average_weight_sum[i] > 1e-20 && mu3 >= 0)
      {
        v.P() +=  samples_average[i] * (mu3 / samples_average_weight_sum[i]);
      }

      if (/*average_weight_sum[i] > 1e-20 && */repulsion_weight_sum[i] > 1e-20 )
      {
        Point3f diff = v.P() - c; 
        double move_error = sqrt(diff.SquaredNorm());
        error_x += move_error; 
      }
    }
    error_x = error_x / samples->vn;
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


void WLOP::recomputePCA_Normal()
{
  CMesh temp_mesh;
  //if (para->getBool("Use Adaptive Sample Neighbor"))
  if (false)
  {

    //double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
    //double radius = para->getDouble("CGrid Radius"); 

    //double radius2 = radius * radius;
    //double iradius16 = -4 / radius2;

    ////CMesh* samples = mesh;
    ////GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);

    //vector<Point3f> normal_sum;
    //vector<float> normal_weight_sum;

    //normal_sum.assign(samples->vert.size(), Point3f(0.,0.,0.));
    //normal_weight_sum.assign(samples->vert.size(), 0);

    //for(int i = 0; i < samples->vert.size(); i++)
    //{
    //  CVertex& v = samples->vert[i];

    //  for (int j = 0; j < v.neighbors.size(); j++)
    //  {
    //    CVertex& t = samples->vert[v.neighbors[j]];

    //    Point3f diff = v.P() - t.P();
    //    double dist2  = diff.SquaredNorm();

    //    double rep; 

    //    Point3f vm(v.N());
    //    Point3f tm(t.N());
    //    Point3f d = vm-tm;
    //    double psi = exp(-pow(1-vm*tm, 2)/pow(max(1e-8,1-cos(sigma/180.0*3.1415926)), 2));
    //    double theta = exp(dist2*iradius16);
    //    rep = psi * theta;
    //    rep = max(rep, 1e-10);

    //    normal_weight_sum[i] += rep;
    //    normal_sum[i] += tm * rep;         
    //  }

    //  if (normal_weight_sum[i] > 1e-6)
    //  {
    //    v.N() = normal_sum[i] / normal_weight_sum[i];
    //  }
    //}





    //for(int i = 0; i < samples->vn; i++)
    //{
    //  mesh_temp[i].P() = samples->vert[i].P();
    //}

    //int knn = global_paraMgr.norSmooth.getInt("PCA KNN");
    //vcg::NormalExtrapolation<vector<CVertex> >::ExtrapolateNormalsWithExistingNeighbor(mesh_temp.begin(), mesh_temp.end());

    //for(int i = 0; i < samples->vn; i++)
    //{
    //  Point3f& new_normal = mesh_temp[i].N();
    //  CVertex& v = samples->vert[i];
    //  if (v.N() * new_normal > 0)
    //  {
    //    v.N() = new_normal;
    //  }
    //  else
    //  {
    //    v.N() = -new_normal;
    //  }
    //  v.recompute_m_render();
    //}




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

void WLOP::computeInitialSampleNeighbor()
{
  
  Timer time;
  time.start("Sample Sample Neighbor Tree");
  GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
  time.end();

}


void WLOP::runSkelWlop()
{
  cout << "runSkelWlop" << endl;

  Timer time;

  initVertexes(true);

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
      computeDensity(true, para->getDouble("CGrid Radius"));
    }
    time.end();
  }

  time.start("Sample Original neighbor");
  GlobalFun::computeBallNeighbors(samples, original, 
    para->getDouble("CGrid Radius"), box);
  time.end();

  time.start("computeAverageTerm");
  computeAverageTerm(samples, original);
  time.end();

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

  double mu_max = para->getDouble("Repulsion Mu");
  double mu_min = para->getDouble("Repulsion Mu2");
  double mu_length = abs(mu_max - mu_min);
  double sigma_length = abs(max_sigma - min_sigma);
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

    double mu = (mu_length / sigma_length) * (v.eigen_confidence - min_sigma) + mu_min;

    if (average_weight_sum[i] > 1e-20)
    {
      v.P() = average[i] / average_weight_sum[i];

    }
    if (repulsion_weight_sum[i] > 1e-20 && mu >= 0)
    {
      v.P() +=  repulsion[i] * (mu / repulsion_weight_sum[i]);
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


void WLOP::runDragWlop()
{
  use_adaptive_mu = para->getBool("Use Adaptive Mu");

  Timer time;

  initVertexes(false);

//   time.start("Sample Original Neighbor Tree!!!");
//   GlobalFun::computeBallNeighbors(samples, original, 
//     para->getDouble("CGrid Radius"), box);
//   time.end();

//   if (!para->getBool("Use Adaptive Sample Neighbor"))
//   {
//     time.start("Sample Sample Neighbor Tree");
//     GlobalFun::computeBallNeighbors(samples, NULL, 
//       para->getDouble("CGrid Radius"), samples->bbox);
//     time.end();
//   }
//   else
//   {
//     updateAdaptiveNeighbor();
//   }


//   if (nTimeIterated == 0) 
//   {
//     if (para->getBool("Need Compute Density"))
//     {
// //       double local_density_para = 0.95;
// //       time.start("Original Original Neighbor Tree");
// //       GlobalFun::computeBallNeighbors(original, NULL, 
// //         para->getDouble("CGrid Radius") * local_density_para, original->bbox);
// //       time.end();
// 
//       time.start("Compute Original Density");
//       original_density.assign(original->vn, 0);
// 
//       computeDensity(true, para->getDouble("CGrid Radius"));
//       time.end();
//     }
// 
//   }

  if (para->getBool("Need Compute Density"))
  {
    time.start("Compute Density For Sample");
    computeDensity(false, para->getDouble("CGrid Radius"));
    time.end();
  }

//   time.start("Sample Original Neighbor Tree!!!");
//   GlobalFun::computeBallNeighbors(samples, original, 
//     para->getDouble("CGrid Radius"), box);
//   time.end();

  time.start("Compute Average Term");
//   if (para->getBool("Original Combine Sample"))
//   {
//     computeAverageAddSampleTerm(samples, original);
//   }
//   else
//   {
    computeAverageTerm(samples, original);
  //}
  time.end();

  time.start("Compute Repulsion Term");
  computeRepulsionTerm(samples);
  time.end();

//   bool need_sample_average_term = false;
//   if (para->getBool("Need Sample Average"))
//   {
//     time.start("Compute Repulsion Term");
//     computeSampleAverageTerm(samples);
//     time.end();
//     need_sample_average_term = true;
//   }

  double mu = para->getDouble("Repulsion Mu");
  double mu3 = para->getDouble("Sample Average Mu3");

  double current_mu = para->getDouble("Repulsion Mu");

  Point3f c;

  vector<Point3f> new_sample_positions;
  vector<float> move_proj_vec;

  double radius = para->getDouble("CGrid Radius");
  double radius2 = radius * radius;
  double iradius16 = -para->getDouble("H Gaussian Para")/radius2;

  if (use_adaptive_mu)
  {
    for (int i = 0; i < samples->vert.size(); i++)
    {
      if (is_sample_close_to_original[i])
      {
        samples->vert[i].is_fixed_sample = true;
      }
      else
      {
        samples->vert[i].is_fixed_sample = false;
      }
    }
  }

  if (para->getBool("Need Averaging Movement"))
  {
    new_sample_positions.assign(samples->vert.size(), Point3f(0.,0.,0.));
    move_proj_vec.assign(samples->vert.size(), 0.);

    for (int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];
      Point3f temp_p = v.P();

      if (average_weight_sum[i] > 1e-20)
      {
        temp_p = average[i] / average_weight_sum[i];
      }

      if (use_adaptive_mu)
      {
        if (is_sample_close_to_original[i])
        {
          current_mu = mu;
        }
        else
        {
          current_mu = mu3;
        }
      }

      if (repulsion_weight_sum[i] > 1e-20 && current_mu > 0)
      {
        temp_p += repulsion[i] * (current_mu / repulsion_weight_sum[i]);
      }

      Point3f move_vector = temp_p - v.P();
      float move_proj = move_vector * v.N();
      move_proj_vec[i] = move_proj;
    }

    for (int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];

      float sum_move_proj = 0.0;
      float sum_w = 0.0;

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
        v.P() += v.N() * avg_move_proj;
      }

    }
  }
  else
  {
    for(int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];
      c = v.P();

      if (average_weight_sum[i] > 1e-20)
      {
        v.P() = average[i] / average_weight_sum[i];
      }

      if (repulsion_weight_sum[i] > 1e-20 && mu > 0)
      {
        v.P() +=  repulsion[i] * (mu / repulsion_weight_sum[i]);
      }

//       if (need_sample_average_term && samples_average_weight_sum[i] > 1e-20 && mu3 >= 0)
//       {
//         v.P() +=  samples_average[i] * (mu3 / samples_average_weight_sum[i]);
//       }

      if (/*average_weight_sum[i] > 1e-20 && */repulsion_weight_sum[i] > 1e-20 )
      {
        Point3f diff = v.P() - c; 
        double move_error = sqrt(diff.SquaredNorm());
        error_x += move_error; 
      }
    }
    error_x = error_x / samples->vn;
  }


  para->setValue("Current Movement Error", DoubleValue(error_x));
  cout << "****finished compute WLOP error:	" << error_x << endl;





  if (para->getBool("Need Compute PCA"))
  {
    for (int i = 0; i < pioneer_points_id.size(); i++)
    {
      CVertex& dual_v = dual_samples->vert[pioneer_points_id[i]];
      Point3f dual_v_last = pioneer_points_position[i];
      Point3f movement = dual_v.P() - dual_v_last;

      if (sqrt(movement.SquaredNorm()) < error_x * 1.5)
      {
        continue;
      }

      for (int j = 0; j < pioneer_points_origininal[i].size(); j++)
      {
        CVertex new_v;
        new_v.bIsOriginal = true;
        new_v.P() = pioneer_points_origininal[i][j] + movement;
        new_v.m_index = original->vert.size();
        original->vert.push_back(new_v);
      }

    }
    original->vn = original->vert.size();


    //time.start("Recompute PCA");
    //recomputePCA_Normal();
    //time.end();
  }
}


void WLOP::runRegularizeSamples()
{
  cout << "runRegularizeSamples" << endl;

  GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
  GlobalFun::computeEigenWithTheta(samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));

//   int branch_KNN = para->getDouble("Branch Search KNN");
//   GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, branch_KNN, false, "void Skeletonization::searchNewBranches()");

  vector<Point3f> new_sample_set;
  new_sample_set.assign(samples->vert.size(), Point3f(0, 0, 0));
  for (int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];
    Point3f direction = v.eigen_vector0.Normalize();

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

//void WLOP::runRegularizeNormals()
//{
//  cout << "runRegularizeSamples" << endl;
//
//  
//
//  GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
//
//  if (samples->vert.size() < 3)
//  {
//    return;
//  }
//  Point3f direction = (samples->vert[0].N() ^ samples->vert[1].N()).Normalize();
//
//  vector<Point3f> new_normals(samples->vert.size());
//  for (int i = 0; i < samples->vert.size(); i++)
//  {
//    CVertex& v = samples->vert[i];
//
//    double min_clockwise_angle = 500;
//    double max_clockwise_angle = -500;
//
//    double max_anticlockwise_angle = -500;
//    double min_anticlockwise_angle = 500;
//
//    Point3f clockwise_normal = v.N();
//    Point3f anticlockwise_normal = v.N();
//
//    Point3f temp_clockwise_normal = v.N();
//    Point3f temp_anticlockwise_normal = v.N();
//
//    for (int j = 0; j < samples->vert.size(); j++)
//    {
//      CVertex& t = samples->vert[j];
//
//      double angle = GlobalFun::computeDirectionalAngleOfTwoVertor(v.N(), t.N(), direction);
//      if ((angle < 0.000001 && angle > -0.000001) /*|| angle > 179.99999 || angle < -179.99999*/)
//      {
//        //t.N() *= -1;
//        //cout << "skip!!" << endl;
//        //break;;
//      }
//      if (angle > 0)
//      {
//        if (angle < min_clockwise_angle)
//        {
//          min_clockwise_angle = angle;
//          clockwise_normal = t.N();
//        }
//        if (angle > max_clockwise_angle)
//        {
//          max_clockwise_angle = angle;
//          temp_clockwise_normal = t.N();
//        }
//      }
//      else
//      {
//        if (angle > max_anticlockwise_angle)
//        {
//          max_anticlockwise_angle = angle;
//          anticlockwise_normal = t.N();
//        }
//        if (angle < min_anticlockwise_angle)
//        {
//          min_anticlockwise_angle = angle;
//          temp_anticlockwise_normal = t.N();
//        }
//      }
//    }
//
//    if (min_clockwise_angle > 360)
//    {
//      cout << "big min_clockwise_angle circle  " << min_clockwise_angle << endl;
//      clockwise_normal = temp_clockwise_normal;
//    }
//    if (max_anticlockwise_angle < -360)
//    {
//      cout << "small max_anticlockwise_angle circle" << max_anticlockwise_angle << endl;
//      anticlockwise_normal = temp_anticlockwise_normal;
//    }
//
//    Point3f new_normal = (clockwise_normal + anticlockwise_normal) / 2.0;
//    if ((min_clockwise_angle > 360) || (max_anticlockwise_angle < -360))
//    {
//      new_normal *= -1;
//    }
//    new_normals[i] = new_normal.Normalize();
//  }
//
//  for (int i = 0 ; i < new_normals.size(); i++)
//  {
//    CVertex& v = samples->vert[i];
//    v.N() = new_normals[i];
//  }
//
//}


//void WLOP::runRegularizeNormals()
//{
//  cout << "runRegularizeSamples" << endl;
//
//  GlobalFun::removeNormalOverlaps(samples);
//
//  GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
//
//  if (samples->vert.size() < 3)
//  {
//    return;
//  }
//  Point3f direction = (samples->vert[0].N() ^ samples->vert[1].N()).Normalize();
//
//  vector<Point3f> new_normals(samples->vert.size());
//  for (int i = 0; i < samples->vert.size(); i++)
//  {
//    CVertex& v = samples->vert[i];
//
//    double min_clockwise_angle = 500;
//    double max_clockwise_angle = -500;
//
//    Point3f min_clockwise_normal = v.N();
//    Point3f max_clockwise_normal = v.N();
//
//    for (int j = 0; j < samples->vert.size(); j++)
//    {
//      CVertex& t = samples->vert[j];
//
//      double angle = GlobalFun::computeDirectionalAngleOfTwoVertor(v.N(), t.N(), direction);
//
//      if (angle < min_clockwise_angle)
//      {
//        min_clockwise_angle = angle;
//        min_clockwise_normal = t.N();
//      }
//      if (angle > max_clockwise_angle)
//      {
//        max_clockwise_angle = angle;
//        max_clockwise_normal = t.N();
//      }
//
////       if ((angle < 0.000001 && angle > -0.000001) /*|| angle > 179.99999 || angle < -179.99999*/)
////       {
////         //t.N() *= -1;
////         //cout << "skip!!" << endl;
////         //break;;
////       }
////       if (angle > 0)
////       {
////         if (angle < min_clockwise_angle)
////         {
////           min_clockwise_angle = angle;
////           min_clockwise_normal = t.N();
////         }
////         if (angle > max_clockwise_angle)
////         {
////           max_clockwise_angle = angle;
////           max_clockwise_normal = t.N();
////         }
////       }
//    }
//
//    Point3f new_normal = (min_clockwise_normal + max_clockwise_normal) / 2.0;
//    if (max_clockwise_angle < 180)
//    {
//      cout << "small circle" << endl;
//      new_normal *= -1;
//    }
//    new_normals[i] = new_normal.Normalize();
//  }
//
//  for (int i = 0 ; i < new_normals.size(); i++)
//  {
//    CVertex& v = samples->vert[i];
//    v.N() = new_normals[i];
//  }
//
//}



void WLOP::runRegularizeNormals()
{
  GlobalFun::computeAnnNeigbhors(samples->vert, dual_samples->vert, 1, false, "WlopParaDlg::runRegularizeNormals()");
  GlobalFun::computeAnnNeigbhors(dual_samples->vert, samples->vert, 1, false, "WlopParaDlg::runRegularizeNormals()");

  for(int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];

	int neighbor_idx = v.neighbors[0];

    CVertex& dual_v = dual_samples->vert[neighbor_idx];

	Point3f dir = (v.P() - dual_v.P()).Normalize();
    v.N() = ((dir+v.N())/2.0).Normalize();
  }
}




void WLOP::computeJointNeighborhood()
{
  GlobalFun::computeBallNeighbors(dual_samples, NULL, para->getDouble("CGrid Radius"), dual_samples->bbox);
  //GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius") * 3, samples->bbox);
  
  int sigma_KNN = 800;
  GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, sigma_KNN, false, "void Skeletonization::eigenThresholdClassification()");

  Timer time;
  time.start("Sample Original Neighbor Tree!!!");
  GlobalFun::computeBallNeighbors(dual_samples, original, 
    para->getDouble("CGrid Radius"), box);
  time.end();

  for (int i = 0; i < dual_samples->vert.size(); i++)
  {
    CVertex& dual_v = dual_samples->vert[i];
    CVertex& v = samples->vert[i];

    v.N() = dual_v.N();
  }

  vector<NeighborDisk> neighbor_disks;
  for (int i = 0; i < dual_samples->vert.size(); i++)
  {
    CVertex& dual_v = dual_samples->vert[i];
    NeighborDisk disk(dual_v.P(), dual_v.N());
    for (int j = 0; j < dual_v.neighbors.size(); j++)
    {
      disk.add_point(dual_samples->vert[dual_v.neighbors[j]]);
    }
    neighbor_disks.push_back(disk);
  }

  double angle = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
  //double angle_threshold = angle * 3.1415926 / 180.;
  vector< vector<int> > new_neighbors_vec;
  vector< set<int> > new_original_neighbors_vec;
  vector<int> update_dual_sample_indexes;

  for (int i = 0; i < dual_samples->vert.size(); i++)
  {
    double occupy_percentage = neighbor_disks[i].getOccupyPercentage();
    if (occupy_percentage > 0.85)
      continue;

    CVertex& dual_v = dual_samples->vert[i];
    CVertex& v = samples->vert[i];

    vector<int> new_neighbors;
    set<int> new_original_neighbors;
    //vector<int> new_original_neighbors;

    for (int j = 0; j < v.neighbors.size(); j++)
    {
      CVertex& t = samples->vert[v.neighbors[j]];
      double normal_diff = GlobalFun::computeRealAngleOfTwoVertor(t.N(), v.N()); 
      
      //new_neighbors.push_back(v.neighbors[j]);
      if (normal_diff > angle)
      {
        continue;
      }

      CVertex& dual_t = dual_samples->vert[v.neighbors[j]];
      //bool is_slot_occupied = neighbor_disks[i].isSlotOccupied(dual_t.P());
      //if (is_slot_occupied)
      //{
      //  continue;
      //}
      bool is_new_point_good = neighbor_disks[i].isNewPointGood(dual_t.P());
      if (!is_new_point_good)
      {
        continue;
      }
      neighbor_disks[i].add_point(dual_t.P());

      new_neighbors.push_back(v.neighbors[j]);
      for (int k = 0; k < dual_t.original_neighbors.size(); k++)
      {
        new_original_neighbors.insert(dual_t.original_neighbors[k]);
        //new_original_neighbors.push_back(dual_t.original_neighbors[k]);
      }   
 
    }

    new_neighbors_vec.push_back(new_neighbors);
    new_original_neighbors_vec.push_back(new_original_neighbors);
    update_dual_sample_indexes.push_back(i);
  }

  for (int i = 0; i < update_dual_sample_indexes.size(); i++)
  {
    pioneer_points_id.push_back(update_dual_sample_indexes[i]);
    CVertex& dual_v = dual_samples->vert[update_dual_sample_indexes[i]];

    pioneer_points_position.push_back(dual_v.P());
    vector<Point3f> pioneer_original;
    for (int j = 0; j < dual_v.original_neighbors.size(); j++)
    {
      pioneer_original.push_back(original->vert[dual_v.original_neighbors[j]]);
    }
    pioneer_points_origininal.push_back(pioneer_original);
  }

  for (int i = 0; i < update_dual_sample_indexes.size(); i++)
  {
    CVertex& dual_v = dual_samples->vert[update_dual_sample_indexes[i]];

    vector<int>& new_neighbors = new_neighbors_vec[i];
    
    if (i < 15)
    {
      cout << "before:  " << dual_v.original_neighbors.size() << "  "<<  dual_v.neighbors.size() << endl;
    }
    set<int>& new_original_neighbors_set = new_original_neighbors_vec[i];
    for (int j = 0; j < dual_v.original_neighbors.size(); j++)
    {
      new_original_neighbors_set.insert(dual_v.original_neighbors[j]);
    }


    set<int> new_neighbors_set;
    for (int j = 0; j < dual_v.neighbors.size(); j++)
    {
      new_neighbors_set.insert(dual_v.neighbors[j]);
    }
    for (int j = 0; j < new_neighbors.size(); j++)
    {
      new_neighbors_set.insert(new_neighbors[j]);
    }

    set<int>::iterator iter;

    dual_v.neighbors.clear();
    for (iter = new_neighbors_set.begin(); iter != new_neighbors_set.end(); ++iter)
    {
      dual_v.neighbors.push_back(*iter);
    }

    dual_v.original_neighbors.clear();
    for (iter = new_original_neighbors_set.begin(); iter != new_original_neighbors_set.end(); ++iter)
    {
      dual_v.original_neighbors.push_back(*iter);
    }
    if (i < 15)
    {
      cout << "after:  " << dual_v.original_neighbors.size() << "  "<<  dual_v.neighbors.size() << endl << endl;
    }
  }
}


void WLOP::runProjection()
{

  //if (global_paraMgr.glarea.getDouble("Picked Index") < 0.1)
  //    return;

  int pick_index = global_paraMgr.glarea.getDouble("Picked Index");

  vector<Point3f> new_new_normals(samples->vert.size());
  //for (pick_index = 0; pick_index < samples->vert.size(); pick_index++)
  {
    GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
    GlobalFun::computeEigenWithTheta(samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));

    CVertex pick_v = samples->vert[pick_index];
    pick_v.N().Normalize();

    vector<CVertex> local_samples;
    local_samples.push_back(pick_v);

    Point3f X_axis = pick_v.N().Normalize();
    Point3f random_dir = pick_v.eigen_vector0.Normalize();
    Point3f Y_axis = (X_axis ^ random_dir).Normalize();
    Point3f Z_axis = (X_axis ^ Y_axis).Normalize();

    for (int j = 0; j < pick_v.neighbors.size(); j++)
    {
      CVertex t = samples->vert[pick_v.neighbors[j]];
      t.P() = pick_v.P();
      float X_proj_dist = t.N() * X_axis;
      float Y_proj_dist = t.N() * Y_axis;
      Point3f new_normal = ((t.P() + X_axis * X_proj_dist + Y_axis * Y_proj_dist) - t.P()).Normalize();
      t.N() = new_normal.Normalize();
      local_samples.push_back(t);
    }

    CMesh local_proj_samples;
    for (int i = 0; i < local_samples.size(); i++)
    {
      CVertex v = local_samples[i];
      local_proj_samples.vert.push_back(v);
    }
    local_proj_samples.vn = local_proj_samples.vert.size();


//        samples->vert.clear();
//        for (int i = 0; i < local_samples.size(); i++)
//        {
//          samples->vert.push_back(local_samples[i]);
//        }
//        samples->vn = samples->vert.size();
// 
//        return;


    GlobalFun::computeBallNeighbors(&local_proj_samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);

    if (local_proj_samples.vert.size() < 3)
    {
      cout << "small neighbor" << endl;
      return;
    }
    Point3f direction = (local_proj_samples.vert[0].N() ^ local_proj_samples.vert[1].N()).Normalize();

    vector<Point3f> new_normals(local_proj_samples.vert.size());
    CVertex test0_min;
    CVertex test1_max;
    test0_min.P() = pick_v.P();
    test1_max.P() = pick_v.P();

    for (int iteration = 0; iteration < 1; iteration++)
    {
      for (int i = 0; i < local_proj_samples.vert.size(); i++)
      {
        CVertex& v = local_proj_samples.vert[i];

        double min_clockwise_angle = 500;
        double max_anticlockwise_angle = -500;
        Point3f min_normal = v.N();
        Point3f max_normal = v.N();

        for (int j = 0; j < local_proj_samples.vert.size(); j++)
        {
          CVertex& t = local_proj_samples.vert[j];

          double angle = GlobalFun::computeDirectionalAngleOfTwoVertor(v.N(), t.N(), direction);
          if (angle > 0)
          {
            if (angle < min_clockwise_angle)
            {
              min_clockwise_angle = angle;
              min_normal = t.N();
            }
          }
          else
          {
            if (angle > max_anticlockwise_angle)
            {
              max_anticlockwise_angle = angle;
              max_normal = t.N();
            }
          }
        }

        Point3f new_normal = (min_normal + max_normal) / 2.0;
        new_normals[i] = new_normal.Normalize();

        if (i == 0)
        {
          test0_min.N() = min_normal;
          test1_max.N() = max_normal;
        }
      }

      for (int i = 0 ; i < new_normals.size(); i++)
      {
        CVertex& v = local_proj_samples.vert[i];
        v.N() = new_normals[i];
      }
    }

    new_new_normals[pick_index] =  new_normals[0];
    //samples->vert[pick_index].N() = new_normals[0];
  }

  for (int i = 0; i < samples->vert.size(); i++)
  {
    samples->vert[i].N() = new_new_normals[i];
  }

    //dual_samples->vert.clear();
    //dual_samples->vert.push_back(pick_v);
    //dual_samples->vert.push_back(test0_min);
    //dual_samples->vert.push_back(test1_max);
    //dual_samples->vn = dual_samples->vert.size();
}

// void WLOP::runProjection()
// {
//   cout << "runRegularizeSamples" << endl;
// 
//   double circle_radius = para->getDouble("CGrid Radius");
//   GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);
//   GlobalFun::computeEigenWithTheta(samples, para->getDouble("CGrid Radius") / sqrt(para->getDouble("H Gaussian Para")));
// 
//   //   int branch_KNN = para->getDouble("Branch Search KNN");
//   //   GlobalFun::computeAnnNeigbhors(samples->vert, samples->vert, branch_KNN, false, "void Skeletonization::searchNewBranches()");
// 
//   vector<Point3f> new_sample_set;
//   vector<CVertex> new_circle_points;
// 
//   float PI = 3.1415926;
// 
//   new_sample_set.assign(samples->vert.size(), Point3f(0, 0, 0));
//   for (int i = 0; i < samples->vert.size(); i++)
//   {
//     CVertex& v = samples->vert[i];
//     Point3f direction = v.eigen_vector0.Normalize();
// 
//     //Point3f front_nearest_p = v.P();
//     //Point3f back_nearest_p = v.P();
//     //double front_nearest_dist = 1000000.;
//     //double back_nearest_dist = -1000000.;
// 
//     //for (int j = 0; j < v.neighbors.size(); j++)
//     //{
//     //  int neighbor_idx = v.neighbors[j];
//     //  CVertex& t = samples->vert[neighbor_idx];
// 
//     //  Point3f diff = t.P() - v.P();
//     //  double proj_dist = diff * direction;
// 
//     //  if (proj_dist > 0)
//     //  {
//     //    if (proj_dist < front_nearest_dist)
//     //    {
//     //      front_nearest_p = t.P();
//     //      front_nearest_dist = proj_dist;
//     //    }
//     //  }
//     //  else
//     //  {
//     //    if (proj_dist > back_nearest_dist)
//     //    {
//     //      back_nearest_p = t.P();
//     //      back_nearest_dist = proj_dist;
//     //    }
//     //  }
//     //}
// 
//     //if (front_nearest_dist > 100000 || back_nearest_dist < -100000)
//     //{
//     //  cout << "end pointssss" << endl;
//     //  continue;
//     //}
// 
//     //v.P() = (front_nearest_p + back_nearest_p) / 2.0;
// 
//     //Point3f axis_direction = (front_nearest_p - v.P()).Normalize();
//     Point3f axis_direction = v.eigen_vector0.Normalize();
// 
//     Point3f random_p(0.1234, 0.1234, 0.1234);
//     Point3f random_direction = (random_p - v.P()).Normalize();
//     Point3f plane_dir0 = (random_direction ^ axis_direction).Normalize();
//     Point3f plane_dir1 = (plane_dir0 ^ axis_direction).Normalize();
//     Point3f axis_X = plane_dir0 * circle_radius;
//     Point3f axis_Y = plane_dir1 * circle_radius;
// 
//     float n = 100;
//     for (int j = 0; j < n; j++)
//     {
//       CVertex new_v;
//       //new_v.P() = v.P() + axis_Y * cos(2*PI/n*j) + axis_X * sin(2*PI/n*j); 
//       new_v.P() = v.P() + axis_Y* cos(2*PI/n*j) + axis_X * sin(2*PI/n*j); 
// 
//       new_v.N() = (new_v.P() - v.N()).Normalize();
//       new_v.m_index = new_circle_points.size();
//       new_circle_points.push_back(new_v);
//     }
//   }
// 
//   //dual_samples->vert.clear();
//   //for (int i = 0; i < new_circle_points.size(); i++)
//   //{
//   //  dual_samples->vert.push_back(new_circle_points[i]);
//   //}
//   //dual_samples->vn = dual_samples->vert.size();
// 
//   samples->vert.clear();
//   for (int i = 0; i < new_circle_points.size(); i++)
//   {
//     samples->vert.push_back(new_circle_points[i]);
//   }
//   samples->vn = samples->vert.size();
// 
//   return;
// }