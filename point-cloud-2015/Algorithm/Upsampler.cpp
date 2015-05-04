#include "Algorithm/Upsampler.h"
#include <cassert>

vector<double> Upsampler::proj_dist;
vector<double> Upsampler::proj_weight;

vector<Point3f> Upsampler::sum_N;
vector<Point3f> Upsampler::sum_Gw;
vector<Point3f> Upsampler::sum_Gf;


Upsampler::Upsampler(RichParameterSet* _para)
{
	samples = NULL;
	para = _para;
	b_first = true;
  also_insert_dual_points = false;
  use_adaptive_upsampling = false;
  do_upsample_on_skeltal_points = false;

	sigma = 25.0;
	old_radius = 0.0;
}

Upsampler::~Upsampler(void)
{
	samples = NULL;
}

void Upsampler::setInput(DataMgr* pData)
{
	if(!pData->isSamplesEmpty())
	{

		input(pData->getCurrentSamples());

    if (para->getBool("Use Upsample On Skeletal Points"))
    {
      input(pData->getCurrentSkelPoints());
    }


    //input(pData->getCurrentDualSamples());

    dual_samples = pData->getCurrentDualSamples();
    skel_points = pData->getCurrentSkelPoints();
	}
	else
	{
		cout << "ERROR: Upsampler:: setInput: empty!!" << endl;
		return;
	}
}

void Upsampler::run()
{
  also_insert_dual_points = global_paraMgr.glarea.getBool("Show Skeltal Points");
  use_adaptive_upsampling = global_paraMgr.upsampling.getBool("Use Adaptive Upsampling");

  if (para->getBool("Use Upsample On Skeletal Points"))
  {
    do_upsample_on_skeltal_points = true;
    also_insert_dual_points = false;
    use_adaptive_upsampling = false;
  }
  else
  {
    do_upsample_on_skeltal_points = false;
  }


	if(samples == NULL )
	{
		cout << "ERROR: Upsampler::run() mesh == NULL || original == NULL";
		return;
	}

	if (para->getBool("Run Projection"))
	{
		cout << "projection" <<endl;
		optimizeProjection();
		return;
  }

  if (para->getBool("Run Predict Constant Threshold"))
  {
    cout << "Run Predict Constant Threshold" << endl;
    runGetConstantPredictThreshold();
    return;
  }

	sigma = para->getDouble("Feature Sigma");
	G_value = para->getDouble("Edge Parameter");
	cout << "Edge Parameter: " << G_value << endl;

	if(b_first)
	{
		radius = para->getDouble("CGrid Radius");
	}

  if (para->getBool("Run Points Extrapolation"))
  {
    cout << "Extrapolation! " << endl;
    recomputeAllNeighbors();
    pointsExtrapoaltion();
    return;
  }


  if (also_insert_dual_points)
  {
    assert(samples->vert.size() == dual_samples->vert.size());
    assert(skel_points->vert.size() == dual_samples->vert.size());
  }

  if (para->getBool("Use Constant Threshold"))
  {
    radius = para->getDouble("CGrid Radius");
		double current_radius = radius;
		if (abs(old_radius - current_radius) > 1e-10 || b_first)
		{
			recomputeAllNeighbors();
			old_radius = current_radius;
		}
    clearAllThresholdFlag();
    runConstantUpsampling();
    return;
  }

	doUpsampling();	
}

// run 
void Upsampler::doUpsampling()
{
	radius = para->getDouble("CGrid Radius");


	// insert new points by dist threshold
	if(para->getBool("Using Threshold Process"))
	{
		double current_radius = para->getDouble("CGrid Radius");
		if( abs( old_radius - current_radius) > 1e-10 || b_first)
		{
			recomputeAllNeighbors();
			old_radius = current_radius;	
		}
		clearAllThresholdFlag();
		insertPointsByThreshold();
		return;
	}
	else
	{

	}


}

void Upsampler::clearAllThresholdFlag()
{
	is_abandon_by_threshold.assign(samples->vert.size(), false);
	CMesh::VertexIterator vi;
	for(vi = samples->vert.begin(); vi != samples->vert.end(); ++vi)
	{
		is_abandon_by_threshold[vi->m_index] = false;
	}
}


void Upsampler::input(CMesh* _mesh)
{
	if(_mesh == NULL)
	{
		cout << "ERROR:  Upsampler::setInput mesh == NULL!!" << endl;
		system("Pause");
		return;
	}
	else
	{
		samples = _mesh;

	}
}



void Upsampler::initVertexes()
{
	is_abandon_by_threshold.assign(samples->vert.size(), false);
	CMesh::VertexIterator vi, vi_end;
	vi_end = samples->vert.end();
	for(vi = samples->vert.begin(); vi != vi_end; ++vi) 
	{
		vi->neighbors.clear();
		is_abandon_by_threshold[vi->m_index] = false;
	}

}



// compute new vertex's neighbors 
void Upsampler::getNewPointNeighbors(CVertex & newv, set<int> & set_nbs, bool isBoth)
{
	double radius2 = radius * radius;

  //if (global_paraMgr.wLop.getBool("Use Adaptive Sample Neighbor"))
  //{
  //  radius2 /= 4;
  //}
	// 
	set<int>::iterator setIter,setEndIter;
	setEndIter = set_nbs.end();
	for (setIter = set_nbs.begin(); setIter != setEndIter; ++setIter)
	{
		CVertex & t = samples->vert[*setIter];
		Point3f diff = t.P() - newv.P();
		double dist2 = diff.SquaredNorm();

		if (dist2 < radius2)
		{
			newv.neighbors.push_back(t.m_index);

			// add it in its neighbor, this is for re-using this function to another purpose
			if(isBoth)
				t.neighbors.push_back(newv.m_index);
		}
	} // for
}


void Upsampler::removeOneNeighbor(vector<int>& neighbors, int remove_idx)
{
  int location = -1;
  for (int i = 0; i < neighbors.size(); i++)
  {
    if (neighbors[i] == remove_idx)
    {
      location = i;
      break;
    }
  }

  if (location < 0)
  {
    cout << "no way" << endl;
    system("Pause");
  }

  neighbors.erase(neighbors.begin() + location);
}

void Upsampler::getLineVertNeighorsIndex(set<int> & neighbors, int firstV, int secV)
{
	neighbors.clear();

	CVertex & v = samples->vert[firstV];
	CVertex & t = samples->vert[secV];

//  if (global_paraMgr.wLop.getBool("Use Adaptive Sample Neighbor"))
//  {
////     removeOneNeighbor(v.neighbors, t.m_index);
////     removeOneNeighbor(t.neighbors, v.m_index);
//  }

	vector<int>::iterator iter;

	// visit v's neighbors
	for (iter = v.neighbors.begin(); iter !=  v.neighbors.end(); ++iter)
	{
		neighbors.insert(*iter);
	}

	for (iter = t.neighbors.begin(); iter != t.neighbors.end(); ++iter)
	{
		neighbors.insert(*iter);
	}

	// line's vertex
	neighbors.insert(firstV);
	neighbors.insert(secV);

} // end




void Upsampler::computeEigenVerctorForRendering()
{
	for(int i = 0; i < samples->vert.size(); i++) 
	{
		CVertex &v = samples->vert[i];
		v.recompute_m_render();
	}
}



// find the max dist of minimum dist between midpoint and their neighbors,return its index using para
double Upsampler::findMaxMidpoint(CVertex & v, int & neighbor_index)
{	
	bool trick_for_closeby_surface = para->getBool("Use Upsample For Closeby Surface");

	double G_value = para->getDouble("Edge Parameter");
	int nb_size = v.neighbors.size();
	double bestDist = -1; //

	Point3f midPoint = Point3f(0.0,0.0,0.0);
	Point3f diff = Point3f(0.0,0.0,0.0);
	double proj_s = 0.0;
	double dist2 = 0.0;
	double minDist; // minimum dist of the current midpoint

	for (int i = 0; i < nb_size; ++i)
	{
		CVertex & t = samples->vert[v.neighbors[i]];

    if (use_adaptive_upsampling && !t.is_fixed_sample)
    {
      continue;
    }

		midPoint = (v.P() + t.P()) / 2.0;

		if (trick_for_closeby_surface)
		{
			//if (GlobalFun::computeRealAngleOfTwoVertor(v.N(), t.N()) > 90 )
			if (v.N() * t.N() < 0 )
			{
				//cout << "skip" << endl ;
				continue;
			}
		}

		Point3f vm = v.N();
		Point3f tm = t.N();
		double dot_produce = pow((2.0 - vm * tm), G_value);

		Point3f diff_md_t = (midPoint - t.P());
		double proj = diff_md_t * t.N();

		minDist = diff_md_t.SquaredNorm() - proj * proj;

		for (int j = 0; j < nb_size; ++j)
		{
			CVertex & s = samples->vert[v.neighbors[j]];
			diff = midPoint - s.P();
			proj_s = diff * s.N();

			dist2 = diff.SquaredNorm();
			double proj_min = dist2 - proj_s * proj_s;

			if(proj_min < minDist)
			{
				minDist = proj_min;
			}
		}
		minDist *= dot_produce;

		if(minDist > bestDist)
		{
			bestDist = minDist;
			neighbor_index = v.neighbors[i];
		}
	} // for

	return bestDist;// for dist threshold
} // end





//double Upsampler::getPredictThreshold()
//{
//	int testNumber = 50;
//	float sumDist = 0.0;
//	for(int i = 0; i < samples->vert.size() && i < testNumber; i++)
//	{
//		CVertex& v = samples->vert[i];
//		int secondVertIndx = -1;
//		sumDist += findMaxMidpoint(v, secondVertIndx);
//	}
//
//	return (sumDist / testNumber) * 0.8;
//}


double Upsampler::getPredictThresholdFirstTime()
{

	if(samples == NULL )
	{
		cout << "ERROR: Upsampler::run() mesh == NULL || original == NULL";
		return -1;
	}

	sigma = para->getDouble("Feature Sigma");
	G_value = para->getDouble("Edge Parameter");
	cout << "Edge Parameter: " << G_value << endl;

	if(b_first)
	{
		radius = para->getDouble("CGrid Radius");
	}

	radius = para->getDouble("CGrid Radius");


	double current_radius = para->getDouble("CGrid Radius");
	if( abs( old_radius - current_radius) > 1e-10 || b_first)
	{
		recomputeAllNeighbors();
		old_radius = current_radius;	
	}
	clearAllThresholdFlag();


	int testNumber = 50;
	double sumDist = 0.0;
	int count = 0;
	for(int i = 0; i < samples->vert.size() && i < testNumber; i++)
	{
		CVertex& v = samples->vert[i];
		int secondVertIndx = -1;
		sumDist += findMaxMidpoint(v, secondVertIndx);
		count++;
	}

	return (sumDist / count);
}
// end


void Upsampler::runGetConstantPredictThreshold()
{
  if (!global_paraMgr.wLop.getBool("Use Adaptive Sample Neighbor"))
  {
    GlobalFun::computeBallNeighbors(samples, NULL, 
      para->getDouble("CGrid Radius"), samples->bbox);
  }

  float sumDist = 0.0;
  int testNumber = 0;

  for(int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];
    int secondVertIndx = -1;
    double density = findMaxMidpoint(v, secondVertIndx);
    if (density < 0)
    {
      continue;
    }
    sumDist += density;
    testNumber++;
  }

  para->setValue("Dist Threshold", DoubleValue((sumDist / testNumber)*2.));
}

double Upsampler::getPredictThreshold()
{
  float sumDist = 0.0;
  int testNumber = 0;

  for(int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];
    int secondVertIndx = -1;
    double dist = findMaxMidpoint(v, secondVertIndx);
    if (dist < 0)
    {
      continue;
    }
    sumDist += dist;
    testNumber++;
  }

  return (sumDist / testNumber) * 0.8;
}

void Upsampler::runConstantUpsampling()
{
  int abandonCounter = 0;
  int loopCounter = 0;
  double addCounter = 0;
  int MAX_LOOP = 1;

  while(1)
  {
    abandonCounter = 0;
    loopCounter++;

    int print_threshold_num = 50;

    double dist_threshold = para->getDouble("Dist Threshold");

    cout << "threshold: " << dist_threshold << endl;

    samples->vn = samples->vert.size();

    for(int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];
      if(is_abandon_by_threshold[v.m_index])
      {
        abandonCounter ++;
        continue;
      }

			if (i < 5)
			{
				cout << "neighbor size: " << v.neighbors.size() << endl;;
			}

      int firstVertInex = v.m_index;
      int secondVertIndx = -1;
      double bestDist = findMaxMidpoint(v, secondVertIndx);
      //bestDist = sqrt(bestDist);

      if(secondVertIndx < 0)
        continue;

      if(bestDist < dist_threshold)
      {
        if(i < print_threshold_num)
          cout << bestDist << endl;

        is_abandon_by_threshold[v.m_index] = true;
        continue;
      }

      //assert(secondVertIndx >= 0 && secondVertIndx < samples->vert.size());
      CVertex & t = samples->vert[secondVertIndx];

      //assert(secondVertIndx == t.m_index);  

      CVertex newv;
      newv.P() = (v.P() + t.P()) / 2.0;  // mass center
      newv.m_index = samples->vert.size();

      if (do_upsample_on_skeltal_points)
      {
        //newv.is_dual_sample = true;
        newv.is_skel_point = true;
      }

      if (also_insert_dual_points)
      {
        CVertex& dual_v = dual_samples->vert[v.m_index];
        CVertex& dual_t = dual_samples->vert[t.m_index];

        CVertex new_dual_v;
        new_dual_v.P() = (dual_v.P() + dual_t.P()) / 2.0;  // mass center
        new_dual_v.N() = (dual_v.N() + dual_t.N()) / 2.0;  // mass center

        new_dual_v.m_index = dual_samples->vert.size();
        new_dual_v.dual_index = new_dual_v.m_index;
        new_dual_v.is_dual_sample = true;

        // add new point
        dual_samples->vert.push_back(new_dual_v);

        CVertex& skel_v = skel_points->vert[v.m_index];
        CVertex& skel_t = skel_points->vert[t.m_index];

        CVertex new_skel_v;
        new_skel_v.P() = (skel_v.P() + skel_t.P()) / 2.0;  // mass center
        new_skel_v.N() = (skel_v.N() + skel_t.N()) / 2.0;  // mass center

        new_skel_v.m_index = skel_points->vert.size();
        new_skel_v.dual_index = new_skel_v.m_index;
        new_skel_v.is_skel_point = true;

        // add new point
        skel_points->vert.push_back(new_skel_v);
      }

      // add new point
      samples->vert.push_back(newv);
      is_abandon_by_threshold.push_back(false);
      addCounter++;

      // get its temporary neighbors to the set
      set<int> setNewVertexNeighors;
      getLineVertNeighorsIndex(setNewVertexNeighors, firstVertInex, secondVertIndx);
      CVertex & pv = samples->vert[samples->vert.size()-1];

      // compute new vertex's temporary neighbors, its will update until project it
      getNewPointNeighbors(pv, setNewVertexNeighors, false);


      //computeNewVertexNormAvgMethod(pv, firstVertInex, secondVertIndx);// 7-11 not sure
      computeNewVertexProjDist_Sigma(pv, firstVertInex, secondVertIndx);


      if(wd != 0.0)
        pv.P() = pv.P() + pv.N() * (d / wd);

			pv.is_fixed_sample = true;

      // update new vertex's neighbors
      pv.neighbors.clear();
      getNewPointNeighbors(pv, setNewVertexNeighors, true);

			// new 2014-12-23
 			computeNewVertexProjDist_Again(pv);
 			if (wd != 0.0)
 				pv.P() = pv.P() + pv.N() * (d / wd);
    }

    if(abandonCounter == samples->vert.size() || loopCounter > MAX_LOOP)
      break;
    else
    {
      cout << "current size: " << samples->vert.size() << endl;
      cout << "abandent : " << abandonCounter << endl << endl; 
    }
  }

  samples->vn = samples->vert.size();


  computeEigenVerctorForRendering();
}


void Upsampler::pointsExtrapoaltion()
{
  int sample_number = samples->vert.size();
  for (int i = 0; i < sample_number; i++)
  {
    CVertex& v = samples->vert[i];
    v.is_boundary = false;
    v.m_index = i;
  }

  CVertex new_v;
  new_v.is_boundary = true;


  if (global_paraMgr.glarea.getDouble("Picked Index") > 0.1)
  {
    int pick_index = global_paraMgr.glarea.getDouble("Picked Index");
    cout  << "pick index: " << pick_index <<endl;

    if (pick_index >= sample_number)
    {
      return;
    }

    CVertex& v = samples->vert[pick_index];

    Point3f movement(0., 0., 0.);
    for (int j = 0; j < v.neighbors.size(); j++)
    {
      CVertex& t = samples->vert[v.neighbors[j]];
      Point3f diff = v.P() - t.P();
      movement += diff;
    }
    movement  = movement / v.neighbors.size();

    new_v.P() = v.P() + movement;
    new_v.N() = v.N();

    

    new_v.m_index = sample_number;
    samples->vert.push_back(new_v);
    samples->vn = samples->vert.size();

    recomputeAllNeighbors();
    projectOneCVertex(samples->vert[samples->vert.size()-1]);
  }
  else
  {


    vector<CVertex> extrapolation_points(sample_number, new_v);
    vector<double> extrapolation_dist(sample_number, 0.);

    for (int i = 0; i < sample_number; i++)
    {
      CVertex& v = samples->vert[i];

      Point3f movement(0., 0., 0.);
      for (int j = 0; j < v.neighbors.size(); j++)
      {
        CVertex& t = samples->vert[v.neighbors[j]];
        Point3f diff = v.P() - t.P();
        movement += diff;
      }

      CVertex& ext_v = extrapolation_points[i];
      ext_v.P() = v.P() + movement;
      ext_v.N() = v.N();
    }

    GlobalFun::computeAnnNeigbhors(samples->vert, extrapolation_points, 2, false, "WlopParaDlg::applyDualConnection()");

    for (int i = 0; i < sample_number; i++)
    {
      CVertex& ext_v = extrapolation_points[i];
      CVertex& ext_nearest_v = samples->vert[ext_v.neighbors[0]];
      double dist = GlobalFun::computeEulerDist(ext_v.P(),ext_nearest_v.P());
      extrapolation_dist[i] = dist;
    }

    int max_index = 0;
    double max_dist = 0.0;
    for (int i = 0; i < sample_number; i++)
    {
      if (extrapolation_dist[i] > max_dist)
      {
        max_dist = extrapolation_dist[i];
        max_dist = i;
      }
    }

    new_v = extrapolation_points[max_index];

    new_v.m_index = sample_number;
    samples->vert.push_back(new_v);
    samples->vn = samples->vert.size();
  }



}


void Upsampler::improveNeighborhood()
{
  for (int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];
    CVertex& dual_v = samples->vert[v.dual_index];

    vector<int> new_neighborhood;
    for (int j = 0; j < v.neighbors.size(); j++)
    {
      CVertex& t = samples->vert[v.neighbors[j]];
      CVertex& dual_t = dual_samples->vert[t.dual_index];

      double dist1 = GlobalFun::computeEulerDist(v.P(), t.P());
      double dist2 = GlobalFun::computeEulerDist(dual_v.P(), dual_t.P());
      if (dist2 != (dist1 * 1.5))
      {
        cout << "not good" << endl;
        continue;
      }

      new_neighborhood.push_back(v.neighbors[j]);
    }

    v.neighbors = new_neighborhood;

  }

}


void Upsampler::insertPointsByThreshold()
{
	cout << "insertPointsByThreshold" << endl;
  int abandonCounter = 0;
  int loopCounter = 0;
  double addCounter = 0;
  int MAX_LOOP = 5;
  int oldSize = samples->vert.size();

  int max_add_number = para->getInt("Number of Add Point");


  if (use_adaptive_upsampling)
  {
    improveNeighborhood();
  }


  while(1)
  {
    abandonCounter = 0;
    loopCounter++;

    int print_threshold_num = 50;

    double dist_threshold = para->getDouble("Dist Threshold");

    cout << "threshold: " << dist_threshold << endl;

    samples->vn = samples->vert.size();

    for(int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];

      if (use_adaptive_upsampling && !v.is_fixed_sample)
      {
        continue;
      }

      if (v.neighbors.empty())
      {
        cout << "empty neighbors!!!" << endl;
        continue;
      }

			if (i < 15)
			{
				cout << "neighbor size: " << v.neighbors.size();
			}

      if(is_abandon_by_threshold[v.m_index])
      {
        abandonCounter++;
        continue;
      }

      int firstVertInex = v.m_index;
      int secondVertIndx = -1;
      double bestDist = findMaxMidpoint(v, secondVertIndx);
      //bestDist = sqrt(bestDist);

      if(secondVertIndx < 0)
        continue;

      if(bestDist < dist_threshold)
      {
        if(i < print_threshold_num)
          cout << bestDist << endl;

        is_abandon_by_threshold[v.m_index] = true;
        continue;
      }

      //assert(secondVertIndx >= 0 && secondVertIndx < samples->vert.size());
      CVertex & t = samples->vert[secondVertIndx];

      //assert(secondVertIndx == t.m_index);  

      CVertex newv;
      newv.P() = (v.P() + t.P()) / 2.0;  // mass center
      newv.m_index = samples->vert.size();
      newv.dual_index = newv.m_index;

      if (use_adaptive_upsampling)
      {
        newv.is_fixed_sample = true;
      }

      if (do_upsample_on_skeltal_points)
      {
        //newv.is_dual_sample = true;
        newv.is_skel_point = true;
      }

      if (also_insert_dual_points)
      {
        CVertex& dual_v = dual_samples->vert[v.m_index];
        CVertex& dual_t = dual_samples->vert[t.m_index];

        CVertex new_dual_v;
        new_dual_v.P() = (dual_v.P() + dual_t.P()) / 2.0;  // mass center
        new_dual_v.N() = (dual_v.N() + dual_t.N()) / 2.0;  // mass center

        new_dual_v.m_index = dual_samples->vert.size();
        new_dual_v.dual_index = new_dual_v.m_index;
        new_dual_v.is_dual_sample = true;

        // add new point
        dual_samples->vert.push_back(new_dual_v);




        CVertex& skel_v = skel_points->vert[v.m_index];
        CVertex& skel_t = skel_points->vert[t.m_index];
        
        CVertex new_skel_v;
        new_skel_v.P() = (skel_v.P() + skel_t.P()) / 2.0;  // mass center
        new_skel_v.N() = (skel_v.N() + skel_t.N()) / 2.0;  // mass center
        
        new_skel_v.m_index = skel_points->vert.size();
        new_skel_v.dual_index = new_skel_v.m_index;
        new_skel_v.is_skel_point = true;

        // add new point
        skel_points->vert.push_back(new_skel_v);
      }



      // add new point
      samples->vert.push_back(newv);
      is_abandon_by_threshold.push_back(false);
      addCounter++;

      // get its temporary neighbors to the set
      set<int> setNewVertexNeighors;
      getLineVertNeighorsIndex(setNewVertexNeighors, firstVertInex, secondVertIndx);
      CVertex & pv = samples->vert[samples->vert.size()-1];

      // compute new vertex's temporary neighbors, its will update until project it
      getNewPointNeighbors(pv, setNewVertexNeighors, false);


      //computeNewVertexNormAvgMethod(pv, firstVertInex, secondVertIndx);// 7-11 not sure
      computeNewVertexProjDist_Sigma(pv, firstVertInex, secondVertIndx);


      if(wd != 0.0)
        pv.P() = pv.P() + pv.N() * (d / wd);

      // update new vertex's neighbors
      pv.neighbors.clear();
      getNewPointNeighbors(pv, setNewVertexNeighors, true);

      if (addCounter > max_add_number)
      {
        break;
      }
    }

    if(abandonCounter == samples->vert.size() || loopCounter > MAX_LOOP || addCounter > max_add_number)
      break;
    else
    {
      cout << "current size: " << samples->vert.size() << endl;
      cout << "add point: " << samples->vert.size() - oldSize << endl;
      oldSize = samples->vert.size();
      cout << "abandent : " << abandonCounter << endl << endl; 
      cout << "dual sample number: " << dual_samples->vert.size() << endl;

    }

  }

  cout << endl  << "dual sample number: " << dual_samples->vert.size() << 
       "skel number" << skel_points->vert.size() << endl;

  samples->vn = samples->vert.size();

  if (addCounter > max_add_number)
  {
    cout << "exeed max add number" << endl;
  }
  para->setValue("Dist Threshold", DoubleValue(getPredictThreshold()));
  cout << "getPredictThreshold" << getPredictThreshold() << endl;

  computeEigenVerctorForRendering();

}





void Upsampler::recomputeAllNeighbors()
{
	double grid_radius = para->getDouble("CGrid Radius");
	radius = para->getDouble("CGrid Radius");

	cout << "recomputeAllNeighbors" << endl;

	initVertexes();

	cout << "radius: " << grid_radius << endl;

	GlobalFun::computeBallNeighbors(samples, NULL, 
		para->getDouble("CGrid Radius"), samples->bbox);

	bool trick_for_closeby_surface = para->getBool("Use Upsample For Closeby Surface");
	
	if (trick_for_closeby_surface)
	{
		for (int i = 0; i < samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];
			vector<int> new_neighbors;
			for (int j = 0; j < v.neighbors.size(); j++)
			{
				CVertex& t = samples->vert[v.neighbors[j]];

				if (v.N() * t.N() > 0)
				{
					new_neighbors.push_back(v.neighbors[j]);
				}
			}

			v.neighbors = new_neighbors;
			if (v.neighbors.empty())
			{
				cout << "become empty neighbor" << endl;
			}
		}
	}

	cout << "recomputeAllNeighbors end "<< endl;
}



void Upsampler::computeNewVertexNormAvgMethod(CVertex & v, int firstV, int secV)// add by wsh 1-5
{
	Point3f sum_theta1_multiply_normal = Point3f(0.0, 0.0, 0.0);
	double sum_theta1 = 0.0;

	double radius2 = radius * radius;
	double radius16 =   -4 / (radius2);

	Point3f & p = v.P();

	for(int i = 0; i < v.neighbors.size(); i++)
	{
		CVertex & t =  samples->vert[v.neighbors[i]];

		Point3f & q = t.P();
		double dist2 = (q - p).SquaredNorm();

		double theta_1 = exp(dist2 * radius16);

		sum_theta1 += theta_1;
		sum_theta1_multiply_normal += t.N() * theta_1;
	}

	Point3f avg_normal = sum_theta1_multiply_normal / sum_theta1;

	v.N() = avg_normal.Normalize();
}


void Upsampler::projectOneCVertex(CVertex& v)
{
  double radius2 = radius * radius;
  double radius16 =   -4 / (radius2);

  int neighbor_number = v.neighbors.size();
  //vector<double> sum_weight(neighbor_number, 0.0);
  //vector<double> sum_proj_dist(neighbor_number, 0.0);
  //vector<Point3f> sum_normal(neighbor_number, Point3f(0, 0, 0));

  double sum_weight = 0.0;
  double sum_proj_dist = 0.0;
  Point3f sum_normal(Point3f(0, 0, 0));

  if (v.neighbors.empty())
  {
    return;
  }

  for (int i = 0; i < neighbor_number; i++)
  {
    CVertex& t = samples->vert[v.neighbors[i]];
    
    double theta_1 = exp( (t.P() - v.P()).SquaredNorm() * radius16);
    double psi = exp(-pow(1- v.N() * t.N(), 2)/pow(max(1e-8,1-cos(sigma/180.0*3.1415926)), 2));

    double diff_proj_vt = (v.P() - t.P()) * t.N();
    double weight = psi * theta_1;

    sum_proj_dist +=  weight * diff_proj_vt;
    sum_weight += weight;
    sum_normal +=  t.N() * weight;
  }

  double proj_dist = sum_proj_dist / sum_weight;
  v.N() = sum_normal / sum_weight;
  v.P() -= v.N() * proj_dist;

  //pv.P() = pv.P() + pv.N() * (d / wd);
}

//our new sigma
//void Upsampler::computeNewVertexProjDist_Sigma(CVertex & v, int firstV, int secV)
//{
//	const int NORMAL_NUM = 2;
//	double radius2 = radius * radius;
//	double radius16 =   -4 / (radius2);
//	Point3f & p = v.P();
//	//Point3f vm(v.N());
//
//
//	// 0:first, 1:second
//	vector<double> projDist(NORMAL_NUM, 0.0);
//	vector<double> sumW(NORMAL_NUM, 0.0);
//	vector<Point3f> v_normals;
//
//	vector<Point3f> grad_f(NORMAL_NUM, Point3f(0., 0., 0.));
//
//	v_normals.push_back(samples->vert[firstV].N());
//	v_normals.push_back(samples->vert[secV].N());
//
//	for(int i = 0; i < v.neighbors.size(); i++)
//	{
//		CVertex & t =  samples->vert[v.neighbors[i]];
//		Point3f & q = t.P();
//
//		Point3f tm(t.N());
//		double theta_1 = exp( (q - p).SquaredNorm() * radius16);
//
//		for(int j = 0; j < NORMAL_NUM; j++)
//		{
//			double psi = exp(-pow(1- v_normals[j] * tm, 2)/pow(max(1e-8,1-cos(sigma/180.0*3.1415926)), 2));
//
//			double diff_proj_vt = (p - q) * tm;
//			double weight = psi * theta_1;
//
//			projDist[j] +=  weight * diff_proj_vt;
//			sumW[j] += weight;
//			grad_f[j] +=  tm * weight;
//		}
//	}
//
//	double min_dist = 10000.0;
//	int min_index = 0;
//	for(int i = 0; i < NORMAL_NUM; i++)
//	{
//		double dist = abs((projDist[i] / sumW[i]));
//		if(dist < min_dist)
//		{
//			min_dist = dist;
//			min_index = i;
//		}
//	}
//
//	d = -projDist[min_index] / sumW[min_index];
//	wd = 1.0;
//
//	v.N() = (grad_f[min_index] / sumW[min_index]).Normalize();
//}


void Upsampler::computeNewVertexProjDist_Again(CVertex & v)
{
	const int NORMAL_NUM = 2;
	double radius2 = radius * radius;
	double radius16 = -4 / (radius2);
	Point3f & p = v.P();

	double projDist(0.0);
	double sumW(0.0);
	Point3f v_normals = v.N();

	Point3f grad_f(Point3f(0., 0., 0.));

	for (int i = 0; i < v.neighbors.size(); i++)
	{
		CVertex & t = samples->vert[v.neighbors[i]];
		Point3f & q = t.P();

		Point3f tm(t.N());
		double theta_1 = exp((q - p).SquaredNorm() * radius16);

		for (int j = 0; j < NORMAL_NUM; j++)
		{
			double psi = exp(-pow(1 - v_normals * tm, 2) / pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2));

			double diff_proj_vt = (p - q) * tm;
			double weight = psi * theta_1;

			projDist += weight * diff_proj_vt;
			sumW += weight;
			grad_f += tm * weight;
		}
	}

	d = -projDist/ sumW;
	wd = 1.0;

	v.N() = (grad_f / sumW).Normalize();
}


void Upsampler::computeNewVertexProjDist_Sigma(CVertex & v, int firstV, int secV)
{
	const int NORMAL_NUM = 2;
	double radius2 = radius * radius;
	double radius16 = -4 / (radius2);
	Point3f & p = v.P();
	//Point3f vm(v.N());


	// 0:first, 1:second
	vector<double> projDist(NORMAL_NUM, 0.0);
	vector<double> sumW(NORMAL_NUM, 0.0);
	vector<Point3f> v_normals;

	vector<Point3f> grad_f(NORMAL_NUM, Point3f(0., 0., 0.));

	v_normals.push_back(samples->vert[firstV].N());
	v_normals.push_back(samples->vert[secV].N());

	for (int i = 0; i < v.neighbors.size(); i++)
	{
		CVertex & t = samples->vert[v.neighbors[i]];
		Point3f & q = t.P();

		Point3f tm(t.N());
		double theta_1 = exp((q - p).SquaredNorm() * radius16);

		for (int j = 0; j < NORMAL_NUM; j++)
		{
			double psi = exp(-pow(1 - v_normals[j] * tm, 2) / pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2));

			double diff_proj_vt = (p - q) * tm;
			double weight = psi * theta_1;

			projDist[j] += weight * diff_proj_vt;
			sumW[j] += weight;
			grad_f[j] += tm * weight;
		}
	}

	double min_dist = 10000.0;
	int min_index = 0;
	for (int i = 0; i < NORMAL_NUM; i++)
	{
		double dist = abs((projDist[i] / sumW[i]));
		if (dist < min_dist)
		{
			min_dist = dist;
			min_index = i;
		}
	}

	d = -projDist[min_index] / sumW[min_index];
	wd = 1.0;

	v.N() = (grad_f[min_index] / sumW[min_index]).Normalize();
}



void Upsampler::optimizeProjection()
{
	proj_dist.assign(samples->vn, 0);
	proj_weight.assign(samples->vn, 0);

	sum_N.assign(samples->vn, Point3f(0, 0, 0));
	sum_Gw.assign(samples->vn, Point3f(0, 0, 0));
	sum_Gf.assign(samples->vn, Point3f(0, 0, 0));

	double radius = para->getDouble("CGrid Radius");
	CGrid mesh_grid;
	CGrid original_grid;
	//mesh_grid.init(samples->vert, samples->bbox, radius);

	double radius2 = radius * radius;
	double iradius16 = -4 / radius2;

	/* updateNormal Before projection */
	//mesh_grid.iterate(self_optimize_normals, other_optimize_normals);
	//mesh_grid.iterate(self_optimize_upsamples_global, other_optimize_upsamples_global);

	GlobalFun::computeBallNeighbors(samples, NULL, para->getDouble("CGrid Radius"), samples->bbox);

	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];

		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& t = samples->vert[v.neighbors[j]];

			Point3f diff = v.P() - t.P();  // 
			double dist2 = diff.SquaredNorm();

			Point3f vm(v.N());
			Point3f tm(t.N());

				double theta_1 = exp(dist2 * iradius16);
				double psi = exp(-pow(1 - vm*tm, 2) / pow(max(1e-18, 1 - cos(sigma / 180.0*3.1415926)), 2));

				double weight = theta_1 * psi;
				updateVT_proj_normal(v, t, weight, radius / 2.);
		}
	}

	for (int i = 0; i < samples->vn; i++)
	{
		CVertex& v = samples->vert[i];
    if (v.neighbors.empty())
    {
      continue;
    }
		if (proj_weight[i] > 1e-20 && !v.neighbors.empty())
		{
			v.N() = sum_N[i] / proj_weight[i];
			v.N().Normalize();
			v.P() -= v.N() * (proj_dist[i] / proj_weight[i]);
		}
	}
}

void Upsampler::self_optimize_upsamples_global(CGrid::iterator start, CGrid::iterator end, double radius)
{
	double radius2 = radius*radius;
	double iradius16 = -4/radius2;  
	const double PI = 3.1415926;
	double delta = 2.0;
	double sigma = global_paraMgr.upsampling.getDouble("Feature Sigma");

	for(CGrid::iterator dest = start; dest != end; dest++) {
		CVertex &v = *(*dest);
		Point3f &p = v.P();
		for(CGrid::iterator origin = dest+1; origin != end; origin++) {
			CVertex &t = *(*origin);
			Point3f &q = t.P();
			Point3f diff = q-p;  // 
			double dist2 = diff.SquaredNorm();

			if(dist2 < radius2) 
			{    
				Point3f vm(v.N());
				Point3f tm(t.N());

				double theta_1 = exp(dist2 * iradius16); 
				double psi = exp(-pow(1-vm*tm, 2)/pow(max(1e-18,1-cos(sigma/180.0*3.1415926)), 2));

				double weight = theta_1 * psi;
				updateVT_proj(v, t, weight);
				updateVT_proj(t, v, weight);			
			}                         
		}
	}
}

void Upsampler::other_optimize_upsamples_global(CGrid::iterator starta, CGrid::iterator enda, 
														CGrid::iterator startb, CGrid::iterator endb, double radius)
{
	double radius2 = radius*radius;
	double iradius16 = -4/radius2;  
	const double PI = 3.1415926;
	double delta = 2.0;
	double sigma = global_paraMgr.upsampling.getDouble("Feature Sigma");

	for(CGrid::iterator dest = starta; dest != enda; dest++) {
		CVertex &v = *(*dest);
		Point3f &p = v.P();
		for(CGrid::iterator origin = startb; origin != endb; origin++) {
			CVertex &t = *(*origin);
			Point3f &q = t.P();
			Point3f diff = q-p;  // 
			double dist2 = diff.SquaredNorm();

			if(dist2 < radius2) {     
				Point3f vm(v.N());
				Point3f tm(t.N());

				double theta_1 = exp(dist2 * iradius16); 
				double psi = exp(-pow(1-vm*tm, 2)/pow(max(1e-18,1-cos(sigma/180.0*3.1415926)), 2));

				double weight = theta_1 * psi;
				updateVT_proj(v, t, weight);
				updateVT_proj(t, v, weight);
			}                         
		}
	}
}

void Upsampler::updateVT_proj(CVertex& v, CVertex& t, double weight)
{
	double diff_proj_vt = v.N() * (t.P() - v.P());
	proj_dist[v.m_index] += diff_proj_vt * weight;
	proj_weight[v.m_index] += weight;
}


void Upsampler::self_optimize_normals(CGrid::iterator start, CGrid::iterator end, double radius)
{
	double radius2 = radius*radius;
	double iradius16 = -4/radius2;  
	const double PI = 3.1415926;
	double delta = 2.0;
	double sigma = global_paraMgr.upsampling.getDouble("Feature Sigma");

	for(CGrid::iterator dest = start; dest != end; dest++) 
  {
		CVertex &v = *(*dest);
		Point3f &p = v.P();
		for(CGrid::iterator origin = dest+1; origin != end; origin++)
		{
			CVertex &t = *(*origin);
			Point3f &q = t.P();
			Point3f diff = q-p;  // 
			double dist2 = diff.SquaredNorm();

			if(dist2 < radius2) 
			{    
				Point3f vm(v.N());
				Point3f tm(t.N());

				double theta_1 = exp(dist2 * iradius16); 
				double psi = exp(-pow(1-vm*tm, 2)/pow(max(1e-18,1-cos(sigma/180.0*3.1415926)), 2));

				double weight = theta_1 * psi;
				updateVT_proj_normal(v, t, weight, radius / 2.);
				updateVT_proj_normal(t, v, weight, radius / 2.);
			}                         
		}
	}
}

void Upsampler::other_optimize_normals(CGrid::iterator starta, CGrid::iterator enda, 
	CGrid::iterator startb, CGrid::iterator endb, double radius)
{
	double radius2 = radius*radius;
	double iradius16 = -4/radius2;  
	const double PI = 3.1415926;
	double delta = 2.0;
	double sigma = global_paraMgr.upsampling.getDouble("Feature Sigma");

	for(CGrid::iterator dest = starta; dest != enda; dest++) {
		CVertex &v = *(*dest);
		Point3f &p = v.P();
		for(CGrid::iterator origin = startb; origin != endb; origin++) 
		{
			CVertex &t = *(*origin);
			Point3f &q = t.P();
			Point3f diff = q-p;  // 
			double dist2 = diff.SquaredNorm();

			if(dist2 < radius2) {   
				Point3f vm(v.N());
				Point3f tm(t.N());

				double theta_1 = exp(dist2 * iradius16); 
				double psi = exp(-pow(1-vm*tm, 2)/pow(max(1e-18,1-cos(sigma/180.0*3.1415926)), 2));

				double weight = theta_1 * psi;
				updateVT_proj_normal(v, t, weight, radius / 2.);
				updateVT_proj_normal(t, v, weight, radius / 2.);
			}                         
		}
	}
}

void Upsampler::updateVT_proj_normal(CVertex& v, CVertex& t, double weight, double radius)
{
	double radius2 =  radius * radius;
	double diff_proj_vt = (v.P() - t.P()) * t.N();
	proj_dist[v.m_index]  += diff_proj_vt * weight;
	proj_weight[v.m_index] += weight;

	Point3f grad_w = (v.P() - t.P()) * weight * (-2.) / radius2;

	sum_Gw[v.m_index] += grad_w;
	sum_Gf[v.m_index] += grad_w * diff_proj_vt;
	sum_N[v.m_index] += t.N() * weight;
}

