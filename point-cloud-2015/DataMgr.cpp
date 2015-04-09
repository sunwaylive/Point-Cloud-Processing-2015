#include "DataMgr.h"


DataMgr::DataMgr(RichParameterSet* _para)
{
	para = _para;
}


DataMgr::~DataMgr(void)
{
}

void DataMgr::clearCMesh(CMesh& mesh)
{
	mesh.face.clear();
	mesh.fn = 0;
	mesh.vert.clear();
	mesh.vn = 0;
	mesh.bbox = Box3f();
}

bool DataMgr::isSamplesEmpty()
{
	return samples.vert.empty();
}


bool DataMgr::isOriginalEmpty()
{
	return original.vert.empty();
}

bool DataMgr::isSkeletonEmpty()
{
  return skeleton.isEmpty();
}

bool DataMgr::isDualSamplesEmpty()
{
	return dual_samples.vert.empty();
}

bool DataMgr::isSkeletalPointsEmpty()
{
	return skel_points.vert.empty();
}


void DataMgr::loadPlyToOriginal(QString fileName)
{
	clearCMesh(original);
	curr_file_name = fileName;

	int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;

	int err = tri::io::Importer<CMesh>::Open(original, curr_file_name.toStdString().data(), mask);  
	if(err) 
	{
		cout << "Failed reading mesh: " << err << "\n";
		return;
	}  
	cout << "points loaded\n";


	CMesh::VertexIterator vi;
	int idx = 0;
	for(vi = original.vert.begin(); vi != original.vert.end(); ++vi)
	{
		vi->bIsOriginal = true;
		vi->m_index = idx++;
		//vi->N() = Point3f(0, 0, 0);
		original.bbox.Add(vi->P());
	}
	original.vn = original.vert.size();
}

void DataMgr::loadPlyToDualSample(QString fileName)
{
	clearCMesh(dual_samples);
	curr_file_name = fileName;

	int mask = tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL;

	int err = tri::io::Importer<CMesh>::Open(dual_samples, curr_file_name.toStdString().data(), mask);
	if (err)
	{
		cout << "Failed reading mesh: " << err << "\n";
		return;
	}
	cout << "points loaded\n";


	CMesh::VertexIterator vi;
	int idx = 0;
	for (vi = dual_samples.vert.begin(); vi != dual_samples.vert.end(); ++vi)
	{
		vi->is_dual_sample = true;
		vi->m_index = idx++;
		//vi->N() = Point3f(0, 0, 0);
		dual_samples.bbox.Add(vi->P());
	}
	dual_samples.vn = dual_samples.vert.size();
}



void DataMgr::loadPlyToSample(QString fileName)
{
	clearCMesh(samples);
	curr_file_name = fileName;

	int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;
	mask += tri::io::Mask::IOM_VERTCOLOR;
	mask += tri::io::Mask::IOM_BITPOLYGONAL;

	int err = tri::io::Importer<CMesh>::Open(samples, curr_file_name.toStdString().data(), mask);  
	if(err) 
	{
		cout << "Failed reading mesh: " << err << "\n";
		return;
	}  

	CMesh::VertexIterator vi;
	int idx = 0;
	for(vi = samples.vert.begin(); vi != samples.vert.end(); ++vi)
	{
		vi->bIsOriginal = false;
		vi->m_index = idx++;
		samples.bbox.Add(vi->P());
	}
	samples.vn = samples.vert.size();
}

void DataMgr::loadXYZN(QString fileName)
{
  clearCMesh(samples);
  ifstream infile;
  infile.open(fileName.toStdString().c_str());

  int i = 0;
  while(!infile.eof())
  {
    CVertex v;
    float temp = 0.;
    for (int j=0; j<3; j++)
    {
      infile >> temp;
      v.P()[j] = temp;
    }


//      for (int j=0; j<3; j++) {
//        infile >> v.N()[j];
//      }

// 		for (int j = 0; j < 3; j++) {
// 			infile >> temp;
// 			//infile >> v.C()[j];
// 		}

    v.m_index = i++;

    samples.vert.push_back(v);
    samples.bbox.Add(v.P());
  }

 // mesh.vert.erase(mesh.vert.end()-1);
  //samples.vert.pop_back();
	samples.vert.erase(samples.vert.end() - 1);
  samples.vn = samples.vert.size();

  infile.close();
}


void DataMgr::tryFixPly(QString fileName)
{
  clearCMesh(samples);
  ifstream infile;
  infile.open(fileName.toStdString().c_str());

  std::string str;
  int num;
  while (true)
  {
    infile >> str;
    if (str == std::string("vertex"))
    {
      infile >> num;
    }
    if (str == std::string("end_header"))
    {
      break;
    }
  }

  int i = 0;
  double x_value = 1.0;
  Point3f zero_pt(0., 0., 0.);
  for (; i < num; i++)
  {
    CVertex v;
    double temp = 0.;
    for (int j = 0; j < 3; j++)
    {
      infile >> temp;
      v.P()[j] = temp;

      if (j < 1)
      {
        x_value = temp;
      }
    }


    for (int j=0; j < 3; j++) {
      infile >> temp;
      v.N()[j] = temp;
    }

   for (int j = 0; j < 4; j++) {
   	infile >> temp;
   	//infile >> v.C()[j];
   }

    v.m_index = i;

    if (x_value > 1e100
       ||  x_value < -1e100
       || (x_value > 0 && x_value < 1e-100)
       || (x_value < 0 && x_value > -1e-100)
       || GlobalFun::computeEulerDistSquare(v.P(), zero_pt) < 1e-100 
       )
    {
      continue;
    }

    samples.vert.push_back(v);
    samples.bbox.Add(v.P());

    if (infile.eof())
    {
      break;
    }
  }


  samples.vn = samples.vert.size();

  infile.close();
}

void DataMgr::loadOFF(QString fileName)
{
  clearCMesh(original);
  ifstream infile;
  infile.open(fileName.toStdString().c_str());

  std::string tem_str;
  int tem_int;

  infile >> tem_str;
  int points_num, faces_num;
  infile >> points_num >> faces_num;
  infile >> tem_int;

  cout << "off points number: " << points_num << endl;
  for (int i = 0; i < points_num; i++)
  {
    CVertex v;
    float temp = 0.;
    for (int j = 0; j < 3; j++)
    {
      infile >> temp;
      v.P()[j] = temp;
    }

    v.m_index = i;

    original.vert.push_back(v);
    original.bbox.Add(v.P());
  }

  original.vn = original.vert.size();

  set<int> valid_id;
  vector< vector<Point3f>> valid_nomals;
  vector<Point3f> temp_nomals;
  valid_nomals.assign(original.vert.size(), temp_nomals);

  for (int i = 0; i < faces_num; i++)
  {
    int temp = 0.;
    vector<Point3f> three_points;
    vector<int> three_indexes;
    for (int j = 0; j < 3; j++)
    {
      infile >> temp;
      valid_id.insert(temp);

      three_points.push_back(samples.vert[temp].P());
      three_indexes.push_back(temp);
    }

    Point3f normal = (three_points[0] - three_points[1]).Normalize() ^ (three_points[1] - three_points[2]).Normalize();

    for (int j = 0; j < 3; j++)
    {
      int index = three_indexes[j];
      valid_nomals[index].push_back(normal.Normalize());
    }
  }

  int id = 0;
  set<int>::iterator iter;
  vector<CVertex> valid_samples;
  for (iter = valid_id.begin(); iter != valid_id.end(); ++iter)
  {
    int index = *iter;
    CVertex v = samples.vert[index];
    v.m_index = id++;

    vector<Point3f> normals = valid_nomals[index];

    if (normals.empty())
    {
      cout << "wrong file" << endl;
      system("Pause");
    }

    Point3f avg_normal = Point3f(0, 0, 0);
    for (int j = 0; j < normals.size(); j++)
    {
      avg_normal += normals[j];
    }
    avg_normal = (avg_normal / normals.size()).Normalize();
    v.N() = avg_normal;
    v.recompute_m_render();
    v.bIsOriginal = true;
    valid_samples.push_back(v);
  }

  original.vert.clear();
  original.bbox.SetNull();
  for (int i = 0; i < valid_samples.size(); i++)
  {
    valid_samples[i].bIsOriginal = true;
    original.vert.push_back(valid_samples[i]);
    original.bbox.Add(valid_samples[i]);
  }
  original.vn = original.vert.size();

  infile.close();

}




void DataMgr::loadImage(QString fileName)
{

	//image = cv::imread(fileName.toStdString().data());

	////cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
	////cv::imshow("image", image);

	//clearCMesh(samples);
	//clearCMesh(original);
	//int cnt = 0;
	//for (int i = 0; i < image.rows; i++)
	//{
	//	for (int j = 0; j < image.cols; j++)
	//	{
	//		cv::Vec3b intensity = image.at<cv::Vec3b>(i, j);
	//		Point3f p;
	//		Color4b c;
	//		c.Z() = 1;
	//		p.X() = c.X() = intensity.val[0];
	//		p.Y() = c.Y() = intensity.val[1];
	//		p.Z() = c.Z() = intensity.val[2];
	//		CVertex new_v;
	//		new_v.P() = p;
	//		new_v.C() = c;
	//		new_v.m_index = cnt++;

	//		samples.vert.push_back(new_v);
	//		samples.bbox.Add(p);

	//		new_v.bIsOriginal = true;
	//		original.vert.push_back(new_v);
	//		original.bbox.Add(p);
	//	}
	//}
	//samples.vn = samples.vert.size();
	//original.vn = samples.vn;

	//cv::waitKey();

}

CMesh*  DataMgr::getCurrentTargetSamples()
{
	if (&target_samples == NULL)
	{
		//cout << "DataMgr::getCurrentSamples samples = NULL!!" <<endl;
		return NULL;
	}

	if (target_samples.vert.empty())
	{
		//cout << "DataMgr::getCurrentSamples samples.vert.empty()!!" <<endl;
		//return NULL;
		return &target_samples;
	}

	return &target_samples;
}

CMesh*  DataMgr::getCurrentTargetDualSamples()
{
	if (&target_dual_samples == NULL)
	{
		//cout << "DataMgr::getCurrentSamples samples = NULL!!" <<endl;
		return NULL;
	}

	if (target_dual_samples.vert.empty())
	{
		//cout << "DataMgr::getCurrentSamples samples.vert.empty()!!" <<endl;
		//return NULL;
		return &target_dual_samples;
	}

	return &target_dual_samples;
}

CMesh* DataMgr::getCurrentSamples()
{
  if(&samples == NULL)
  {
    //cout << "DataMgr::getCurrentSamples samples = NULL!!" <<endl;
    return NULL;
  }

	if(samples.vert.empty())
	{
		//cout << "DataMgr::getCurrentSamples samples.vert.empty()!!" <<endl;
		//return NULL;
    return & samples;
	}

	return & samples;
}

CMesh* DataMgr::getCurrentDualSamples()
{
  if(&dual_samples == NULL)
  {
    //cout << "DataMgr::getCurrentSamples samples = NULL!!" <<endl;
    return NULL;
  }

  if(dual_samples.vert.empty())
  {
    //cout << "DataMgr::getCurrentSamples samples.vert.empty()!!" <<endl;
    //return NULL;
    return & dual_samples;
  }

  return & dual_samples;
}

CMesh* DataMgr::getCurrentSkelPoints()
{
	if (&skel_points == NULL)
	{
		return NULL;
	}

	if (skel_points.vert.empty())
	{
		return &skel_points;
	}

	return &skel_points;
}

CMesh* DataMgr::getCurrentOriginal()
{
  if(&original == NULL)
  {
    //cout << "DataMgr::getCurrentOriginal() samples = NULL!!" <<endl;
    return NULL;
  }

	if(original.vert.empty())
	{
		//cout << "DataMgr::getCurrentOriginal() original.vert.empty()!!" <<endl;
		return & original;
	}

	return & original;
}

CMesh* DataMgr::getCurrentEllipsoid()
{
	if (&ellipsoid_mesh == NULL)
	{
		//cout << "DataMgr::getCurrentOriginal() samples = NULL!!" <<endl;
		return NULL;
	}

	if (ellipsoid_mesh.vert.empty())
	{
		//cout << "DataMgr::getCurrentOriginal() original.vert.empty()!!" <<endl;
		return &ellipsoid_mesh;
	}

	return &ellipsoid_mesh;
}

Skeleton* DataMgr::getCurrentSkeleton()
{
	return & skeleton;
}


void DataMgr::recomputeBox()
{
	samples.bbox.SetNull();
	original.bbox.SetNull();
	dual_samples.bbox.SetNull();
	skel_points.bbox.SetNull();

	CMesh::VertexIterator vi;
	for(vi = samples.vert.begin(); vi != samples.vert.end(); ++vi) 
	{
		if (vi->is_skel_ignore)
		{
			continue;
		}
		samples.bbox.Add(vi->P());
	}

	for(vi = original.vert.begin(); vi != original.vert.end(); ++vi) 
	{
		original.bbox.Add(vi->P());
	}

	for (vi = dual_samples.vert.begin(); vi != dual_samples.vert.end(); ++vi)
	{
		dual_samples.bbox.Add(vi->P());
	}

	for (vi = skel_points.vert.begin(); vi != skel_points.vert.end(); ++vi)
	{
		skel_points.bbox.Add(vi->P());
	}
}

double DataMgr::getInitRadiuse()
{
	double init_para = para->getDouble("Init Radius Para");
	if (!isOriginalEmpty())
	{
		Box3f box = original.bbox;
		if ( abs(box.min.X() - box.max.X()) < 1e-5 ||   
			abs(box.min.Y() - box.max.Y()) < 1e-5 ||   
			abs(box.min.Z() - box.max.Z()) < 1e-5 )
		{
			double diagonal_length = sqrt((box.min - box.max).SquaredNorm());
			double original_size = sqrt(double(original.vn));
			init_radius = 2 * init_para * diagonal_length / original_size;
		}
		else
		{
			double diagonal_length = sqrt((box.min - box.max).SquaredNorm());
			double original_size = pow(double(original.vn), 0.333);
			init_radius = init_para * diagonal_length / original_size;
		}
	}

	double current_radius = global_paraMgr.wLop.getDouble("CGrid Radius");
	if (current_radius < 0 || global_paraMgr.glarea.getBool("Show Skeleton"))
	{
		global_paraMgr.setGlobalParameter("CGrid Radius", DoubleValue(init_radius));
	}

  global_paraMgr.setGlobalParameter("Initial Radius", DoubleValue(init_radius));

	return init_radius;
}


// void DataMgr::downSamplesByNum(bool use_random_downsample)
// {
// 	if (isOriginalEmpty() && !isSamplesEmpty())
// 	{
// 		subSamples();
// 		return;
// 	}
// 
// 	if (isOriginalEmpty())
// 	{
// 		return;
// 	}
// 
// 	double radius = para->getDouble("CGrid Radius");
// 	double sigma = global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma");
// 	//GlobalFun::computeBilateralConfidence(&original, radius, sigma);
// 	cout << "GlobalFun::computeBilateralConfidence" << endl;
// 
// 	int want_sample_num = para->getDouble("Down Sample Num");
// 
// 	if (want_sample_num > original.vn)
// 	{
// 		want_sample_num = original.vn;
// 	}
// 
// 	clearCMesh(samples);
// 	samples.vn = want_sample_num;
// 
//  	vector<int> nCard = GlobalFun::GetRandomCards(original.vert.size());
// 
// 	int inserted_number = 0;
// 	int i = 0;
// //  	while (want_sample_num > inserted_number)
// //  	{
// //  		if (i >= nCard.size())
// //  		{
// //  			break;
// //  		}
// //  
// //  		int index = nCard[i++];
// //  
// //  		CVertex& v = original.vert[index];
// //  		double probability = 1 - v.eigen_confidence;
// // 		//double probability = v.eigen_confidence;
// // 
// //  		double r = (rand() % 1000) * 0.001;
// //  
// //  		if (r < probability)
// //  		{
// //  			samples.vert.push_back(v);
// //  			samples.bbox.Add(v.P());
// //  
// //  			inserted_number++;
// //  		}
// //  	}
// 
//    	for (int i = 0; i < samples.vn; i++)
//    	{
//    		int index = nCard[i]; //could be not random!
//    
//    		if (!use_random_downsample)
//    		{
//    			index = i;
//    		}
//    
//    		CVertex& v = original.vert[index];
//    		v.dual_index = i;
//    		samples.vert.push_back(v);
//    		samples.bbox.Add(v.P());
//    	}
//  
//  	CMesh::VertexIterator vi;
//  	for (vi = samples.vert.begin(); vi != samples.vert.end(); ++vi)
//  	{
//  		vi->bIsOriginal = false;
//  	}
//  
//  	dual_samples.vert.clear();
//  	for (int i = 0; i < samples.vert.size(); i++)
//  	{
//  		CVertex v = samples.vert[i];
//  		v.is_dual_sample = true;
//  		v.dual_index = i;
//  		dual_samples.vert.push_back(v);
//  	}
//  	dual_samples.bbox = samples.bbox;
//  	dual_samples.vn = samples.vn;
// 
// 	cout << "compute density confidence" << endl;
// 	GlobalFun::computeBallNeighbors(&original, NULL,
// 		para->getDouble("CGrid Radius"), original.bbox);
// 
// 	getInitRadiuse();
// }


void DataMgr::downSamplesByNum(bool use_random_downsample)
{
	if (isOriginalEmpty() && !isSamplesEmpty())
	{
		subSamples();
		return;
	}

	if (isOriginalEmpty())
	{
		return;
	}

	int want_sample_num = para->getDouble("Down Sample Num");

	if (want_sample_num > original.vn)
	{
		want_sample_num = original.vn;
	}

	clearCMesh(samples);
	samples.vn = want_sample_num;

	vector<int> nCard = GlobalFun::GetRandomCards(original.vert.size());
	for(int i = 0; i < samples.vn; i++) 
	{
		int index = nCard[i]; //could be not random!

    if (!use_random_downsample)
    {
      index = i;
    }

		CVertex& v = original.vert[index];
    v.dual_index = i;
		samples.vert.push_back(v);
		samples.bbox.Add(v.P());
	}

	CMesh::VertexIterator vi;
	for(vi = samples.vert.begin(); vi != samples.vert.end(); ++vi)
	{
		vi->bIsOriginal = false;
	}

  dual_samples.vert.clear();
  for (int i = 0; i < samples.vert.size(); i++)
  {
    CVertex v = samples.vert[i];
    v.is_dual_sample = true;
    v.dual_index = i;
    dual_samples.vert.push_back(v);
  }
  dual_samples.bbox = samples.bbox;
  dual_samples.vn = samples.vn;

  getInitRadiuse();
}


void DataMgr::subSamples()
{
	clearCMesh(original);

	CMesh::VertexIterator vi;
	original.vn = samples.vert.size();
	original.bbox.SetNull();
	for(vi = samples.vert.begin(); vi != samples.vert.end(); ++vi)
	{
		CVertex v = (*vi);
		v.bIsOriginal = true;
		original.vert.push_back(v);
		original.bbox.Add(v.P());
	}



	downSamplesByNum();
  getInitRadiuse();
}


void DataMgr::savePly(QString fileName, CMesh& mesh)
{
	int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;
	//mask += tri::io::Mask::IOM_VERTCOLOR;
	mask += tri::io::Mask::IOM_BITPOLYGONAL;

	if (fileName.endsWith("ply"))
		tri::io::ExporterPLY<CMesh>::Save(mesh, fileName.toStdString().data(), mask, false);
}

//void DataMgr::normalizeROSA_Mesh(CMesh& mesh)
//{
//	if (mesh.vert.empty())
//	{
//		return;
//	}
//	Box3f box = mesh.bbox;
//	mesh.bbox.SetNull();
//	float max_x = abs((box.min - box.max).X());
//	float max_y = abs((box.min - box.max).Y());
//	float max_z = abs((box.min - box.max).Z());
//	float max_length = max_x > max_y ? max_x : max_y;
//	max_length = max_length > max_z ? max_length : max_z;
//
//	for(int i = 0; i < mesh.vert.size(); i++)
//	{
//		Point3f& p = mesh.vert[i].P();
//
//		p -= box.min;
//		p /= max_length;
//
//		p = (p - Point3f(0.5, .5, .5));
//
//		mesh.vert[i].N().Normalize(); 
//		mesh.bbox.Add(p);
//	}
//}



void DataMgr::normalizeROSA_MeshForOriginal(CMesh& mesh, Point3f mid_point)
{
	if (mesh.vert.empty()) return;

	Box3f box = mesh.bbox;

	mesh.bbox.SetNull();
	float max_x = abs((box.min - box.max).X());
	float max_y = abs((box.min - box.max).Y());
	float max_z = abs((box.min - box.max).Z());
	float max_length = max_x > max_y ? max_x : max_y;
	max_length = max_length > max_z ? max_length : max_z;

	Box3f box_temp;
	for (int i = 0; i < mesh.vert.size(); i++)
	{
		Point3f& p = mesh.vert[i].P();

		p /= max_length;

		mesh.vert[i].N().Normalize();
		box_temp.Add(p);
	}

// 	Point3f mid_point = (box_temp.min + box_temp.max) / 2.0;
// 	mid_point = (box.min + box.max) / 2.0;


	mesh.bbox.SetNull();
	for (int i = 0; i < mesh.vert.size(); i++)
	{
		Point3f& p = mesh.vert[i].P();
		p -= mid_point;
		p *= 2.0;
		mesh.bbox.Add(p);
	}
}


void DataMgr::normalizeROSA_UsingKnownCondition(CMesh& mesh, Point3f center, double length)
{
  mesh.bbox.SetNull();
  for (int i = 0; i < mesh.vert.size(); i++)
  {
    Point3f& p = mesh.vert[i].P();

    p /= length;
    p -= center;
    p *= 2.0;

    mesh.vert[i].N().Normalize();
    mesh.bbox.Add(p);
  }


}


void DataMgr::saveNomalization(QString fileName)
{
  ofstream outfile;
  outfile.open(fileName.toStdString().c_str());

  Point3f center_point = global_paraMgr.data.getPoint3f("Normalization Center Point");
  double length = global_paraMgr.data.getDouble("Normalization Length");

  outfile << length << "  " << center_point.X() << "  " << center_point.Y() << "  " << center_point.Z() << endl;
  outfile.close();
}



void DataMgr::loadNomalization(QString fileName)
{
  ifstream infile;
  infile.open(fileName.toStdString().c_str());

  Point3f center_point;
  double length;

  infile >> length;
  infile >> center_point.X() >> center_point.Y() >> center_point.Z();

  normalizeROSA_UsingKnownCondition(original, center_point, length);
  normalizeROSA_UsingKnownCondition(samples, center_point, length);
  normalizeROSA_UsingKnownCondition(dual_samples, center_point, length);
  normalizeROSA_UsingKnownCondition(skel_points, center_point, length);

  getInitRadiuse();
}



Point3f DataMgr::normalizeROSA_Mesh(CMesh& mesh)
{
	Box3f box = mesh.bbox;

	mesh.bbox.SetNull();
	float max_x = abs((box.min - box.max).X());
	float max_y = abs((box.min - box.max).Y());
	float max_z = abs((box.min - box.max).Z());
	float max_length = max_x > max_y ? max_x : max_y;
	max_length = max_length > max_z ? max_length : max_z;

	Box3f box_temp;
	for (int i = 0; i < mesh.vert.size(); i++)
	{
		Point3f& p = mesh.vert[i].P();

		p /= max_length;

		mesh.vert[i].N().Normalize();
		box_temp.Add(p);
	}

	Point3f mid_point = (box_temp.min + box_temp.max) / 2.0;

  para->setValue("Normalization Center Point", Point3fValue(mid_point));
  para->setValue("Normalization Length", DoubleValue(max_length));


	mesh.bbox.SetNull();
	for (int i = 0; i < mesh.vert.size(); i++)
	{
		Point3f& p = mesh.vert[i].P();
		p -= mid_point;
		p *= 2.0;
		mesh.bbox.Add(p);
	}

	return mid_point;
}

//Box3f DataMgr::normalizeAllMesh()
//{
//	Box3f box;
//	if (!isSamplesEmpty())
//	{
//		for (int i = 0; i < samples.vert.size(); i++)
//		{
//			box.Add(samples.vert[i].P());
//		}
//	}
//	if (!isOriginalEmpty())
//	{
//		for (int i = 0; i < original.vert.size(); i++)
//		{
//			box.Add(original.vert[i].P());
//		}
//		original.bbox =box;
//	}
//
//	samples.bbox = box;
//	dual_samples.bbox = box;
//
//	normalizeROSA_Mesh(samples);
//  normalizeROSA_Mesh(dual_samples);
//	normalizeROSA_Mesh(original);
//
//	recomputeBox();
//	getInitRadiuse();
//
//	return samples.bbox;
//}

Box3f DataMgr::normalizeAllMesh()
{
	Box3f box;
	if (!isSamplesEmpty())
	{
		for (int i = 0; i < samples.vert.size(); i++)
		{
			box.Add(samples.vert[i].P());
		}
	}
	if (!dual_samples.vert.empty())
	{
		for (int i = 0; i < dual_samples.vert.size(); i++)
		{
			box.Add(dual_samples.vert[i].P());
		}
	}
	if (!isOriginalEmpty())
	{
		for (int i = 0; i < original.vert.size(); i++)
		{
			box.Add(original.vert[i].P());
		}
	}

	samples.bbox = box;
	dual_samples.bbox = box;
	original.bbox = box;
  skel_points.bbox = box;

	Point3f mid = normalizeROSA_Mesh(samples);
  //normalizeROSA_Mesh(dual_samples);
	normalizeROSA_MeshForOriginal(dual_samples, mid);
	normalizeROSA_MeshForOriginal(skel_points, mid);
	normalizeROSA_MeshForOriginal(original, mid);

	recomputeBox();
	getInitRadiuse();

	return samples.bbox;
}


void DataMgr::eraseRemovedSamples()
{
	int cnt = 0;
	vector<CVertex> temp_mesh;
	for (int i = 0; i < samples.vert.size(); i++)
	{
		CVertex& v = samples.vert[i];
		if (!v.is_skel_ignore)
		{
			temp_mesh.push_back(v);
		}
	}

	samples.vert.clear();
	samples.vn = temp_mesh.size();
	for (int i = 0; i < temp_mesh.size(); i++)
	{
		temp_mesh[i].m_index = i;
		samples.vert.push_back(temp_mesh[i]);
	}

}

void DataMgr::clearData()
{
	clearCMesh(original);
	clearCMesh(samples);
  clearCMesh(dual_samples);
	clearCMesh(skel_points);

	clearCMesh(target_samples);
	clearCMesh(target_dual_samples);

	skeleton.clear();
}

void DataMgr::recomputeQuad()
{
	for (int i = 0; i < samples.vert.size(); i++)
	{
		samples.vert[i].recompute_m_render();
	}

	for (int i = 0; i < dual_samples.vert.size(); i++)
  {
    dual_samples.vert[i].recompute_m_render();
  }

  for (int i = 0; i < original.vert.size(); i++)
  {
    original.vert[i].recompute_m_render();
  }
}


void DataMgr::saveTargetSkeletonAsSkel(QString fileName)
{
	ofstream outfile;
	outfile.open(fileName.toStdString().c_str());

	ostringstream strStream;

	strStream << "ON " << original.vert.size() << endl;
	for (int i = 0; i < original.vert.size(); i++)
	{
		CVertex& v = original.vert[i];
		strStream << v.P()[0] << "	" << v.P()[1] << "	" << v.P()[2] << "	";
		strStream << v.N()[0] << "	" << v.N()[1] << "	" << v.N()[2] << "	" << endl;
	}
	strStream << endl;

	strStream << "SN " << target_samples.vert.size() << endl;
	for (int i = 0; i < target_samples.vert.size(); i++)
	{
		CVertex& v = target_samples.vert[i];
		strStream << v.P()[0] << "	" << v.P()[1] << "	" << v.P()[2] << "	";
		strStream << v.N()[0] << "	" << v.N()[1] << "	" << v.N()[2] << "	" << endl;
	}
	strStream << endl;

	strStream << "CN " << skeleton.branches.size() << endl;
	for (int i = 0; i < skeleton.branches.size(); i++)
	{
		Branch& branch = skeleton.branches[i];
		strStream << "CNN " << branch.curve.size() << endl;
		for (int j = 0; j < branch.curve.size(); j++)
		{
			strStream << branch.curve[j][0] << "	" << branch.curve[j][1] << "	" << branch.curve[j][2] << "	" << endl;
		}
	}
	strStream << endl;

	strStream << "EN " << 0 << endl;
	strStream << endl;


	strStream << "BN " << 0 << endl;
	strStream << endl;


	strStream << "S_onedge " << target_samples.vert.size() << endl;
	for (int i = 0; i < target_samples.vert.size(); i++)
	{
		CVertex& v = target_samples.vert[i];
		strStream << v.is_fixed_sample << "	";
	}
	strStream << endl;

	strStream << "GroupID " << target_samples.vert.size() << endl;
	for (int i = 0; i < target_samples.vert.size(); i++)
	{
		int id = target_samples.vert[i].m_index;//group_id no use now
		strStream << id << "	";
	}
	strStream << endl;

	//strStream << "SkelRadius " << 0 << endl;
	//strStream << endl;

	strStream << "SkelRadius " << skeleton.size << endl;
	for (int i = 0; i < skeleton.branches.size(); i++)
	{
		for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
		{
			double skel_radius = skeleton.branches[i].curve[j].skel_radius;
			strStream << skel_radius << "	";
		}
	}
	strStream << endl;

	strStream << "Confidence_Sigma	" << target_samples.vert.size() << endl;
	for (int i = 0; i < target_samples.vert.size(); i++)
	{
		double sigma = target_samples.vert[i].eigen_confidence;
		strStream << sigma << "	";
	}
	strStream << endl;

	strStream << "SkelRadius2 " << 0 << endl;
	strStream << endl;

	strStream << "Alpha " << 0 << endl;
	strStream << endl;

	strStream << "Sample_isVirtual " << target_samples.vert.size() << endl;
	for (int i = 0; i < target_samples.vert.size(); i++)
	{
		CVertex& v = target_samples.vert[i];
		strStream << v.is_skel_virtual << "	";
	}
	strStream << endl;

	strStream << "Sample_isBranch " << target_samples.vert.size() << endl;
	for (int i = 0; i < target_samples.vert.size(); i++)
	{
		CVertex& v = target_samples.vert[i];
		strStream << v.is_skel_branch << "	";
	}
	strStream << endl;

	strStream << "Sample_radius " << target_samples.vert.size() << endl;
	for (int i = 0; i < target_samples.vert.size(); i++)
	{
		CVertex& v = target_samples.vert[i];
		strStream << 0 << "	";
	}
	strStream << endl;

	strStream << "Skel_isVirtual " << skeleton.size << endl;
	for (int i = 0; i < skeleton.branches.size(); i++)
	{
		for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
		{
			bool is_virtual = skeleton.branches[i].curve[j].is_skel_virtual;
			strStream << is_virtual << "	";
		}
	}
	strStream << endl;


	strStream << "Corresponding_sample_index " << skeleton.size << endl;
	for (int i = 0; i < skeleton.branches.size(); i++)
	{
		for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
		{
			int index = skeleton.branches[i].curve[j].m_index;
			strStream << index << "	";
		}
	}
	strStream << endl;

	strStream << "DSN " << target_dual_samples.vert.size() << endl;
	for (int i = 0; i < target_dual_samples.vert.size(); i++)
	{
		CVertex& v = target_dual_samples.vert[i];
		strStream << v.P()[0] << "	" << v.P()[1] << "	" << v.P()[2] << "	";
		strStream << v.N()[0] << "	" << v.N()[1] << "	" << v.N()[2] << "	" << endl;
	}
	strStream << endl;

	strStream << "dual_corresponding_index " << target_samples.vert.size() << endl;
	for (int i = 0; i < target_dual_samples.vert.size(); i++)
	{
		CVertex& v = target_dual_samples.vert[i];
		strStream << target_dual_samples.vert[i].dual_index << endl;
	}
	strStream << endl;


	outfile.write(strStream.str().c_str(), strStream.str().size());
	outfile.close();

}




void DataMgr::saveSkeletonAsSkel(QString fileName)
{
	ofstream outfile;
	outfile.open(fileName.toStdString().c_str());

	ostringstream strStream; 

	strStream << "ON " << original.vert.size() << endl;
	for(int i = 0; i < original.vert.size(); i++)
	{
		CVertex& v = original.vert[i];
		strStream << v.P()[0] << "	" << v.P()[1] << "	" << v.P()[2] << "	";
		strStream << v.N()[0] << "	" << v.N()[1] << "	" << v.N()[2] << "	" << endl;
	}
	strStream << endl;

	strStream << "SN " << samples.vert.size() << endl;
	for(int i = 0; i < samples.vert.size(); i++)
	{
		CVertex& v = samples.vert[i];
		strStream << v.P()[0] << "	" << v.P()[1] << "	" << v.P()[2] << "	";
		strStream << v.N()[0] << "	" << v.N()[1] << "	" << v.N()[2] << "	" << endl;
	}
	strStream << endl;

	strStream << "CN " << skeleton.branches.size() << endl;
	for (int i = 0; i < skeleton.branches.size(); i++)
	{
		Branch& branch = skeleton.branches[i];
		strStream << "CNN " << branch.curve.size() << endl;
		for (int j = 0; j < branch.curve.size(); j++)
		{
			strStream << branch.curve[j][0] << "	" << branch.curve[j][1] << "	" << branch.curve[j][2] << "	" << endl;
		}
	}
	strStream << endl;

	strStream << "EN " << 0 << endl;
	strStream << endl;


	strStream << "BN " << 0 << endl;
	strStream << endl;


	strStream << "S_onedge " << samples.vert.size() << endl;
	for(int i = 0; i < samples.vert.size(); i++)
	{
		CVertex& v = samples.vert[i];
		strStream << v.is_fixed_sample << "	"; 
	}
	strStream << endl;

	strStream << "GroupID " << samples.vert.size() << endl;
	for(int i = 0; i < samples.vert.size(); i++)
	{
		int id = samples.vert[i].m_index;//group_id no use now
		strStream << id << "	"; 
	}
	strStream << endl;

	//strStream << "SkelRadius " << 0 << endl;
	//strStream << endl;
  
  strStream << "SkelRadius " << skeleton.size << endl;
  for (int i = 0; i < skeleton.branches.size(); i++)
  {
    for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
    {
      double skel_radius = skeleton.branches[i].curve[j].skel_radius;
      strStream << skel_radius << "	"; 
    }
  }
  strStream << endl;

	strStream << "Confidence_Sigma	" << samples.vert.size() << endl;
	for(int i = 0; i < samples.vert.size(); i++)
	{
		double sigma = samples.vert[i].eigen_confidence;
		strStream << sigma << "	"; 
	}
	strStream << endl;

	strStream << "SkelRadius2 " << 0 << endl;
	strStream << endl;

	strStream << "Alpha " << 0 << endl;
	strStream << endl;

	strStream << "Sample_isVirtual " << samples.vert.size() << endl;
	for(int i = 0; i < samples.vert.size(); i++)
	{
		CVertex& v = samples.vert[i];
		strStream << v.is_skel_virtual << "	"; 
	}
	strStream << endl;

	strStream << "Sample_isBranch " << samples.vert.size() << endl;
	for(int i = 0; i < samples.vert.size(); i++)
	{
		CVertex& v = samples.vert[i];
		strStream << v.is_skel_branch << "	"; 
	}
	strStream << endl;

  strStream << "Sample_radius " << samples.vert.size() << endl;
  for(int i = 0; i < samples.vert.size(); i++)
  {
    CVertex& v = samples.vert[i];
    strStream << 0 << "	"; 
  }
  strStream << endl;

	strStream << "Skel_isVirtual " << skeleton.size << endl;
	for (int i = 0; i < skeleton.branches.size(); i++)
	{
		for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
		{
			bool is_virtual = skeleton.branches[i].curve[j].is_skel_virtual;
			strStream << is_virtual << "	"; 
		}
	}
	strStream << endl;


  strStream << "Corresponding_sample_index " << skeleton.size << endl;
  for (int i = 0; i < skeleton.branches.size(); i++)
  {
    for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
    {
      int index = skeleton.branches[i].curve[j].m_index;
      strStream << index << "	"; 
    }
  }
  strStream << endl;

  strStream << "DSN " << dual_samples.vert.size() << endl;
  for(int i = 0; i < dual_samples.vert.size(); i++)
  {
    CVertex& v = dual_samples.vert[i];
    strStream << v.P()[0] << "	" << v.P()[1] << "	" << v.P()[2] << "	";
    strStream << v.N()[0] << "	" << v.N()[1] << "	" << v.N()[2] << "	" << endl;
  }
  strStream << endl;

  strStream << "dual_corresponding_index " << samples.vert.size() << endl;
  for (int i = 0; i < samples.vert.size(); i++)
  {
    CVertex& v = samples.vert[i];
    strStream << samples.vert[i].dual_index << endl;
  }
  strStream << endl;

	strStream << "skel_radius " << samples.vert.size() << endl;
	for (int i = 0; i < samples.vert.size(); i++)
	{
		CVertex& v = samples.vert[i];
		strStream << samples.vert[i].skel_radius << endl;

	}
	strStream << endl;

	strStream << "SkelPN " << skel_points.vert.size() << endl;
	for (int i = 0; i < skel_points.vert.size(); i++)
	{
		CVertex& v = skel_points.vert[i];
		strStream << v.P()[0] << "	" << v.P()[1] << "	" << v.P()[2] << "	";
		strStream << v.N()[0] << "	" << v.N()[1] << "	" << v.N()[2] << "	" << endl;
	}
	strStream << endl;

	outfile.write( strStream.str().c_str(), strStream.str().size() ); 
	outfile.close();
}


void DataMgr::loadTargetSkeletonFromSkel(QString fileName)
{
	clearCMesh(target_samples);
	clearCMesh(original);
	clearCMesh(target_dual_samples);

	ifstream infile;
	infile.open(fileName.toStdString().c_str());

	stringstream sem;
	sem << infile.rdbuf();

	string str;
	int num;
	int num2;

	sem >> str;
	if (str == "ON")
	{
		sem >> num;
		bool is_same_original = false;
		if (num == original.vn)
		{
			is_same_original = true;
		}
		if (is_same_original)
		{
			double temp;
			for (int i = 0; i < num * 6; i++)
			{
				sem >> temp;
			}
		}
		else
		{
			for (int i = 0; i < num; i++)
			{
				CVertex v;
				v.bIsOriginal = true;
				v.m_index = i;
				sem >> v.P()[0] >> v.P()[1] >> v.P()[2];
				sem >> v.N()[0] >> v.N()[1] >> v.N()[2];
				original.vert.push_back(v);
				original.bbox.Add(v.P());
			}
			original.vn = original.vert.size();
		}
	}

	sem >> str;
	if (str == "SN")
	{
		sem >> num;
		for (int i = 0; i < num; i++)
		{
			CVertex v;
			v.bIsOriginal = false;
			v.m_index = i;
			sem >> v.P()[0] >> v.P()[1] >> v.P()[2];
			sem >> v.N()[0] >> v.N()[1] >> v.N()[2];
			target_samples.vert.push_back(v);
			target_samples.bbox.Add(v.P());
		}
		target_samples.vn = target_samples.vert.size();
	}


	sem >> str;
	if (str == "CN")
	{
		sem >> num;
		for (int i = 0; i < num; i++)
		{
			Branch branch;
			sem >> str;
			sem >> num2;
			for (int j = 0; j < num2; j++)
			{
				Point3f p;
				CVertex v;
				sem >> p[0] >> p[1] >> p[2];
				v.P() = p;
				branch.curve.push_back(v);
			}
			skeleton.branches.push_back(branch);
		}
	}

	sem >> str;
	if (str == "EN")
	{
		sem >> num;
		for (int i = 0; i < num; i++)
		{
			int a, b;
			sem >> a >> b;
		}
	}

	sem >> str;
	if (str == "BN")
	{
		sem >> num;
		for (int i = 0; i < num; i++)
		{
			sem >> str;
			sem >> num2;

			for (int j = 0; j < num2; j++)
			{
				int id;
				sem >> id;

			}
		}

		if (!sem.eof())
		{
			sem >> str;
			if (str == "S_onedge")
			{
				sem >> num;
				for (int i = 0; i < num; i++)
				{
					bool b;
					sem >> b;
					target_samples.vert[i].is_fixed_sample = b;
				}
			}
		}

		sem >> str;
		if (str == "GroupID")
		{
			sem >> num;
			for (int i = 0; i < num; i++)
			{
				int id;
				sem >> id;

			}
		}
	}

	sem >> str;
	if (str == "SkelRadius")
	{
		sem >> num;

		if (num > 1)
		{
			double radius;
			for (int i = 0; i < skeleton.branches.size(); i++)
			{
				for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
				{
					sem >> radius;
					skeleton.branches[i].curve[j].skel_radius = radius;
				}
			}
		}

	}

	sem >> str;
	if (str == "Confidence_Sigma")
	{
		sem >> num;
		for (int i = 0; i < num; i++)
		{
			double sigma;
			sem >> sigma;
			target_samples.vert[i].eigen_confidence = sigma;
		}
	}

	sem >> str;
	if (str == "SkelRadius2")
	{
		sem >> num;

		if (num > 1)
		{
			double radius;
			for (int i = 0; i < skeleton.branches.size(); i++)
			{
				for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
				{
					sem >> radius;
					//skeleton.branches[i].curve[j].skel_radius = radius;
				}
			}
		}

	}

	sem >> str;
	if (str == "Alpha")
	{
		sem >> num;
		double Alpha;
		if (num > 1)
		{
			for (int i = 0; i < skeleton.branches.size(); i++)
			{
				for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
				{
					sem >> Alpha;
					//skeleton.curves[i][j].alpha = Alpha;
				}
			}
		}

	}

	if (!sem.eof())
	{
		sem >> str;
		if (str == "Sample_isVirtual")
		{
			sem >> num;
			for (int i = 0; i < num; i++)
			{
				bool b;
				sem >> b;
				target_samples.vert[i].is_skel_virtual = b;
			}
		}
	}

	if (!sem.eof())
	{
		sem >> str;
		if (str == "Sample_isBranch")
		{
			sem >> num;
			for (int i = 0; i < num; i++)
			{
				bool b;
				sem >> b;
				target_samples.vert[i].is_skel_branch = b;
			}
		}
	}

	sem >> str;
	if (str == "Sample_radius")
	{
		sem >> num;
		for (int i = 0; i < num; i++)
		{
			double temp;
			sem >> temp;
			//target_samples.vert[i].saved_radius = temp;
		}
	}

	sem >> str;
	if (str == "Skel_isVirtual")
	{
		cout << "load Skel_isVirtual" << endl;

		sem >> num;
		bool temp;
		for (int i = 0; i < skeleton.branches.size(); i++)
		{
			for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
			{
				sem >> temp;
				skeleton.branches[i].curve[j].is_skel_virtual = temp;
			}
		}
	}

	sem >> str;
	if (str == "Corresponding_sample_index")
	{
		cout << "load Corresponding_sample_index" << endl;

		sem >> num;
		int temp;
		for (int i = 0; i < skeleton.branches.size(); i++)
		{
			for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
			{
				sem >> temp;
				skeleton.branches[i].curve[j].m_index = temp;
			}
		}
	}

	sem >> str;
	if (str == "DSN")
	{
		cout << "load DSN" << endl;

		sem >> num;
		for (int i = 0; i < num; i++)
		{
			CVertex v;
			v.bIsOriginal = false;
			v.is_dual_sample = true;
			v.m_index = i;
			sem >> v.P()[0] >> v.P()[1] >> v.P()[2];
			sem >> v.N()[0] >> v.N()[1] >> v.N()[2];
// 			target_dual_samples.vert.push_back(v);
// 			target_dual_samples.bbox.Add(v.P());
		}
		//target_dual_samples.vn = target_dual_samples.vert.size();
	}

	if (!sem.eof())
	{
		sem >> str;
		if (str == "dual_corresponding_index")
		{
			cout << "load dual_corresponding_index" << endl;

			sem >> num;
			for (int i = 0; i < num; i++)
			{
				int index;
				sem >> index;
				target_dual_samples.vert[i].dual_index = index;
			}
		}
	}

	if (!sem.eof())
	{
		sem >> str;
		if (str == "skel_radius")
		{
			cout << "load skel_radius" << endl;
			sem >> num;
			for (int i = 0; i < num; i++)
			{
				int radius;
				sem >> radius;
				dual_samples.vert[i].skel_radius = radius;
				samples.vert[i].skel_radius = radius;

				if (i < 20)
				{
					cout << "skel_radius: " << radius;
				}
				//if (radius > 0)
				//{

				//}
			}
		}
	}

	

	for (int i = 0; i < target_samples.vert.size(); i++)
	{
		target_samples.vert[i].N().Normalize();
	}
}


void DataMgr::loadSkeletonFromSkel(QString fileName)
{
	clearCMesh(samples);
	clearCMesh(original);
  clearCMesh(dual_samples);
	clearCMesh(skel_points);

	skeleton.clear();

	ifstream infile;
	infile.open(fileName.toStdString().c_str());

	stringstream sem; 
	sem << infile.rdbuf(); 

	string str;
	int num;
	int num2;

	sem >> str;
	if (str == "ON")
	{
		sem >> num;
		bool is_same_original = false;
		if (num == original.vn)
		{
			is_same_original = true;
		}
		if (is_same_original)
		{
			double temp;
			for (int i = 0; i < num * 6; i++)
			{
				sem >> temp;
			}
		}
		else
		{
			for (int i = 0; i < num; i++)
			{
				CVertex v;
				v.bIsOriginal = true;
				v.m_index = i;
				sem >> v.P()[0] >> v.P()[1] >> v.P()[2];
				sem >> v.N()[0] >> v.N()[1] >> v.N()[2];
				original.vert.push_back(v);
				original.bbox.Add(v.P());
			}
			original.vn = original.vert.size();
		}
	}

	sem >> str;
	if (str == "SN")
	{
		sem >> num;
		for (int i = 0; i < num; i++)
		{
			CVertex v;
			v.bIsOriginal = false;
			v.m_index = i;
			v.dual_index = i;

			sem >> v.P()[0] >> v.P()[1] >> v.P()[2];
			sem >> v.N()[0] >> v.N()[1] >> v.N()[2];
			samples.vert.push_back(v);
			samples.bbox.Add(v.P());

		}
		samples.vn = samples.vert.size();
	}


	sem >> str;
	if (str == "CN")
	{
		sem >> num;
		for (int i = 0; i < num; i++)
		{
			Branch branch;
			sem >> str;
			sem >> num2;
			for(int j = 0; j < num2; j++)
			{
				Point3f p;
				CVertex v;
				sem >> p[0] >> p[1] >> p[2];
				v.P() = p;
				branch.curve.push_back(v);
			}
			skeleton.branches.push_back(branch);
		}
	}

	sem >> str;
	if (str == "EN")
	{
    sem >> num;
    for (int i = 0; i < num; i++)
    {
      int a, b;
      sem >> a >> b;
    }
	}

	sem >> str;
	if (str == "BN")
	{
    sem >> num;
    for (int i = 0; i < num; i++)
    {
      sem >> str;
      sem >> num2;

      for(int j = 0; j < num2; j++)
      {
        int id;
        sem >> id;

    }
	}

	if (!sem.eof())
	{
		sem >> str;
		if (str == "S_onedge")
		{
			sem >> num;
			for (int i = 0; i < num; i++)
			{
				bool b;
				sem >> b;
				samples.vert[i].is_fixed_sample = b;
			}
		}
	}

	sem >> str;
	if (str == "GroupID")
	{
		sem >> num;
		for (int i = 0; i < num; i++)
		{
			int id;
			sem >> id;

      }
		}
	}

	sem >> str;
	if (str == "SkelRadius")
	{
		sem >> num;

    if (num > 1)
    {
      double radius;
      for (int i = 0; i < skeleton.branches.size(); i++)
      {
        for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
        {
          sem >> radius;
          skeleton.branches[i].curve[j].skel_radius = radius;
        }
      }
    }

	}

	sem >> str;
	if (str == "Confidence_Sigma")
	{
		sem >> num;
		for (int i = 0; i < num; i++)
		{
			double sigma;
			sem >> sigma;
			samples.vert[i].eigen_confidence = sigma;
		}
	}

	sem >> str;
	if (str == "SkelRadius2")
	{
		sem >> num;

    if (num > 1)
    {
      double radius;
      for (int i = 0; i < skeleton.branches.size(); i++)
      {
        for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
        {
          sem >> radius;
          //skeleton.branches[i].curve[j].skel_radius = radius;
        }
      }
    }

	}

	sem >> str;
	if (str == "Alpha")
	{
		sem >> num;
    double Alpha;
    if (num > 1)
    {
      for (int i = 0; i < skeleton.branches.size(); i++)
      {
        for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
        {
          sem >> Alpha;
          //skeleton.curves[i][j].alpha = Alpha;
        }
      }
    }

	}

	if (!sem.eof())
	{
		sem >> str;
		if (str == "Sample_isVirtual")
		{
			sem >> num;
			for (int i = 0; i < num; i++)
			{
				bool b;
				sem >> b;
				samples.vert[i].is_skel_virtual = b;
			}
		}
	}

	if (!sem.eof())
	{
		sem >> str;
		if (str == "Sample_isBranch")
		{
			sem >> num;
			for (int i = 0; i < num; i++)
			{
				bool b;
				sem >> b;
				samples.vert[i].is_skel_branch = b;
			}
		}
	}

	sem >> str;
	if (str == "Sample_radius")
	{
    sem >> num;
    for (int i = 0; i < num; i++)
    {
      double temp;
      sem >> temp;
      //samples.vert[i].saved_radius = temp;
    }
	}

	sem >> str;
	if (str == "Skel_isVirtual")
	{
		cout << "load Skel_isVirtual" << endl;

		sem >> num;
		bool temp;
		for (int i = 0; i < skeleton.branches.size(); i++)
		{
			for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
			{
				sem >> temp;
				skeleton.branches[i].curve[j].is_skel_virtual = temp;
			}
		}
	}

	sem >> str;
	if (str == "Corresponding_sample_index")
	{
		

		sem >> num;

		cout << "load Corresponding_sample_index " << num <<endl;

		if (num > 1)
		{

			int temp;
			for (int i = 0; i < skeleton.branches.size(); i++)
			{
				for (int j = 0; j < skeleton.branches[i].curve.size(); j++)
				{
					sem >> temp;
					skeleton.branches[i].curve[j].m_index = temp;
				}
			}
		}

	}

	sem >> str;
	if (str == "DSN")
	{
		cout << "load DSN" << endl;

		sem >> num;
		for (int i = 0; i < num; i++)
		{
			CVertex v;
			v.bIsOriginal = false;
			v.is_dual_sample = true;
			v.m_index = i;
			sem >> v.P()[0] >> v.P()[1] >> v.P()[2];
			sem >> v.N()[0] >> v.N()[1] >> v.N()[2];
			dual_samples.vert.push_back(v);
			dual_samples.bbox.Add(v.P());
		}
		dual_samples.vn = dual_samples.vert.size();
	}

	if (!sem.eof())
	{
		sem >> str;
		if (str == "dual_corresponding_index")
		{
			cout << "load dual_corresponding_index" << endl;

			sem >> num;
			for (int i = 0; i < num; i++)
			{
				int index;
				sem >> index;
				samples.vert[i].dual_index = index;
			}
		}
	}

	if (!sem.eof())
	{
		sem >> str;
		if (str == "skel_radius")
		{
			sem >> num;
			cout << "load skel_radius !!! " << num << endl;

			for (int i = 0; i < num; i++)
			{
				float sk_radius;
				sem >> sk_radius;

				//cout << num << "  @@@@" << sk_radius << endl;

				//dual_samples.vert[i].skel_radius = sk_radius;
				samples.vert[i].skel_radius = sk_radius;

				//cout << num << "  %%%% " << sk_radius << endl;

				if (i < 20)
				{
					cout << num << "  skel_radius: " << sk_radius;
				}
				//if (radius > 0)
				//{

				//}
			}
		}
	}

	sem >> str;
	if (str == "SkelPN")
	{
		cout << "load SkelPN" << endl;

		sem >> num;
		for (int i = 0; i < num; i++)
		{
			CVertex v;
			v.bIsOriginal = false;
			v.is_dual_sample = false;
			v.is_skel_point = true;
			v.m_index = i;
			sem >> v.P()[0] >> v.P()[1] >> v.P()[2];
			sem >> v.N()[0] >> v.N()[1] >> v.N()[2];
			skel_points.vert.push_back(v);
			skel_points.bbox.Add(v.P());
		}
		skel_points.vn = skel_points.vert.size();
	}

 	clearCMesh(target_dual_samples);
 	clearCMesh(target_samples);

	cout << original.vert.size() << " "
		<< samples.vert.size() << " "
		<< dual_samples.vert.size() << " "
		<< target_samples.vert.size() << " "
		<< target_dual_samples.vert.size() << " " << endl;

	skeleton.generateBranchSampleMap();
}


void DataMgr::replaceMesh(CMesh& src_mesh, CMesh& target_mesh, bool isOriginal)
{
	clearCMesh(target_mesh);
	for (int i = 0; i < src_mesh.vert.size(); i++)
	{
		CVertex v = src_mesh.vert[i];
		v.bIsOriginal = isOriginal;
		v.m_index = i;
		target_mesh.vert.push_back(v);
	}
	target_mesh.vn = src_mesh.vn;
	target_mesh.bbox = src_mesh.bbox;
}


void DataMgr::replaceMeshDual(CMesh& src_mesh, CMesh& target_mesh, bool is_dual)
{
  clearCMesh(target_mesh);
  for(int i = 0; i < src_mesh.vert.size(); i++)
  {
    CVertex v = src_mesh.vert[i];
    v.m_index = i;
    v.is_dual_sample = is_dual;
    target_mesh.vert.push_back(v);
  }
  target_mesh.vn = src_mesh.vn;
  target_mesh.bbox = src_mesh.bbox;
}


void DataMgr::switchSampleOriginal()
{
	CMesh temp_mesh;
	replaceMesh(original, temp_mesh, false);
	replaceMesh(samples, original, true);
	replaceMesh(temp_mesh, samples, false);
}

void DataMgr::switchSampleDualSample()
{
  CMesh temp_mesh;
   replaceMeshDual(dual_samples, temp_mesh, false);
   replaceMeshDual(samples, dual_samples, true);
   replaceMeshDual(temp_mesh, samples, false);

   for (int i = 0; i < samples.vert.size(); i++)
   {
     CVertex& v = samples.vert[i];
      v.is_dual_sample = false;
//      v.is_fixed_sample = false;
   }

   for (int i = 0; i < dual_samples.vert.size(); i++)
   {
     CVertex& v = dual_samples.vert[i];
     v.is_dual_sample = true;
//      v.is_fixed_sample = false;
   }

//   replaceMeshDual(original, temp_mesh, false);
//   replaceMeshDual(samples, original, true);
//   replaceMeshDual(temp_mesh, samples, false);
}



void DataMgr::loadDefaultSphere()
{
	loadPlyToSample("sphere_wlop.ply");
	
	Point3f center = Point3f(0., 0., 0.);
	default_sphere.clear();
	for (int i = 0; i < samples.vert.size(); i++)
	{
		CVertex& v = samples.vert[i];
		Point3f direction = (v.P() - center).Normalize();

		SphereSlot slot(direction, 0.0);
		default_sphere.addSlot(slot);
	}

	cout << "default slot size: " << default_sphere.size() << endl;
}
