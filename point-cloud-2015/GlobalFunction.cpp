#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <assert.h>


#include "grid.h"
//#include "LAP_Others/eigen.h"
#include "GlobalFunction.h"
#include "SparseICP.h"


using namespace vcg;
using namespace std;
using namespace tri;


void GlobalFun::find_original_neighbors(CGrid::iterator starta, CGrid::iterator enda, 
	CGrid::iterator startb, CGrid::iterator endb, double radius) 
{	

	double radius2 = radius*radius;
	double iradius16 = -4/radius2;
	const double PI = 3.1415926;

	for(CGrid::iterator dest = starta; dest != enda; dest++) 
	{
		CVertex &v = *(*dest);

		Point3f &p = v.P();
		for(CGrid::iterator origin = startb; origin != endb; origin++)
		{
			CVertex &t = *(*origin);

			Point3f &q = t.P();
			Point3f diff = p-q;

			double dist2 = diff.SquaredNorm();

			if(dist2 < radius2) 
			{                          
				v.original_neighbors.push_back((*origin)->m_index);
			}
		}
	}
}



// get neighbors
void GlobalFun::self_neighbors(CGrid::iterator start, CGrid::iterator end, double radius)
{
	double radius2 = radius*radius;
	for(CGrid::iterator dest = start; dest != end; dest++)
	{
		CVertex &v = *(*dest);
		Point3f &p = v.P();


		for(CGrid::iterator origin = dest+1; origin != end; origin++)
		{
			CVertex &t = *(*origin);
			Point3f &q = t.P();
			Point3f diff = p-q;
			double dist2 = diff.SquaredNorm();
			if(dist2 < radius2) 
			{   
				v.neighbors.push_back((*origin)->m_index);
				t.neighbors.push_back((*dest)->m_index);
			}
		}
	}
}

void GlobalFun::other_neighbors(CGrid::iterator starta, CGrid::iterator enda, 
	CGrid::iterator startb, CGrid::iterator endb, double radius)
{
	double radius2 = radius*radius;
	for(CGrid::iterator dest = starta; dest != enda; dest++)
	{
		CVertex &v = *(*dest);
		Point3f &p = v.P();

		for(CGrid::iterator origin = startb; origin != endb; origin++)
		{
			CVertex &t = *(*origin);
			Point3f &q = t.P();
			Point3f diff = p-q;
			double dist2 = diff.SquaredNorm();
			if(dist2 < radius2) 
			{   
				v.neighbors.push_back((*origin)->m_index);
				t.neighbors.push_back((*dest)->m_index);
			}
		}
	}
}


void GlobalFun::computeBallNeighbors(CMesh* mesh0, CMesh* mesh1, double radius, vcg::Box3f& box)
{
	if (radius < 0.0001)
	{
		cout << "too small grid!!" << endl; 
		return;
	}
	//mesh1 should be original

	//cout << "compute_Bll_Neighbors" << endl;
	//cout << "radius: " << radius << endl;

	CGrid samples_grid;
	samples_grid.init(mesh0->vert, box, radius);
	//cout << "finished init" << endl;

	if (mesh1 != NULL)
	{
		for (int i = 0; i < mesh0->vn; i++)
		{
			mesh0->vert[i].original_neighbors.clear();
		}

		CGrid original_grid;
		original_grid.init(mesh1->vert, box, radius); // This can be speed up
		samples_grid.sample(original_grid, find_original_neighbors);
	}
	else
	{
		for (int i = 0; i < mesh0->vn; i++)
		{
			mesh0->vert[i].neighbors.clear();
		}

		samples_grid.iterate(self_neighbors, other_neighbors);
	}

}


void GlobalFun::computeAnnNeigbhors(vector<CVertex> &datapts, vector<CVertex> &querypts, int knn, bool need_self_included = false, QString purpose = "?_?")
{
	cout << endl <<"Compute ANN for:	 " << purpose.toStdString() << endl;
	int numKnn = knn + 1;

	if (querypts.size() <= numKnn+2)
	{
		vector<CVertex>::iterator vi;
		for(vi = datapts.begin(); vi != datapts.end(); ++vi)
		{
			for(int j = 0; j < 3; j++)
			{
				vi->neighbors.clear();
			}
		}
		return;
	}

	int					nPts;					// actual number of data points
	ANNpointArray		dataPts;				// data points
	ANNpoint			queryPt;				// query point
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		dists;					// near neighbor distances
	ANNkd_tree*			kdTree;					// search structure

	int				k				= numKnn;			// number of nearest neighbors
	int				dim				= 3;			// dimension
	double			eps				= 0;			// error bound
	int				maxPts			= numKnn + 3000000;			// maximum number of data points

	if (datapts.size() >= maxPts)
	{
		cout << "Too many data" << endl;
		return;
	}


	queryPt = annAllocPt(dim);					// allocate query point
	dataPts = annAllocPts(maxPts, dim);			// allocate data points
	nnIdx = new ANNidx[k];						// allocate near neigh indices
	dists = new ANNdist[k];						// allocate near neighbor dists

	nPts = datapts.size();									// read data points

	vector<CVertex>::iterator vi;
	int index = 0;
	for(vi = datapts.begin(); vi != datapts.end(); ++vi)
	{
		for(int j = 0; j < 3; j++)
		{
			dataPts[index][j] = double(vi->P()[j]); 
		}
		index++;
	}


	knn++;

	kdTree = new ANNkd_tree(					// build search structure
		dataPts,					// the data points
		nPts,						// number of points
		dim);						// dimension of space

	knn--;

	for (vi = querypts.begin(); vi != querypts.end(); ++vi) 
	{
		vi->neighbors.clear();
		for (int j = 0; j < 3; j++) 
		{
			queryPt[j] = vi->P()[j];
		}

		kdTree->annkSearch(						// search
			queryPt,						// query point
			k,								// number of near neighbors
			nnIdx,							// nearest neighbors (returned)
			dists,							// distance (returned)
			eps);							// error bound

		for (int k = 1; k < numKnn; k++)
		{
			vi->neighbors.push_back(nnIdx[k]);
		}
	}

	annDeallocPt(queryPt);					// deallocate query point
	annDeallocPts(dataPts);			// deallocate data points
	delete [] nnIdx;							// clean things up
	delete [] dists;
	delete kdTree;
	annClose();									// done with ANN
}


void GlobalFun::computeKnnNeigbhors(vector<CVertex> &datapts, vector<CVertex> &querypts, int numKnn, bool need_self_included = false, QString purpose = "?_?")
{
	if (querypts.size() <= numKnn+1)
	{
		vector<CVertex>::iterator vi;
		for(vi = datapts.begin(); vi != datapts.end(); ++vi)
		{
			for(int j = 0; j < 3; j++)
			{
				vi->neighbors.clear();
			}
		}
		return;
	}

	bool isComputingOriginalNeighbor = false;
	//if (!datapts.empty() && datapts[0].bIsOriginal)
	//{
	//	isComputingOriginalNeighbor = true;
	//}

	int starttime, stoptime, timeused;
	starttime = clock();

	cout << endl;
	cout << "compute KNN Neighbors for: " << purpose.toStdString() << endl;


	ofstream outfile1;
	ofstream outfile2;
	float val;

	outfile1.open("point_cloud.txt", ofstream::binary);
	outfile2.open("query.txt", ofstream::binary);

	val = datapts.size();
	outfile1.write((char *)(&val), sizeof(float));
	val = querypts.size();
	outfile2.write((char *)(&val), sizeof(float));
	val = 3;
	outfile1.write((char *)(&val), sizeof(float));
	val = 4;
	outfile2.write((char *)(&val), sizeof(float));


	vector<CVertex>::iterator vi;
	for(vi = datapts.begin(); vi != datapts.end(); ++vi)
	{
		for(int j = 0; j < 3; j++)
		{
			val = vi->P()[j];
			outfile1.write((char *)(&val), sizeof(float));
		}
	}


	for (vi = querypts.begin(); vi != querypts.end(); ++vi) 
	{
		for (int j = 0; j < 3; j++) 
		{
			val = vi->P()[j];
			outfile2.write((char *)(&val), sizeof(float));
		}
		val = 0;
		outfile2.write((char *)(&val), sizeof(float));
	}

	outfile1.close();
	outfile2.close();

	char mycmd[100];
	sprintf(mycmd, "ANN32//RG_NearestNeighbors.exe point_cloud.txt query.txt result.txt %d", numKnn+1);
	//sprintf(mycmd, "RG_NearestNeighbors.exe point_cloud.txt query.txt result.txt", numKnn+1);

	//cout << mycmd;

	system(mycmd); 

	//cout << "knn_neighbor file saved\n";

	//clean querypts 
	for (vi = querypts.begin(); vi != querypts.end(); ++vi)
	{
		if (isComputingOriginalNeighbor)
		{
			vi->original_neighbors.clear();
		}
		else
		{
			vi->neighbors.clear();
		}
	}

	ifstream infile;
	float size[2];
	int row,col;
	float *data;

	infile.open ("result.txt", ifstream::binary);
	infile.read((char*)size, 2*sizeof(float));
	row = (int)size[0];
	col = (int)size[1];
	data = new float [row*col];
	infile.read((char*)data,row*col*sizeof(float));
	infile.close();

	for (int idx = 0; idx < row; idx++)
	{

		CVertex &v = querypts[(int)data[idx*col+1]-1];
		if (isComputingOriginalNeighbor)
		{
			v.original_neighbors.push_back((int)data[idx*col]-1);
		}
		else
		{
			v.neighbors.push_back((int)data[idx*col]-1);
		}
	}

	if (!need_self_included)// slow solution...
	{
		for(int i = 0; i < querypts.size(); i++)
		{
			CVertex& v = querypts[i];
			v.neighbors.erase(v.neighbors.begin());
		}
	}


	delete[] data;
	//cout << "compute_knn_neighbor end." << endl << endl;

	stoptime = clock();
	timeused = stoptime - starttime;
	cout << "KNN time used:  " << timeused/double(CLOCKS_PER_SEC) << " seconds." << endl;
	cout << endl;
}		


vector<int> GlobalFun::GetRandomCards(int Max)
{
	vector<int> nCard(Max, 0);
	srand(time(NULL));
	for(int i=0; i < Max; i++)
	{
		nCard[i] = i;
	}
	random_shuffle(nCard.begin(), nCard.begin() + Max);


	return nCard;
}


void GlobalFun::computeEigenIgnoreBranchedPoints(CMesh* _samples)
{
	vector<vector<int> > neighborMap;

	typedef vector<CVertex>::iterator VertexIterator;

	VertexIterator begin = _samples->vert.begin();
	VertexIterator end = _samples->vert.end();

	neighborMap.assign(end - begin, vector<int>());

	int curr_index = 0;
	for (VertexIterator iter=begin; iter!=end; ++iter, curr_index++)
	{
    if(iter->neighbors.size() <= 3)
    {
      iter->eigen_confidence = 0.5;
      continue;
    }

		//neighborMap[curr_index].push_back(curr_index);
		for(int j = 0; j < iter->neighbors.size(); j++)
		{
			CVertex& t = _samples->vert[iter->neighbors[j]];
			if (t.is_skel_branch || t.is_skel_ignore)
			{
				continue;
			}
			neighborMap[curr_index].push_back(iter->neighbors[j]);
		}
	}


	int currIndex = 0;
	for (VertexIterator iter=begin; iter!=end; iter++, currIndex++)
	{
		int neighbor_size = neighborMap[currIndex].size();

		if (neighbor_size < 3)
		{
			iter->eigen_confidence = 0.95;
			iter->eigen_vector0 = Point3f(0, 0, 0);

			continue;
		}

		Matrix33d covariance_matrix;
		Point3f diff;
		covariance_matrix.SetZero();
		int neighborIndex = -1;

		for (unsigned int n=0; n<neighbor_size; n++)
		{
			neighborIndex = neighborMap[currIndex][n];
			if(neighborIndex < 0)
				break;
			VertexIterator neighborIter = begin + neighborIndex;

			diff = iter->P() - neighborIter->P();

			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
					covariance_matrix[i][j] += diff[i]*diff[j];
		}

		Point3f   eigenvalues;
		Matrix33d	eigenvectors;
		int required_rotations;
		vcg::Jacobi< Matrix33d, Point3f >(covariance_matrix, eigenvalues, eigenvectors, required_rotations);
		vcg::SortEigenvaluesAndEigenvectors< Matrix33d, Point3f >(eigenvalues, eigenvectors);

		double sum_eigen_value = (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
		iter->eigen_confidence = eigenvalues[0] / sum_eigen_value;

		for (int d=0; d<3; d++)
			iter->eigen_vector0[d] = eigenvectors[d][0];
		for (int d=0; d<3; d++)
			iter->eigen_vector1[d] = eigenvectors[d][1];
		for (int d=0; d<3; d++)
			iter->N()[d] = eigenvectors[d][2];

		iter->eigen_vector0.Normalize();
		iter->eigen_vector1.Normalize();
		iter->N().Normalize();
	}
}

void GlobalFun::computeUndirectedNormal(CMesh* _samples, float radius)
{
  float radius2 = radius * radius;
  radius2 /= 4.0;

  for (int i = 0; i < _samples->vert.size(); i++)
  {
    CVertex& v = _samples->vert[i];

    std::vector<Point3f> ptVec;
    ptVec.push_back(v.P());

    for (int j = 0; j < v.neighbors.size(); j++)
    {
      CVertex& t = _samples->vert[v.neighbors[j]];

      float dist2 = GlobalFun::computeEulerDistSquare(v.P(), t.P());

      if (dist2 < radius2)
      {
        ptVec.push_back(t.P());
      }
    }

    Plane3f plane;
    vcg::FitPlaneToPointSet(ptVec, plane);
    v.N()=plane.Direction();
    v.N().normalized();
  }
}

void GlobalFun::computeEigen(CMesh* _samples)
{
	vector<vector<int> > neighborMap;

	typedef vector<CVertex>::iterator VertexIterator;

	VertexIterator begin = _samples->vert.begin();
	VertexIterator end = _samples->vert.end();

	int curr_index = 0;

	int currIndex = 0;
	for (VertexIterator iter=begin; iter!=end; iter++, currIndex++)
	{
		Matrix33d covariance_matrix;
		Point3f diff;
		covariance_matrix.SetZero();
		int neighbor_size = iter->neighbors.size();
		for (unsigned int n=0; n<neighbor_size; n++)
		{
			Point3f& tP =_samples->vert[iter->neighbors[n]].P();
			diff = iter->P() - tP;

			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
					covariance_matrix[i][j] += diff[i]*diff[j];
		}


		Point3f   eigenvalues;
		Matrix33d	eigenvectors;
		int required_rotations;
		vcg::Jacobi< Matrix33d, Point3f >(covariance_matrix, eigenvalues, eigenvectors, required_rotations);
		vcg::SortEigenvaluesAndEigenvectors< Matrix33d, Point3f >(eigenvalues, eigenvectors);


		double sum_eigen_value = (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);


		iter->eigen_confidence = eigenvalues[0] / sum_eigen_value;

		for (int d=0; d<3; d++)
			iter->eigen_vector0[d] = eigenvectors[d][0];
		for (int d=0; d<3; d++)
			iter->eigen_vector1[d] = eigenvectors[d][1];
		for (int d=0; d<3; d++)
			iter->N()[d] = eigenvectors[d][2];

		iter->eigen_vector0.Normalize();
		iter->eigen_vector1.Normalize();
		iter->N().Normalize();
	}

}




void GlobalFun::computeEigenWithTheta(CMesh* _samples, double radius)
{
	vector<vector<int> > neighborMap;

	typedef vector<CVertex>::iterator VertexIterator;

	VertexIterator begin = _samples->vert.begin();
	VertexIterator end = _samples->vert.end();

	neighborMap.assign(end - begin, vector<int>());

	int curr_index = 0;

	for (VertexIterator iter=begin; iter!=end; iter++, curr_index++)
	{
		if(iter->neighbors.size() <= 3)
		{
			iter->eigen_confidence = 0.5;
			continue;
		}

		for(int j = 0; j < iter->neighbors.size(); j++)
		{
			neighborMap[curr_index].push_back(iter->neighbors[j]);
		}
	}

	double radius2 = radius*radius;
	double iradius16 = -1/radius2; 

	int currIndex = 0;
	for (VertexIterator iter=begin; iter!=end; iter++, currIndex++)
	{
    if(iter->neighbors.size() <= 3)
    {
      iter->eigen_confidence = 0.5;
      continue;
    }

		Matrix33d covariance_matrix;
		Point3f diff;
		covariance_matrix.SetZero();
		int neighborIndex = -1;
		int neighbor_size = iter->neighbors.size();
		for (unsigned int n=0; n<neighbor_size; n++)
		{
			neighborIndex = neighborMap[currIndex][n];
			if(neighborIndex < 0)
				break;
			VertexIterator neighborIter = begin + neighborIndex;

			diff = iter->P() - neighborIter->P();

			Point3f vm = iter->N();
			Point3f tm = neighborIter->N();
			double dist2 = diff.SquaredNorm();
			double theta = exp(dist2*iradius16);

			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
					covariance_matrix[i][j] += diff[i]*diff[j] * theta;
		}

		Point3f   eigenvalues;
		Matrix33d	eigenvectors;
		int required_rotations;
		vcg::Jacobi< Matrix33d, Point3f >(covariance_matrix, eigenvalues, eigenvectors, required_rotations);
		vcg::SortEigenvaluesAndEigenvectors< Matrix33d, Point3f >(eigenvalues, eigenvectors);


		double sum_eigen_value = (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);

		iter->eigen_confidence = eigenvalues[0] / sum_eigen_value;

		for (int d=0; d<3; d++)
			iter->eigen_vector0[d] = eigenvectors[d][0];
		for (int d=0; d<3; d++)
			iter->eigen_vector1[d] = eigenvectors[d][1];
// 		for (int d=0; d<3; d++)
// 			iter->N()[d] = eigenvectors[d][2];

		iter->eigen_vector0.Normalize();
		iter->eigen_vector1.Normalize();
//		iter->N().Normalize();
	}
}



double GlobalFun::computeEulerDist(Point3f& p1, Point3f& p2)
{
	double dist2 = (p1-p2).SquaredNorm();
	if (dist2 < 1e-8 || dist2 > 1e8)
	{
		return 0;
	}
	return sqrt(dist2);
}

double GlobalFun::computeEulerDistSquare(Point3f& p1, Point3f& p2)
{
	return (p1-p2).SquaredNorm();
}

double GlobalFun::computeProjDist(Point3f& p1, Point3f& p2, Point3f& normal_of_p1)
{
	return (p2-p1) * normal_of_p1.Normalize();
}



double GlobalFun::computeProjDistSquare(Point3f& p1, Point3f& p2, Point3f& normal_of_p1)
{
	double proj_dist = computeProjDist(p1, p2, normal_of_p1);
	return proj_dist * proj_dist;
}

double GlobalFun::computePerpendicularDistSquare(Point3f& p1, Point3f& p2, Point3f& normal_of_p1)
{
	//Point3f v_p2_p1 = p1-p2;
	//double proj_dist = computeProjDist(p1, p2, normal_of_p1);
	//Point3f v_proj = /*p1 + */normal_of_p1 * proj_dist;
	//   return (v_p2_p1 + v_proj).SquaredNorm();
	double proj_dist = computeProjDist(p1, p2, normal_of_p1);
	Point3f proj_p = p1 + normal_of_p1 * proj_dist;
	return (proj_p - p2).SquaredNorm();
}



double GlobalFun::computePerpendicularDist(Point3f& p1, Point3f& p2, Point3f& normal_of_p1)
{
	return sqrt(computePerpendicularDistSquare(p1, p2, normal_of_p1));
}

double GlobalFun::computeProjPlusPerpenDist(Point3f& p1, Point3f& p2, Point3f& normal_of_p1)
{
	normal_of_p1.Normalize();
	double proj_dist = GlobalFun::computeProjDist(p1, p2, normal_of_p1);

	if (proj_dist <= 0)
	{
		return -1.;
	}

	Point3f proj_p = p1 + normal_of_p1 * proj_dist;
	double perpend_dist = sqrt((proj_p - p2).SquaredNorm());
	double eular_dist = GlobalFun::computeEulerDist(p1, p2);
	return eular_dist + perpend_dist;
	/*return proj_dist  * 0.5 + perpend_dist;*/
}



double GlobalFun::getDoubleMAXIMUM()
{  
	return (numeric_limits<double>::max)();
}






bool GlobalFun::isTwoPoint3fTheSame(Point3f& v0, Point3f& v1)
{
	if (abs(v0[0] - v1[0]) < 1e-7 &&  
		abs(v0[1] - v1[1]) < 1e-7 && 
		abs(v0[2] - v1[2]) < 1e-7)
	{
		return true;
	}

	return false;

}

bool GlobalFun::isTwoPoint3fOpposite(Point3f& v0, Point3f& v1)
{
	if (abs(-v0[0] - v1[0]) < 1e-7 &&  
		abs(-v0[1] - v1[1]) < 1e-7 && 
		abs(-v0[2] - v1[2]) < 1e-7)
	{
		return true;
	}

	return false;
}

double GlobalFun::computeDirectionalAngleOfTwoVertor(Point3f v0, Point3f v1, Point3f normal)
{
  v0.Normalize();
  v1.Normalize();

  double absolute_angle = computeRealAngleOfTwoVertor(v0, v1);
  
  Point3f dir0_1 = (v0 ^ v1).Normalize();

  if ( (dir0_1 * normal.Normalize()) < 0)
  {
    //return (-absolute_angle);
    return (360-absolute_angle);
  }
  else
  {
    return absolute_angle;
  }
}

double GlobalFun::computeRealAngleOfTwoVertor(Point3f v0, Point3f v1)
{
	v0.Normalize();
	v1.Normalize();

  const double epsilon = 1.0e-6;
  const double nyPI = acos(-1.0);

  double angle_cos = v0 * v1;
  if (fabs(angle_cos-1.0) < epsilon)
  {
    cout << "watch out this angle: if (fabs(angle_cos-1.0) < epsilon)" << endl;
    cout << "return angle 0.0" << endl;
    return 0.0;
  }
  else if(fabs(angle_cos+1.0) < epsilon)
  {
    cout << "watch out this angle else if(fabs(angle_cos+1.0) < epsilon)" << endl;

    return 180;
  }

	if (isTwoPoint3fTheSame(v0, v1))
	{
    cout << "watch out this angle if (isTwoPoint3fTheSame(v0, v1))" << endl;

		return 0;
	}

	if (isTwoPoint3fOpposite(v0, v1))
	{
    cout << "watch out this angle if (isTwoPoint3fOpposite(v0, v1))" << endl;

		return 180;
	}

	if (angle_cos > 1)
	{
    cout << "watch out this angle if (angle_cos > 1)" << endl;

		angle_cos = 0.99;
	}
	if (angle_cos < -1)
	{
    cout << "watch out this angle if (angle_cos < -1)" << endl;

		angle_cos = -0.99;
	}
	if (angle_cos > 0 && angle_cos < 1e-8)
	{
    cout << "watch out this angle if (angle_cos > 0 && angle_cos < 1e-8)" << endl;

		return 90;
	}

	double angle = acos(angle_cos) * 180. / 3.1415926 ;

	if (angle < 0 || angle > 180)
	{
    cout << "watch out this angle iif (angle < 0 || angle > 180)" << endl;

		//cout << "compute angle wrong!!" << endl;
		//system("Pause");
		return 180;
	}



	return angle;
}

Point3f GlobalFun::getTangentVector(Point3f& diff_vecotr, Point3f& normal)
{
  double proj_dist = diff_vecotr * normal.Normalize();
  Point3f proj_vector = normal * proj_dist;
  Point3f tangent_vector = diff_vecotr - proj_vector;

  return tangent_vector;
}

void GlobalFun::deleteIgnore(CMesh* mesh)
{
  vector<CVertex> temp_vert;
  for (int i = 0; i < mesh->vert.size(); i++)
  {
    CVertex& v = mesh->vert[i];
    if (!v.is_skel_ignore)
    {
      temp_vert.push_back(v);
    }
  }

  mesh->vert.clear();
  for (int i = 0; i < temp_vert.size(); i++)
  {
    CVertex& v = temp_vert[i];
    v.m_index = i;
    mesh->vert.push_back(v);
    mesh->bbox.Add(v.P());
  }
  mesh->vn = mesh->vert.size();
}

void GlobalFun::removeNormalOverlaps(CMesh* mesh)
{

  bool has_normal_overlap = true;
  while (has_normal_overlap)
  {
    GlobalFun::computeBallNeighbors(mesh, NULL, 0.2, mesh->bbox);

    for (int i = 0; i < mesh->vert.size(); i++)
    {
      CVertex& v = mesh->vert[i];
      has_normal_overlap = false;
      for (int j = 0; j < v.neighbors.size(); j++)
      {
        CVertex& t = mesh->vert[v.neighbors[j]];

        double angle = computeRealAngleOfTwoVertor(v.N(), t.N());
        if (angle < 0.001)
        {
          cout << "delete angle!!" << endl;
          if (v.neighbors[j] != 0)
          {
            t.is_skel_ignore = true;
            has_normal_overlap = true;
            break;
          }
          else
          {
            cout << "the first point overlap!!" << endl;
          }
        }
      }

      if (has_normal_overlap)
      {
        break;
      }
    }

    if (has_normal_overlap)
    {
      deleteIgnore(mesh);
    }
  }

}


void GlobalFun::normalizeConfidence(vector<CVertex>& vertexes, float delta)
{
	float min_confidence = GlobalFun::getDoubleMAXIMUM();
	float max_confidence = 0;
	for (int i = 0; i < vertexes.size(); i++)
	{
		CVertex& v = vertexes[i];
		min_confidence = min_confidence < v.eigen_confidence ? min_confidence : v.eigen_confidence;
		max_confidence = max_confidence > v.eigen_confidence ? max_confidence : v.eigen_confidence;
	}
	float space = max_confidence - min_confidence;

	for (int i = 0; i < vertexes.size(); i++)
	{
		CVertex& v = vertexes[i];
		v.eigen_confidence = (v.eigen_confidence - min_confidence) / space;

		v.eigen_confidence += delta;

		if (!(v.eigen_confidence > 0 || v.eigen_confidence <= 1.0))
		{
			v.eigen_confidence = 0.0;
		}
	}

}


void GlobalFun::addOutliers(CMesh *mesh, int add_num, double max_move_dist)
{
	assert(mesh != NULL);
	for (int i = 0; i < add_num; ++i){
		int idx = rand() % mesh->vert.size();
		CVertex v = mesh->vert[idx];

		Point3f move_dir = Point3f(rand() * 1.0f / RAND_MAX, rand() * 1.0f / RAND_MAX, rand() * 1.0f / RAND_MAX)
			- Point3f(0.5f, 0.5f, 0.5f);
		v.P() += move_dir * max_move_dist * (rand() / RAND_MAX + 0.5);
		mesh->vert.push_back(v);
	}

	for (int i = 0; i < mesh->vert.size(); i++)
	{
		mesh->vert[i].m_index = i;
	}
	mesh->vn = mesh->vert.size();
}

void GlobalFun::addOutliers(CMesh *mesh, double outlier_percent, double max_move_dist)
{
	cout << "add outliers" << endl;
	assert(mesh != NULL);
	int outlier_num = outlier_percent * mesh->vert.size();
	GlobalFun::addOutliers(mesh, outlier_num, max_move_dist);
}


vector<double> GlobalFun::computeDensityConfidence(CMesh *mesh, double radius)
{
	double radius2 = radius * radius;
	double iradius16 = -4 / radius2;

	for (int i = 0; i < mesh->vert.size(); i++)
	{
		CVertex& v = mesh->vert[i];
		vector<int>* neighbors = &v.neighbors;
		v.eigen_confidence = 1.0;

		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& t = mesh->vert[(*neighbors)[j]];
			double dist2 = (v.P() - t.P()).SquaredNorm();
			double den = exp(dist2*iradius16);

			v.eigen_confidence += den;
		}
	}

	normalizeConfidence(mesh->vert, 0.0);

	vector<double> confidences(mesh->vert.size());
	for (int i = 0; i < mesh->vert.size(); i++)
	{
		confidences[i] = mesh->vert[i].eigen_confidence;
	}

	return confidences;
}

vector<double> GlobalFun::computeNormalDifference(CMesh *mesh, double radius, double sigma)
{
	vector<double> confidences(mesh->vert.size());
	for (int i = 0; i < mesh->vert.size(); i++)
	{
		confidences[i] = mesh->vert[i].eigen_confidence;
	}

	return confidences;
}

vector<double> GlobalFun::computeBilateralConfidence(CMesh *mesh, double radius, double sigma)
{
	double sigma_threshold = pow(max(1e-8, 1 - cos(sigma / 180.0*3.1415926)), 2);

	double radius2 = radius * radius;
	double iradius16 = -4 / radius2;

	for (int i = 0; i < mesh->vert.size(); i++)
	{
		CVertex& v = mesh->vert[i];
		vector<int>* neighbors = &v.neighbors;
		v.eigen_confidence = 0.0;

		for (int j = 0; j < v.neighbors.size(); j++)
		{
			CVertex& t = mesh->vert[(*neighbors)[j]];
			double dist2 = (v.P() - t.P()).SquaredNorm();
			double den = exp(dist2*iradius16);
			double normal_diff = exp(-pow(1 - v.N()*t.N(), 2) / sigma_threshold);

			v.eigen_confidence += den * normal_diff;
		}
	}

	normalizeConfidence(mesh->vert, 0.0);

	vector<double> confidences(mesh->vert.size());
	for (int i = 0; i < mesh->vert.size(); i++)
	{
		confidences[i] = mesh->vert[i].eigen_confidence;
	}
	return confidences;
}


//  void
//  GlobalFun::computeICP(CMesh *dst, CMesh *src)
//  {
//  	if (dst->vert.empty() || src->vert.empty())
//  	{
//  		cout << "compute ICP Error : " << "Empty meshes!" << endl;
//  		return;
//  	}
//  
//  	int verbose = 0;
//  	bool do_scale = false;
//  	bool do_affine = false;
//  	bool bulkmode = false;
//  
//  	int c;
//  
//  	princeton::TriMesh::set_verbose(verbose);
//  
//  	//******************** for test *****************************
//  	//const char *filename1 = "1.ply", *filename2 = "2.ply";
//  	//princeton::TriMesh *raw1 = princeton::TriMesh::read(filename1);
//  	//princeton::TriMesh *raw2 = princeton::TriMesh::read(filename2);
//  
//  	//CMesh *cm1 = new CMesh;
//  	//CMesh *cm2 = new CMesh;
//  	//int index = 0;
//  	//for (int i = 0; i < raw1->vertices.size(); ++i)
//  	//{
//  	//  point &p = raw1->vertices[i];
//  	//  CVertex t;
//  	//  t.m_index = index++;
//  	//  t.P()[0] = p[0];
//  	//  t.P()[1] = p[1];
//  	//  t.P()[2] = p[2];
//  	//  cm1->vert.push_back(t);
//  	//  cm1->bbox.Add(t.P());
//  	//}
//  	//cm1->vn = cm1->vert.size();
//  
//  	//for (int j = 0; j < raw2->vertices.size(); ++j)
//  	//{
//  	//  point &p = raw2->vertices[j];
//  	//  CVertex t;
//  	//  t.P()[0] = p[0];
//  	//  t.P()[1] = p[1];
//  	//  t.P()[2] = p[2];
//  	//  cm2->vert.push_back(t);
//  	//  cm2->bbox.Add(t.P());
//  	//}
//  	//cm2->vn = cm2->vert.size();
//  	//*********************************************************
//  
//  	princeton::TriMesh *mesh1 = new princeton::TriMesh;
//  	princeton::TriMesh *mesh2 = new princeton::TriMesh;
//  
//  	//set two tri-meshes
//  	for (int i = 0; i < dst->vert.size(); ++i)
//  	{
//  		CVertex &v = dst->vert[i];
//  		point p(v.P()[0], v.P()[1], v.P()[2]);
//  		mesh1->vertices.push_back(p);
//  	}
//  
//  	for (int j = 0; j < src->vert.size(); ++j)
//  	{
//  		CVertex &v = src->vert[j];
//  		point p(v.P()[0], v.P()[1], v.P()[2]);
//  		mesh2->vertices.push_back(p);
//  	}
//  
//  	xform xf1, xf2;
//  
//  	KDtree *kd1 = new KDtree(mesh1->vertices);
//  	KDtree *kd2 = new KDtree(mesh2->vertices);
//  	vector<float> weights1, weights2;
//  
//  	if (bulkmode) {
//  		float area1 = mesh1->stat(princeton::TriMesh::STAT_TOTAL, princeton::TriMesh::STAT_FACEAREA);
//  		float area2 = mesh2->stat(princeton::TriMesh::STAT_TOTAL, princeton::TriMesh::STAT_FACEAREA);
//  		float overlap_area, overlap_dist;
//  		find_overlap(mesh1, mesh2, xf1, xf2, kd1, kd2,
//  			overlap_area, overlap_dist);
//  		float frac_overlap = overlap_area / min(area1, area2);
//  		if (frac_overlap < 0.1f) {
//  			printf("Insufficient overlap\n");
//  			exit(1);
//  		}
//  		else {
//  			printf("%.1f%% overlap\n",
//  				frac_overlap * 100.0);
//  		}
//  	}
//  
//  	float err = ICP(mesh1, mesh2, xf1, xf2, kd1, kd2, weights1, weights2,
//  		verbose, do_scale, do_affine);
//  	if (err >= 0.0f)
//  		err = ICP(mesh1, mesh2, xf1, xf2, kd1, kd2, weights1, weights2,
//  		verbose, do_scale, do_affine);
//  
//  	if (err < 0.0f) {
//  		printf("ICP failed\n");
//  		exit(1);
//  	}
//  
//  	printf("ICP succeeded - distance = %f\n", err);
//  
//  	//add new points to dst
//  	int index = (dst->vert.empty()) ? 0 : (dst->vert.back()).m_index;
//  	cout << "original index: " << index << endl;;
//  	for (int i = 0; i < mesh2->vertices.size(); ++i)
//  	{
//  		point p = xf2 * mesh2->vertices[i];
//  		CVertex& t = src->vert[i];
//  
//  		t.P()[0] = p[0];
//  		t.P()[1] = p[1];
//  		t.P()[2] = p[2];
//  		//CVertex t;
//  		//t.m_index = ++index;
//  		//t.P()[0] = p[0];
//  		//t.P()[1] = p[1];
//  		//t.P()[2] = p[2];
//  		//dst->vert.push_back(t);
//  		//dst->bbox.Add(t.P());
//  	}
//  	dst->vn = dst->vert.size();
//  }
