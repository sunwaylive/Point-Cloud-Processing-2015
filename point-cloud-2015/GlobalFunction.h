

#pragma once
// #include "TriMesh.h"
// #include "TriMesh_algo.h"
// #include "ICP.h"


#include <vector>
#include "CMesh.h"
#include "grid.h"


// #define WIN32_LEAN_AND_MEAN
// # include <windows.h>
// #include <WindowsX.h>





//#include "LAP_Others/eigen.h"
#include <fstream>
#include <float.h>
#include <QString>
#include <iostream>
#include <time.h>
#include <string>
#include <ctime>
#include<algorithm>
#include <math.h>
#include "ANN/ANN.h"
#include <vcg/space/fitting3.h>
#include <eigenlib/Eigen/Dense>



#define EIGEN_DEFAULT_TO_ROW_MAJOR
#define EIGEN_EXCEPTIONS

//#include <Eigen/Dense>

using namespace std;
using namespace vcg;


//typedef Eigen::MatrixXd Matrix;

#define MyMax(a,b) (((a) > (b)) ? (a) : (b))  
#define MyMin(a,b) (((a) < (b)) ? (a) : (b))  


namespace GlobalFun
{
	void computeKnnNeigbhors(vector<CVertex> &datapts, vector<CVertex> &querypts, int numKnn, bool need_self_included, QString purpose);
	void computeEigen(CMesh* _samples);
	void computeEigenIgnoreBranchedPoints(CMesh* _samples);
	void computeEigenWithTheta(CMesh* _samples, double radius);
  
  void computeUndirectedNormal(CMesh* _samples, float radius);

	void computeAnnNeigbhors(vector<CVertex> &datapts, vector<CVertex> &querypts, int numKnn, bool need_self_included, QString purpose);
	void computeBallNeighbors(CMesh* mesh0, CMesh* mesh1, double radius, vcg::Box3f& box);

	void static  __cdecl self_neighbors(CGrid::iterator start, CGrid::iterator end, double radius);
	void static  __cdecl other_neighbors(CGrid::iterator starta, CGrid::iterator enda, 
		CGrid::iterator startb, CGrid::iterator endb, double radius);
	void static __cdecl find_original_neighbors(CGrid::iterator starta, CGrid::iterator enda, 
		CGrid::iterator startb, CGrid::iterator endb, double radius); 

	double computeEulerDist(Point3f& p1, Point3f& p2);
	double computeEulerDistSquare(Point3f& p1, Point3f& p2);
	double computeProjDist(Point3f& p1, Point3f& p2, Point3f& normal_of_p1);
	double computeProjDistSquare(Point3f& p1, Point3f& p2, Point3f& normal_of_p1);
	double computePerpendicularDistSquare(Point3f& p1, Point3f& p2, Point3f& normal_of_p1);
	double computePerpendicularDist(Point3f& p1, Point3f& p2, Point3f& normal_of_p1);
	double computeProjPlusPerpenDist(Point3f& p1, Point3f& p2, Point3f& normal_of_p1);
	double getDoubleMAXIMUM();
	vector<int> GetRandomCards(int Max);

	double computeRealAngleOfTwoVertor(Point3f v0, Point3f v1);
  double computeDirectionalAngleOfTwoVertor(Point3f v0, Point3f v1, Point3f normal);

	bool isTwoPoint3fTheSame(Point3f& v0, Point3f& v1);
	bool isTwoPoint3fOpposite(Point3f& v0, Point3f& v1);

  void deleteIgnore(CMesh* mesh);
  void removeNormalOverlaps(CMesh* mesh);

  Point3f getTangentVector(Point3f& diff_vecotr, Point3f& normal);

	void normalizeConfidence(vector<CVertex>& vertexes, float delta);

	void addOutliers(CMesh *mesh, int add_num, double max_move_dist);
	void addOutliers(CMesh *mesh, double outlier_percent, double max_move_dist);
	//void computeICP(CMesh *dst, CMesh *src);

	vector<double> computeDensityConfidence(CMesh *mesh, double radius);


	vector<double> computeNormalDifference(CMesh *mesh, double radius, double sigma);
	vector<double> computeBilateralConfidence(CMesh *mesh, double radius, double sigma);

	vector<double> smoothConfidences(CMesh *mesh, double radius);
}

class Timer
{
public:

	void start(const string& str)
	{
		cout << endl;
		starttime = clock();
		mid_start = clock();
		cout << "@@@@@ Time Count Strat For: " << str << endl;

		_str = str;
	}

	void insert(const string& str)
	{
		mid_end = clock();
		timeused = mid_end - mid_start;
		cout << "##" << str << "  time used:  " << timeused / double(CLOCKS_PER_SEC) << " seconds." << endl;
		mid_start = clock();
	}

	void end()
	{
		stoptime = clock();
		timeused = stoptime - starttime;
		cout << /*endl <<*/ "@@@@ finish	" << _str << "  time used:  " << timeused / double(CLOCKS_PER_SEC) << " seconds." << endl;
		cout << endl;
	}

private:
	int starttime, mid_start, mid_end, stoptime, timeused;
	string _str;
};


class SphereSlot
{
public:
	Point3f slot_direction;
	double slot_value;

	SphereSlot()
	{
		slot_direction = Point3f(0., 0., 0.);
		slot_value = 0.0;
	}

	SphereSlot(Point3f dir, double val)
	{
		slot_direction = dir;
		slot_value = val;
	}

	SphereSlot(const SphereSlot& sphere_slot)
	{
		slot_direction = sphere_slot.slot_direction;
		slot_value = sphere_slot.slot_value;
	}

	SphereSlot& SphereSlot::operator = (const SphereSlot& sphere_slot)
	{
		if (&sphere_slot != this)
		{
			slot_direction = sphere_slot.slot_direction;
			slot_value = sphere_slot.slot_value;
		}
		return *this;
	}

	bool SphereSlot::operator < (const SphereSlot& sphere_slot)
	{
		return slot_value < sphere_slot.slot_value;
	}
};

// bool cmpare_slots(SphereSlot &a, SphereSlot &b)
// {
// 	return a.slot_value < b.slot_value;
// }

class SphereSlots
{
public:
	SphereSlots(){}
	~SphereSlots(){}
	vector<SphereSlot> getSphereSlots(){ return sphere_slots; }
	void clear(){ sphere_slots.clear(); }
	void addSlot(SphereSlot& slot){ sphere_slots.push_back(slot); }
	SphereSlot* getSlot(int i) { return &sphere_slots[i]; };
	int size(){ return sphere_slots.size(); }



	void rearrangeSlots(Point3f main_dir)
	{
		for (int i = 0; i < sphere_slots.size(); i++)
		{
			SphereSlot& slot = sphere_slots[i];
			slot.slot_value = GlobalFun::computeEulerDistSquare(main_dir, slot.slot_direction);
		}
		//sort(sphere_slots.begin(), sphere_slots.end(), cmpare_slots);
		sort(sphere_slots.begin(), sphere_slots.end());
	}

	double computeDistance(SphereSlots& another_slots)
	{
		double sum_dist = 0.0;
		for (int i = 0; i < sphere_slots.size(); i++)
		{
			SphereSlot& this_slot = sphere_slots[i];
			SphereSlot& other_slot = *another_slots.getSlot(i);

			double dist = abs(this_slot.slot_value - other_slot.slot_value);
			sum_dist += dist;
		}
		return sum_dist;
	}

	double computeL1Distance(SphereSlots& another_slots)
	{
		vector<double> dist_vec;
		for (int i = 0; i < sphere_slots.size(); i++)
		{
			SphereSlot& this_slot = sphere_slots[i];
			SphereSlot& other_slot = *another_slots.getSlot(i);

			double dist = abs(this_slot.slot_value - other_slot.slot_value);
			dist_vec.push_back(dist);
		}
		sort(dist_vec.begin(), dist_vec.end());
		double mid_val = dist_vec[int(dist_vec.size() / 2)];
		return mid_val;
	}


public:
	vector<SphereSlot> sphere_slots;
};


//typedef vector<SphereSlot> SphereSlots;

class NeighborDisk
{
public:
  NeighborDisk(){}
  NeighborDisk(Point3f c, Point3f n)
  {
    center = c;
    normal = n;
    Point3f temp(0, 0, 1.11111111111);
    zero_axis = (normal ^ temp).Normalize();
    slot_num = 12;
    disk_slots.assign(slot_num, false);
  }

  NeighborDisk(const NeighborDisk& d)
  {
    center = d.center;
    normal = d.normal;
    zero_axis = d.zero_axis;
    disk_slots = d.disk_slots;
    projected_points = d.projected_points;
    slot_num = d.slot_num;
  }

  NeighborDisk& NeighborDisk::operator = (const NeighborDisk& d)
  {
    if (&d != this)
    {
      center = d.center;
      normal = d.normal;
      zero_axis = d.zero_axis;
      disk_slots = d.disk_slots;
      projected_points = d.projected_points;
      slot_num = d.slot_num;
    }
    return *this;
  }

  Point3f projectPointToDisk(Point3f new_p)
  {
    Point3f temp1 = (new_p - center).Normalize();
    Point3f temp2 = temp1 ^ normal;
    Point3f temp3 = normal ^ temp2;

    if (temp1 * temp3 < 0)
    {
      temp3 *= -1;
    }

    double proj_dist = (new_p - center) * temp3;
    Point3f new_projected_point = center + temp3 * proj_dist;
    return new_projected_point;
  }

  int getSlotIDbyPoint(Point3f new_p)
  {
    Point3f temp1 = (new_p - center).Normalize();
    Point3f temp2 = temp1 ^ normal;
    Point3f temp3 = normal ^ temp2;

    if (temp1 * temp3 < 0)
    {
      temp3 *= -1;
    }

    double proj_dist = (new_p - center) * temp3;
    Point3f new_projected_point = center + temp3 * proj_dist;

    //update slots
    double angle = GlobalFun::computeRealAngleOfTwoVertor(temp3, zero_axis);
    Point3f temp_test = (temp3.Normalize() ^ zero_axis).Normalize();
    if (temp_test * normal < 0)
    {
      angle = 360. - angle;
    }

    int new_slot_id = angle / 10.0;
    if (new_slot_id >= slot_num)
    {
      new_slot_id = slot_num-1;
    }
    return new_slot_id;
  }

  void add_point(Point3f new_p)
  {
    Point3f temp1 = (new_p - center).Normalize();
    Point3f temp2 = temp1 ^ normal;
    Point3f temp3 = normal ^ temp2;

    if (temp1 * temp3 < 0)
    {
      temp3 *= -1;
    }

    double proj_dist = (new_p - center) * temp3;
    Point3f new_projected_point = center + temp3 * proj_dist;

    //update slots
    double angle = GlobalFun::computeRealAngleOfTwoVertor(temp3, zero_axis);
    Point3f temp_test = (temp3.Normalize() ^ zero_axis).Normalize();
    if (temp_test * normal < 0)
    {
      angle = 360. - angle;
    }

    int new_slot_id = angle / 10.0;
    if (new_slot_id >= slot_num)
    {
      new_slot_id = slot_num-1;
    }
    projected_points.push_back(new_projected_point);
    disk_slots[new_slot_id] = true;
  }

  bool isSlotOccupied(Point3f new_p)
  {
    int new_slot_id = getSlotIDbyPoint(new_p);
    return disk_slots[new_slot_id];
  }

  bool isNewPointGood(Point3f new_p)
  {
    Point3f new_projected_point = projectPointToDisk(new_p);
    Point3f new_direction = (new_projected_point - center).normalized();
    double min_angle = 180;
    for (int i = 0; i < projected_points.size(); i++)
    {
      Point3f proj_p = projected_points[i];
      Point3f direction = (proj_p - center).normalized();
      double angle = GlobalFun::computeRealAngleOfTwoVertor(new_direction, direction);
      if (angle < min_angle)
      {
        min_angle = angle;
      }
    }

    if (min_angle > (360. / slot_num))
    {
      return true;
    }
    else
    {
      return false;
    }

  }

  double getOccupyPercentage()
  {
    int occ_number = 0;
    for (int i = 0; i < disk_slots.size(); i++)
    {
      if (disk_slots[i])
      {
        occ_number++;
      }
    }
    return double(occ_number) / disk_slots.size();
  }

  void printSlots()
  {
    for (int i = 0; i < disk_slots.size(); i++)
    {
      if (disk_slots[i])
      {
        cout << "1  ";
      }
      else
      {
        cout << "0  ";
      }
    }
    cout << endl;
  }



public:
  Point3f center;
  Point3f normal;
  Point3f zero_axis;
  vector<bool> disk_slots;
  vector<Point3f> projected_points;
  int slot_num;
};


/* Useful code template

(1)
for(int i = 0; i < samples->vert.size(); i++)
{
CVertex& v = samples->vert[i];

for (int j = 0; j < v.neighbors.size(); j++)
{
CVertex& t = samples->vert[v.neighbors[j]];
}
}

(2)
int count = 0;
time.start("Test 2");
CMesh::VertexIterator vi;
Point3f p0 = Point3f(0,0,0);
for(vi = original->vert.begin(); vi != original->vert.end(); ++vi)
{
count += GlobalFun::computeEulerDistSquare(p0, vi->P());
}
cout << count << endl;
time.end();


time.start("Test 1");
for(int i = 0; i < original->vert.size(); i++)
{
CVertex& v = original->vert[i];
count += (p0 - v.P()).SquaredNorm();
}
cout << count << endl;
time.end();



*/

