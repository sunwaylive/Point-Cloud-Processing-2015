#include "GLDrawer.h"


GLDrawer::GLDrawer(RichParameterSet* _para)
{
	para = _para;
  generateRandomColorList();
}


GLDrawer::~GLDrawer(void)
{
}

void GLDrawer::updateDrawer(vector<int>& pickList)
{
	bCullFace = para->getBool("Need Cull Points");
	bUseIndividualColor = para->getBool("Show Individual Color");
	useNormalColor = para->getBool("Use Color From Normal");
	useFeatureColor = para->getBool("Show Feature Color");
  useDifferBranchColor = para->getBool("Use Differ Branch Color");

	original_draw_width = para->getDouble("Original Draw Width");
	sample_draw_width = para->getDouble("Sample Draw Width");
	dual_sample_draw_width = para->getDouble("Dual Sample Draw Width");

	original_color = para->getColor("Original Point Color");
	sample_color = para->getColor("Sample Point Color");
	normal_width = para->getDouble("Normal Line Width");
	normal_length  = para->getDouble("Normal Line Length");
	curr_radius = global_paraMgr.data.getDouble("CGrid Radius");

	sample_dot_size = para->getDouble("Sample Dot Size");
	dual_sample_dot_size = para->getDouble("Dual Sample Dot Size");

	original_dot_size = para->getDouble("Original Dot Size");
	normal_color = para->getColor("Normal Line Color");
	dlink_color = para->getColor("DLink Color");
	backface_color = para->getColor("Backface Color");

	feature_color = para->getColor("Feature Color");
	pick_color = para->getColor("Pick Point Color");

	bUseConfidenceColor = global_paraMgr.drawer.getBool("Show Confidence Color");
	bUseSegmentaionColor = global_paraMgr.drawer.getBool("Show Segmentation Color");

	sample_cofidence_color_scale = global_paraMgr.glarea.getDouble("Sample Confidence Color Scale");
	iso_value_shift = global_paraMgr.glarea.getDouble("Point ISO Value Shift");

	skel_bone_color = para->getColor("Skeleton Bone Color");
	skel_node_color = para->getColor("Skeleton Node Color");
	skel_branch_color = para->getColor("Skeleton Branch Color");
	skel_bone_width = para->getDouble("Skeleton Bone Width") / 10000.;
	skel_node_size = para->getDouble("Skeleton Node Size") / 10000.;

	if (!pickList.empty())
	{
		curr_pick_indx = pickList[0];
	}
	

}



void GLDrawer::draw(DrawType type, CMesh* _mesh)
{
	if (!_mesh)
	{
		return;
	}

	bool doPick = para->getBool("Doing Pick");

	int qcnt = 0;
	CMesh::VertexIterator vi;
	for(vi = _mesh->vert.begin(); vi != _mesh->vert.end(); ++vi) 
	{
		if (doPick)
		{
			glLoadName(qcnt);
		}

		if (vi->is_skel_ignore)
		{
			continue;
		}

		Point3f& p = vi->P();      
		Point3f& normal = vi->N();

		CVertex temp_v = (*vi);

		if (type == CIRCLE)
		{
			if (!(bCullFace /*&& !vi->is_dual_sample*/ /*&& !vi->bIsOriginal*/) )
			{
				if (isCanSee(p, normal))
				{
					drawCircle(*vi, false);
				}
				else
				{
					temp_v.N() = vi->N() * -1;
					temp_v.P() = vi->P() + vi->N() * -1e-5;
					drawCircle(temp_v, true);
				}

				/*drawCircle(*vi, false);*/
			}
			else
			{
				if (isCanSee(p, normal))
				{
					drawCircle(*vi, false);
				}
			}
		}

    if(!(bCullFace && !vi->is_dual_sample /*&& !vi->bIsOriginal*/) || isCanSee(p, normal))		
    {
			switch(type)
			{
			case DOT:
				drawDot(*vi);
				break;
// 			case CIRCLE:
// 				drawCircle(*vi, false);
// 				
// 				temp_v.N() = vi->N() * -1;
// 				temp_v.P() = vi->P() + vi->N() * -1e-5;
// 				drawCircle(temp_v, true);
				

				break;
			case QUADE:
				drawQuade(*vi);
				break;
			case NORMAL:
				drawNormal(*vi);
				break;
			case SPHERE:
				drawSphere(*vi);
				break;
			default:
				break;
			}
		}

		if (doPick) 
		{
			qcnt++;
		}
	}

	para->setValue("Doing Pick", BoolValue(false));
}

bool GLDrawer::isCanSee(const Point3f& pos, const Point3f& normal)
{
	return  ( (view_point - pos) * normal >= 0 );
}

GLColor GLDrawer::getColorByType(const CVertex& v)
{
// 	if (bUseConfidenceColor)
// 	{
// 		if (v.m_index < 10 && !v.is_dual_sample && !v.bIsOriginal)
// 		{
// 			cout << v.eigen_confidence << endl;
// 		}
// 		return isoValue2color(v.eigen_confidence, sample_cofidence_color_scale, iso_value_shift, true);
// 	}

	if (v.is_boundary && v.is_fixed_original && bUseSegmentaionColor)
	{
		//return cBlue;
		return cBlack;
	}

	if (v.is_skel_branch && bUseSegmentaionColor)
  	{
  		return cBlue;
  	}

	if (v.is_skel_virtual && bUseSegmentaionColor)
  	{
  		return cGray;
  	}

	if (v.is_boundary && bUseSegmentaionColor)
	{
		//return cBlue;
		return cOrange;
	}

	if (v.is_fixed_original && bUseSegmentaionColor)
	{
		//return cBlue;
		return cGreen;
	}

	if (v.is_boundary && bUseConfidenceColor)
   {
     //return cBlue;
		return cOrange;
   }

	if (v.bIsOriginal)
	{
		return original_color;
	}

	if (!v.is_dual_sample && v.is_fixed_sample)
	{
		return cBlue;
	}

	if (v.is_dual_sample && (/*v.is_skel_branch ||*/ v.is_fixed_sample) /*&& bUseConfidenceColor*/)
	{
		return cGreen;
	}

  if (v.is_dual_sample /*|| v.is_fixed_sample*/ )
  {
    return feature_color;
  }


 	if (bUseConfidenceColor)
 	{
 		return isoValue2color(v.eigen_confidence, sample_cofidence_color_scale, iso_value_shift, true);
 	}


	if (bUseIndividualColor)
	{
		Color4b c = v.C();
		return GLColor(c.X() / 255., c.Y() / 255., c.Z() / 255., 1.);
	}



	if (useNormalColor)
	{
		if (RGB_normals.empty())
		{
			//return GLColor( v.cN()[0] , v.cN()[1] , v.cN()[2] );
			return GLColor( (v.cN()[0] + 1) / 2.0 , (v.cN()[1] + 1) / 2.0 , (v.cN()[2]+ 1) / 2.0 );
		}
		else
		{
			Point3f normal = Point3f(v.cN());
			vector<double> color(3, 0.0);
			for(int i = 0; i < 3; i++)
			{
				color[i] = normal * RGB_normals[i];
				if(color[i] < 0)
				{
					color[i] = 0.;		
				}
				//color[i] = (normal * RGB_normals[i] + 1.5) / 2.0;
			}
			return GLColor( color[0], color[1], color[2]);
		}

	}


// 	else if (v.is_fixed_sample)
// 	{
// 		return feature_color;
// 	}



	return sample_color;
}

void GLDrawer::drawDot(const CVertex& v)
{
	int size;
	if (v.bIsOriginal)
	{
		size = original_dot_size;
	}
	else if (v.is_dual_sample)
	{
		size = dual_sample_dot_size;
	}
	else
	{
		size = sample_dot_size;
	}

	glPointSize(size);
	glBegin(GL_POINTS);

	GLColor color = getColorByType(v);
	glColor4f(color.r, color.g, color.b, 1);
	Point3f p = v.P();
	glVertex3f(p[0], p[1], p[2]);
	glEnd(); 
}

void GLDrawer::drawSphere(const CVertex& v)
{
	double radius;
	if (v.bIsOriginal)
	{
		radius = original_draw_width;
	}
	else if (v.is_dual_sample)
	{
		radius = dual_sample_draw_width;
	}
	else
	{
		radius = sample_draw_width;
	}
	
	GLColor color = getColorByType(v);
	glColor4f(color.r, color.g, color.b, 1);

	Point3f p = v.P();

	glDrawSphere(p, color, radius, 20);
}

void GLDrawer::drawCircle(const CVertex& v, bool is_back)
{
	double radius = sample_draw_width;
	if (v.bIsOriginal)
	{
		radius = original_draw_width;
	}
	else if (v.is_dual_sample)
	{
		radius = dual_sample_draw_width;
	}
	else
	{
		radius = sample_draw_width;
	}

	GLColor color = getColorByType(v);

	Point3f p = v.P(); 
	Point3f normal = v.cN();
	Point3f V1 = v.eigen_vector0;
	Point3f V0 = v.eigen_vector1;

	Point3f temp_V = V1 ^ normal;
	float temp = temp_V * V0;

	glNormal3f(normal[0], normal[1], normal[2]);
	glBegin(GL_POLYGON);

	GLColor color2(backface_color);

	if (is_back)
	{
		glColor4f(color2.r, color2.g, color2.b, 1);
	}
	else
	{
		glColor4f(color.r, color.g, color.b, 1);
	}

	int nn = 25;
	for (int i = 0; i < nn; i++)
	{
		float tt1, tt2;
		float rad = 2*3.1415926/nn*i;
		float a, b; 
		a = radius;
		b = radius;

		tt2 = b * sin(rad);
		tt1 = a * cos(rad);

		Point3f vv;

		if (temp < 0)
		{
			vv = V1 * tt1 + V0 * tt2;
		}
		else
		{
			vv = V1 * tt1 - V0 * tt2;
		}
		glVertex3d(p[0]+vv[0], p[1]+vv[1], p[2]+vv[2]);
	}

	glEnd();
}

void GLDrawer::drawQuade(const CVertex& v)
{
	double h = sample_draw_width;
	if (v.bIsOriginal)
	{
		h = original_draw_width;
	}
	else if (v.is_dual_sample)
	{
		h = dual_sample_draw_width;
	}
	else
	{
		h = sample_draw_width;
	}
	
	GLColor color = getColorByType(v);
	glColor4f(color.r, color.g, color.b, 1);

	Point3f p = v.P(); 
	Point3f normal = v.cN();
	Point3f V1 = v.eigen_vector0;
	Point3f V0 = v.eigen_vector1;

	//const double *m = v.m;
	glBegin(GL_QUADS);
	
	Point3f temp_V = V1 ^ normal;
	float temp = temp_V * V0;

	if (temp > 0)
	{
		Point3f p0 = p + V0 * h;
		Point3f p1 = p + V1 * h;
		Point3f p2 = p + (-V0) * h;
		Point3f p3 = p + (-V1) * h;

		glNormal3d(normal.X(), normal.Y(), normal.Z());
		glVertex3d(p0.X(), p0.Y(), p0.Z());
		glVertex3d(p1.X(), p1.Y(), p1.Z());
		glVertex3d(p2.X(), p2.Y(), p2.Z());
		glVertex3d(p3.X(), p3.Y(), p3.Z());
	}
	else
	{
		//cout << "negtive" << endl;
		Point3f p0 = p + (-V1) * h;
		Point3f p1 = p + (-V0) * h;
		Point3f p2 = p + V1 * h;
		Point3f p3 = p + V0 * h;

		glNormal3d(normal.X(), normal.Y(), normal.Z());
		glVertex3d(p0.X(), p0.Y(), p0.Z());
		glVertex3d(p1.X(), p1.Y(), p1.Z());
		glVertex3d(p2.X(), p2.Y(), p2.Z());
		glVertex3d(p3.X(), p3.Y(), p3.Z());
	}

	glEnd();  
}

void GLDrawer::drawNormal(const CVertex& v)
{
	double width = normal_width;
	double length = normal_length;
	QColor qcolor = normal_color;

	glDisable(GL_LIGHTING);

	glLineWidth(width); 
	GLColor color(qcolor);
	glColor4f(color.r, color.g, color.b, 1);  

	Point3f p = v.P(); 
	Point3f m = v.cN();
//   if (v.eigen_confidence > 0)
//   {
//     m = v.eigen_vector0;
//   }

	glBegin(GL_LINES);	
	glVertex3d(p[0], p[1], p[2]);
	glVertex3f(p[0] + m[0]*length, p[1]+m[1]*length, p[2]+m[2]*length);
	glEnd(); 


	//glBegin(GL_LINES);	
	//glVertex3d(p[0], p[1], p[2]);
	//glVertex3f(p[0] - m[0]*length, p[1] - m[1]*length, p[2] - m[2]*length);
	//glEnd(); 

	glEnable(GL_LIGHTING);
}

void GLDrawer::drawPickedDisk(CMesh* dual_samples, NeighborDisk* disk)
{
  double width = para->getDouble("Sample Draw Width") * 1.2;
  //GLColor dnn_color = para->getColor("Pick Point DNN Color");
  GLColor dnn_color = cGreen;

  glColor3f(dnn_color.r, dnn_color.g, dnn_color.b);

  vector<Point3f> projected_points = disk->projected_points;
  for (int i = 0; i < projected_points.size(); i++)
  {
    Point3f p = projected_points[i];

    glPushMatrix();      
    glTranslatef(p[0], p[1], p[2]);
    glutSolidSphere(width, 10, 10);
    glPopMatrix();
  }

}

void GLDrawer::drawPickedPointOriginalNeighbor(CMesh* samples, CMesh* original, vector<int>& pickList)
{
  double width = para->getDouble("Sample Draw Width") * 2.5;
  //GLColor dnn_color = cOrange;
	GLColor dnn_color = cBlack;

  glColor3f(dnn_color.r, dnn_color.g, dnn_color.b);

  for(int ii = 0; ii < pickList.size(); ii++) 
  {
    int i = pickList[ii];

    if(i < 0 )
      continue;

    CVertex &v = samples->vert[i];
    Point3f &p = v.P();  

    vector<int>::iterator vi = v.original_neighbors.begin();
    int nsize = 6;
    int j = 0;
    for(; vi != v.original_neighbors.end()/* && j < nsize*/; ++vi, ++j)
    {
      if( (*vi) > original->vert.size() || (*vi) < 0)
        break;

// 			if (j < 20)
// 			{
// 				cout << "vi " << *vi << endl;
// 			}

      CVertex& v2 = original->vert[*vi];
      Point3f p2 = v2.P();

      glPushMatrix();      
      glTranslatef(p2[0], p2[1], p2[2]);
      glutSolidSphere(width * 2, 10, 10);
      glPopMatrix();
    }
  } 
}


void GLDrawer::drawPickedPointNeighbor(CMesh* samples, vector<int>& pickList)
{
  double width = para->getDouble("Sample Draw Width") * 2.6;
  //GLColor dnn_color = para->getColor("Pick Point DNN Color");
	GLColor dnn_color = cBlack;

  glColor3f(dnn_color.r, dnn_color.g, dnn_color.b);

  for(int ii = 0; ii < pickList.size(); ii++) 
  {
    int i = pickList[ii];

    if(i < 0 || i >= samples->vert.size())
      continue;

    CVertex &v = samples->vert[i];
    Point3f &p = v.P();  

		if (v.neighbors.empty())
		{
			continue;
		}
    //cout << "Neighbor number: " << v.neighbors.size() << endl;

    vector<int>::iterator vi = v.neighbors.begin();
    int nsize = 6;
    int j = 0;
    for(; vi != v.neighbors.end()/* && j < nsize*/; ++vi, ++j)
    {
      if( (*vi) > samples->vert.size() || (*vi) < 0)
        break;

      CVertex& v2 = samples->vert[*vi];
      Point3f p2 = v2.P();

      glPushMatrix();      
      glTranslatef(p2[0], p2[1], p2[2]);
      glutSolidSphere(width, 10, 10);
      glPopMatrix();
    }


    //Point3d & p3 = v.nb_center / v.neighbors.size();
//     Point3d & p3 = samples->vert[v.fatherIndex].P();
//     glPushMatrix();      
//     glTranslatef(p3[0], p3[1], p3[2]);
//     glutSolidCube(width*1.5);
//     glPopMatrix();
    //cout << "density: " << v.density << endl;
  } 

}

void GLDrawer::drawPickPoint(CMesh* samples, vector<int>& pickList, bool bShow_as_dot)
{
	double width = para->getDouble("Sample Draw Width");
	//GLColor pick_color = para->getColor("Pick Point Color");
  GLColor pick_color = cYellow;

  glColor3f(pick_color.r, pick_color.g, pick_color.b);

	for(int ii = 0; ii < pickList.size(); ii++) 
	{
		int i = pickList[ii];

		if(i < 0 || i >= samples->vert.size())
			continue;

		CVertex &v = samples->vert[i];
		Point3f &p = v.P();     

		if(bShow_as_dot)
		{
			glPointSize(sample_dot_size * 4);
			glBegin(GL_POINTS);

			GLColor color = pick_color;
			glColor4f(color.r, color.g, color.b, 1);

			glVertex3d(p[0], p[1], p[2]);

			glEnd(); 
		}
		else
		{
			glDrawSphere(p, pick_color, sample_draw_width * 2, 40);
		}
	}    
}

void GLDrawer::glDrawLine(Point3f& p0, Point3f& p1, GLColor color, double width)
{
	glColor3f(color.r, color.g, color.b);
	glLineWidth(width);
	glBegin(GL_LINES);
	glVertex3f(p0[0], p0[1], p0[2]);
	glVertex3f(p1[0], p1[1], p1[2]);
	glEnd();
}


void RenderBone(float x0, float y0, float z0, float x1, float y1, float z1, double width = 20)  
{  
	GLdouble  dir_x = x1 - x0;  
	GLdouble  dir_y = y1 - y0;  
	GLdouble  dir_z = z1 - z0;  
	GLdouble  bone_length = sqrt( dir_x*dir_x + dir_y*dir_y + dir_z*dir_z );  
	static GLUquadricObj *  quad_obj = NULL;  
	if ( quad_obj == NULL )  
		quad_obj = gluNewQuadric();  
	gluQuadricDrawStyle( quad_obj, GLU_FILL );  
	gluQuadricNormals( quad_obj, GLU_SMOOTH );  
	glPushMatrix();  
	// 平移到起始点   
	glTranslated( x0, y0, z0 );  
	// 计算长度   
	double  length;  
	length = sqrt( dir_x*dir_x + dir_y*dir_y + dir_z*dir_z );  
	if ( length < 0.0001 ) {   
		dir_x = 0.0; dir_y = 0.0; dir_z = 1.0;  length = 1.0;  
	}  
	dir_x /= length;  dir_y /= length;  dir_z /= length;  
	GLdouble  up_x, up_y, up_z;  
	up_x = 0.0;  
	up_y = 1.0;  
	up_z = 0.0;  
	double  side_x, side_y, side_z;  
	side_x = up_y * dir_z - up_z * dir_y;  
	side_y = up_z * dir_x - up_x * dir_z;  
	side_z = up_x * dir_y - up_y * dir_x;  
	length = sqrt( side_x*side_x + side_y*side_y + side_z*side_z );  
	if ( length < 0.0001 ) {  
		side_x = 1.0; side_y = 0.0; side_z = 0.0;  length = 1.0;  
	}  
	side_x /= length;  side_y /= length;  side_z /= length;  
	up_x = dir_y * side_z - dir_z * side_y;  
	up_y = dir_z * side_x - dir_x * side_z;  
	up_z = dir_x * side_y - dir_y * side_x;  
	// 计算变换矩阵   
	GLdouble  m[16] = { side_x, side_y, side_z, 0.0,  
		up_x,   up_y,   up_z,   0.0,  
		dir_x,  dir_y,  dir_z,  0.0,  
		0.0,    0.0,    0.0,    1.0 };  
	glMultMatrixd( m );  
	// 圆柱体参数   
	GLdouble radius= width;        // 半径   
	GLdouble slices = 40.0;      //  段数   
	GLdouble stack = 3.0;       // 递归次数   
	gluCylinder( quad_obj, radius, radius, bone_length, slices, stack );   
	glPopMatrix();  
}  

void GLDrawer::glDrawCylinder(Point3f& p0, Point3f& p1, GLColor color, double width)
{
	glColor3f(color.r, color.g, color.b);
	RenderBone(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], width);
}


void GLDrawer::glDrawSphere(Point3f& p, GLColor color, double radius, int slide)
{
	if (radius < 0.0001)
		radius = 0.01;

	glColor3f(color.r, color.g, color.b);
	glPushMatrix();      
	glTranslatef(p[0], p[1], p[2]);
	glutSolidSphere(radius, slide, slide);
	glPopMatrix();
}


void GLDrawer::cleanPickPoint()
{
	curr_pick_indx = 0;
	prevPickIndex = 0;
}


void GLDrawer::glDrawBranches(vector<Branch>& branches, GLColor gl_color)
{
  if (useDifferBranchColor && random_color_list.size() < branches.size())
  {
    generateRandomColorList(branches.size() * 2);
  }

	for(int i = 0; i < branches.size(); i++)
	{
		Curve& curve = branches[i].curve;

		if (curve.empty())
		{
			cout << "curve empty!" << endl;
			continue;
		}

		if (curve.size() == 1)
		{
			Point3f p0 = curve[0];
			glDrawSphere(p0, GLColor(cBlue), skel_node_size, 40);

		}
		Point3f p0, p1;
		double size_scale = 0.5;
		//draw nodes
		for (int j = 0; j < curve.size(); j++)
		{
			Point3f p0 = p0 = curve[j].P();

			if (curve[j].is_skel_virtual)
			{
				glDrawSphere(p0, GLColor(feature_color), skel_node_size * 1.1, 40);
			}
      else if (useDifferBranchColor)
      {
        glDrawSphere(p0, GLColor(random_color_list[i]), skel_node_size, 40);
      }
			else
			{
				glDrawSphere(p0, GLColor(skel_node_color), skel_node_size, 40);
			}
		}

		for (int j = 0; j < curve.size()-1; j++)
		{
			p0 = curve[j].P();
			p1 = curve[j+1].P();

			if (skel_bone_width > 0.0003)
			{
        if (useDifferBranchColor)
        {
          glDrawCylinder(p0, p1, random_color_list[i], skel_bone_width);
        }
        else
        {
          glDrawCylinder(p0, p1, gl_color, skel_bone_width);
        }
			}
		}

	}
}

//void GLDrawer::glDrawCurves(vector<Curve>& curves, GLColor gl_color)
//{
//
//
//	for(int i = 0; i < curves.size(); i++)
//	{
//		Curve& curve = curves[i];
//
//		if (curve.empty())
//		{
//			cout << "curve empty!" << endl;
//			continue;
//		}
//
//		if (curve.size() == 1)
//		{
//			Point3f p0 = curve[0];
//			glDrawSphere(p0, GLColor(cBlue), skel_node_size, 40);
//
//		}
//		Point3f p0, p1;
//		double size_scale = 0.5;
//
//		//draw nodes
//		for (int j = 0; j < curve.size(); j++)
//		{
//			Point3f p0 = p0 = curve[j].P();
//
//			if (curve[j].is_skel_virtual)
//			{
//				glDrawSphere(p0, GLColor(feature_color), skel_node_size * 1.1, 40);
//			}
//      else if (useDifferBranchColor)
//      {
//        glDrawSphere(p0, GLColor(random_color_list[i]), skel_node_size, 40);
//      }
//			else
//			{
//				glDrawSphere(p0, GLColor(skel_node_color), skel_node_size, 40);
//			}
//		}
//
//		for (int j = 0; j < curve.size()-1; j++)
//		{
//			p0 = curve[j].P();
//			p1 = curve[j+1].P();
//
//			if (skel_bone_width > 0.0003)
//			{
//        if (useDifferBranchColor)
//        {
//          glDrawCylinder(p0, p1, random_color_list[i], skel_bone_width);
//        }
//        else
//        {
//          glDrawCylinder(p0, p1, gl_color, skel_bone_width);
//        }
//			}
//		}
//
//	}
//}


void GLDrawer::drawCurveSkeleton(Skeleton& skeleton)
{
	if (para->getBool("Skeleton Light"))
	{
		glEnable(GL_LIGHTING);
	}

	glDrawBranches(skeleton.branches, skel_bone_color);

	if (para->getBool("Skeleton Light"))
	{
		glDisable(GL_LIGHTING);
	}
}


void GLDrawer::generateRandomColorList(int num)
{
  if (num <= 1000)
    num = 1000;

  random_color_list.clear();

  srand(time(NULL));
  for(int i = 0; i < num; i++)
  {
    double r = (rand() % 1000) * 0.001;
    double g = (rand() % 1000) * 0.001;
    double b = (rand() % 1000) * 0.001;
    random_color_list.push_back(GLColor(r, g, b));
  }
}

void GLDrawer::drawEigenDirectionsOfone(CVertex& v)
{
	double width = normal_width;
//double length = normal_length;
	double length = curr_radius;

	double eigin_para2 = global_paraMgr.wLop.getDouble("Eigen Neighborhood Para2");
	double eigin_para1 = global_paraMgr.wLop.getDouble("Eigen Neighborhood Para1");

	if (bCullFace && !isCanSee(v.P(), v.N()))
	{
		return;
	}

	Point3f x_end = v.P() + v.eigen_vector0 * length*eigin_para1 + v.eigen_vector0 * v.eigen_value0 * length*eigin_para2;
	Point3f y_end = v.P() + v.eigen_vector1 * length*eigin_para1 + v.eigen_vector1 * v.eigen_value1 * length*eigin_para2;
	Point3f z_end = v.P() + v.eigen_vector2 * length*eigin_para1 + v.eigen_vector2 * v.eigen_value2 * length*eigin_para2;

	//Point3f x_end = v.P() + v.eigen_vector0 * v.eigen_value0 * length;
	//Point3f y_end = v.P() + v.eigen_vector1 * v.eigen_value1 * length;
 	//Point3f z_end = v.P() + v.eigen_vector2 * v.eigen_value2 * length;


	glDrawLine(v.P(), x_end, cRed, width);
	glDrawLine(v.P(), y_end, cGreen, width);
	glDrawLine(v.P(), z_end, cBlue, width);
}

void GLDrawer::drawEigenDirections(CMesh* samples)
{
	double width = normal_width;
	double length = normal_length;
	
	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];

		drawEigenDirectionsOfone(v);


// 		if (bUseConfidenceColor)
// 		{
// 			Point3f diff = v.P() - dual_v.P();
// 			double proj_dist = diff * v.N();
// 			Point3f proj_point = v.P() - v.N() * proj_dist;
// 			glDrawLine(v.P(), proj_point, cBrown, width);
// 		}
// 		else
// 		{
// 			glDrawLine(v.P(), dual_v.P(), cBrown, width);
// 		}


		
	}
}

void GLDrawer::drawDualSampleRelations(CMesh* samples, CMesh* dual_samples)
{
  //if (samples->vn != dual_samples->vn)
  //{
  //  //cout << "No dual sampels!!!" << endl;
  //  return;
  //}

  double width = normal_width;

//   for (int i = 0; i < dual_samples->vert.size(); i++)
//   {
//     CVertex& v = dual_samples->vert[i];
//     //CVertex& dual_v = dual_samples->vert[i];
//     CVertex& dual_v = samples->vert[i];
// 
//     glDrawLine(v.P(), dual_v.P(), cBrown, width);
//   }

  for (int i = 0; i < samples->vert.size(); i++)
  {
    CVertex& v = samples->vert[i];

		if (bCullFace && !isCanSee(v.P(), v.N()))
		{
			continue;
		}

		int index = i;
		if (global_paraMgr.glarea.getBool("Show Cloest Dual Connection"))
		{
			if (v.dual_index > 0 && v.dual_index < dual_samples->vert.size())
			{
				index = v.dual_index;
			}
		}

		CVertex& dual_v = dual_samples->vert[index];
    //CVertex& dual_v = samples->vert[v.dual_index];
 	
		if (bUseConfidenceColor)
		{
			Point3f diff = v.P() - dual_v.P();
			double proj_dist = diff * v.N();
			Point3f proj_point = v.P() - v.N() * proj_dist;
			glDrawLine(v.P(), proj_point, dlink_color, width);
		}
		else
		{
			glDrawLine(v.P(), dual_v.P(), dlink_color, width);
		}


		//glDrawLine(v.P(), dual_v.P(), cBrown, width);
  }

}



GLColor GLDrawer::isoValue2color(double iso_value,
	                               double scale_threshold,
	                               double shift,
	                               bool need_negative)
{
	iso_value += shift;
	//if (!bShowSlice)
	//{
	//  iso_value += shift;
	//}

	if (scale_threshold <= 0)
	{
		scale_threshold = 1;
	}

	iso_value = (std::min)((std::max)(iso_value, -scale_threshold), scale_threshold);

	bool is_inside = true;
	if (iso_value < 0)
	{
		if (!need_negative)
		{
			iso_value = 0;
		}
		else
		{
			iso_value /= -(scale_threshold+1e-10);
		}
	}
	else
	{
		iso_value /= scale_threshold;
		is_inside = false;
	}

	vector<Point3f> base_colors(5);
	double step_size = 1.0 / (base_colors.size() - 1);
	int base_id = iso_value / (step_size + 1e-10);
	if (base_id == base_colors.size() - 1)
	{
		base_id = base_colors.size() - 2;
	}
	Point3f mixed_color;

 	if (is_inside && need_negative)
 	{
 		base_colors[4] = Point3f(0.0, 0.0, 1.0);
 		base_colors[3] = Point3f(0.0, 0.7, 1.0);
 		base_colors[2] = Point3f(0.0, 1.0, 1.0);
 		base_colors[1] = Point3f(0.0, 1.0, 0.7);
 		base_colors[0] = Point3f(0.0, 1.0, 0.0);
 	}
 	else
 	{
 		base_colors[4] = Point3f(1.0, 0.0, 0.0);
 		base_colors[3] = Point3f(1.0, 0.7, 0.0);
 		base_colors[2] = Point3f(1.0, 1.0, 0.0);
 		base_colors[1] = Point3f(0.7, 1.0, 0.0);
 		base_colors[0] = Point3f(0.0, 1.0, 0.0);
 	}


//  	if (is_inside && need_negative)
//  	{
//  		base_colors[4] = Point3f(0.0, 0.0, 1.0);
//  		base_colors[3] = Point3f(0.75, 0.0, 1.0);
//  		base_colors[2] = Point3f(0.5, 0.0, 1.0);
//  		base_colors[1] = Point3f(0.25, 0.0, 1.0);
//  		base_colors[0] = Point3f(1.0, 0.0, 1.0);
//  	}
//  	else
//  	{
//  		base_colors[4] = Point3f(1.0, 0.0, 0.0);
//  		base_colors[3] = Point3f(1.0, 0.0, 0.75);
//  		base_colors[2] = Point3f(1.0, 0.0, 0.5);
//  		base_colors[1] = Point3f(1.0, 0.0, 0.25);
//  		base_colors[0] = Point3f(1.0, 0.0, 1.0);
//  	}

	mixed_color = base_colors[base_id] * (base_id * step_size + step_size - iso_value)
		+ base_colors[base_id + 1] * (iso_value - base_id * step_size);

	return GLColor(mixed_color.X() / float(step_size),
		mixed_color.Y() / float(step_size),
		mixed_color.Z() / float(step_size));
}