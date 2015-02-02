#include "GLArea.h"

GLArea::GLArea(QWidget *parent): QGLWidget(/*QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer |QGL::SampleBuffers),*/ parent),
	               para(global_paraMgr.getGlareaParameterSet()),
								 glDrawer(global_paraMgr.getDrawerParameterSet()),
								 dataMgr(global_paraMgr.getDataParameterSet()),
								 wlop(global_paraMgr.getWLopParameterSet()),
								 norSmoother(global_paraMgr.getNormalSmootherParameterSet()),
								 skeletonization(global_paraMgr.getSkeletonParameterSet()),
								 upsampler(global_paraMgr.getUpsamplingParameterSet()),
								 paintMutex(QMutex::NonRecursive)
{
	setMouseTracking(true); 
	isDragging = false;
	isRightPressed = false;

	trackball_light.center=Point3f(0, 0, 0);
	trackball_light.radius= 1;
	activeDefaultTrackball=true;

	fov = 60;
	clipRatioFar = 1;
	clipRatioNear = 1;
	nearPlane = .2f;
	farPlane = 5.f;
	takeSnapTile = false;
	default_snap_path = QString(".\\snapshot\\");
	current_snap_path = default_snap_path;
	snapDrawScal = 1;
	is_paintGL_locked = false;
	RGB_counter = 0;

	cout << "GLArea constructed" << endl;

	CVertex v;
	cout << sizeof(v) << endl;
	loadDefaultModel();

	show_something = true;
	something = 0;

	need_rotate = false;
	rotate_angle = 0.;
	rotate_pos = Point3f(0, 0, 0);
	rotate_normal = Point3f(1, 0, 0.);
	rotate_delta = 1;
	rotate_angle = 0.;
}

GLArea::~GLArea(void)
{
}


void GLArea::initializeGL()
{
	cout << "initializeGL" << endl;

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0); 
	glEnable(GL_LIGHT1);
	glEnable(GL_COLOR_MATERIAL);  

	glShadeModel (GL_SMOOTH);
	glShadeModel(GL_FLAT);

	//glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST ); 
 	GLfloat mat_ambient[4] = {0.1745, 0.01175, 0.01175,1.0}; 
 	GLfloat mat_diffuse[4] = {0.61424, 0.04136, 0.04136, 1.0 };
 	GLfloat mat_specular[] = {0.727811, 0.626959, 0.626959, 1.0 };

	//GLfloat mat_ambient[4] = { 0.1745, 0.51175, 0.51175, 1.0 };
	//GLfloat mat_diffuse[4] = { 0.61424, 0.54136, 0.54136, 1.0 };
	//GLfloat mat_specular[] = { 0.727811, 0.626959, 0.626959, 1.0 };

	GLfloat shininess = 0.6*128;

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse); 
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, shininess);
	
	glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT );
	glColorMaterial( GL_FRONT_AND_BACK, GL_DIFFUSE );
	//glColorMaterial(GL_FRONT, GL_AMBIENT);
	//glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	glDisable(GL_BLEND);
	glEnable(GL_NORMALIZE);

	glDisable(GL_CULL_FACE);
	glColor4f(1, 1, 1, 1);   

	glEnable(GL_LIGHTING);

	initLight();


	trackball.center = Point3f(0, 0, 0);
	trackball.radius = 1;

  // force to open anti
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_POLYGON_SMOOTH);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);


	glLoadIdentity(); 
}

void GLArea::initLight()
{
	GLColor color_ambient(para->getColor("Light Ambient Color"));
	float ambient[4] = {color_ambient.r, color_ambient.g, color_ambient.b, 1.0};
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);

	GLColor color_diffuse(para->getColor("Light Diffuse Color"));
	float diffuse[4] = {color_diffuse.r, color_diffuse.g, color_diffuse.b, 1.0};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);

	GLColor color_specular(para->getColor("Light Specular Color"));
	float specular[4] = {color_specular.r, color_specular.g, color_specular.b, 1.0};
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular);

}

void GLArea::resizeGL(int w, int h)
{
	//cout << "resizeGL" << endl;

	glViewport(0, 0, (GLint)w, (GLint)h);  
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();  

	float r = w/(float)h;
	gluPerspective(60, r, 0.1, 10);
	glMatrixMode(GL_MODELVIEW);  
}



void GLArea::paintGL() 
{
	if (show_something)
	{
		cout << "show something: 1 " << something++ << endl;
	}
	paintMutex.lock();{

	if (is_paintGL_locked)
	{
		goto PAINT_RETURN;
	}

	lightOnOff(para->getBool("Light On or Off"));

	GLColor color(global_paraMgr.drawer.getColor("Background Color"));
	glClearColor(color.r, color.g, color.b, 1); 

	Point3f lightPos = para->getPoint3f("Light Position");
	float lpos[4];
	lpos[0] = lightPos[0];
	lpos[1] = lightPos[1];
	lpos[2] = lightPos[2];
	lpos[3] = 0;
	glLightfv(GL_LIGHT0, GL_POSITION, lpos);       
	lpos[0] = -lightPos[0];
	lpos[1] = -lightPos[1];
	lpos[2] = -lightPos[2];
	lpos[3] = 0;
	glLightfv(GL_LIGHT1, GL_POSITION, lpos);       
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);   

	double SnapResolutionScale = global_paraMgr.glarea.getDouble("Snapshot Resolution");
	if (takeSnapTile)
	{
		double normal_width = global_paraMgr.drawer.getDouble("Normal Line Width");
		double dot_size = global_paraMgr.drawer.getDouble("Sample Dot Size");
		double original_dot_size = global_paraMgr.drawer.getDouble("Original Dot Size");
		double dual_dot_size = global_paraMgr.drawer.getDouble("Dual Sample Dot Size");

		global_paraMgr.drawer.setValue("Normal Line Width", DoubleValue(normal_width * SnapResolutionScale * SnapResolutionScale * snapDrawScal));
		global_paraMgr.drawer.setValue("Sample Dot Size", DoubleValue(dot_size * SnapResolutionScale * SnapResolutionScale * snapDrawScal));

		global_paraMgr.drawer.setValue("Dual Sample Dot Size", DoubleValue(dual_dot_size * SnapResolutionScale * SnapResolutionScale * snapDrawScal));

		global_paraMgr.drawer.setValue("Original Dot Size", DoubleValue(original_dot_size * SnapResolutionScale * SnapResolutionScale * snapDrawScal));
	}

	glLoadIdentity();
	gluLookAt(0, 0, -3, 0, 0, 0, 0, 1, 0); 

	if (takeSnapTile)
	{
		setView();// high resolution snapshot
		something = false;
	}

	drawLightBall();

	//Drawing the scene 
	glPushMatrix();
	trackball.GetView();

	Point3f c = -gl_box.Center();
	double radius = 2.0f/gl_box.Diag();

	trackball.Apply(false);

	glPushMatrix();
	glScalef(radius, radius, radius);
	glTranslatef(c[0], c[1], c[2]);

	glTranslatef(rotate_pos[0], rotate_pos[1], rotate_pos[2]);
	glRotatef(rotate_angle, rotate_normal[0], rotate_normal[1], rotate_normal[2]);
	glTranslatef(-rotate_pos[0], -rotate_pos[1], -rotate_pos[2]);


	View<float> view;
	view.GetView();
	Point3f viewpoint = view.ViewPoint();
	glDrawer.setViewPoint(viewpoint);

	if (dataMgr.isSamplesEmpty() && dataMgr.isOriginalEmpty())
	{
		if (show_something)
		{
			cout << "empty！";
		}
		show_something = false;

		goto PAINT_RETURN;
	}

	glDrawer.updateDrawer(pickList);


  // have problems...
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   if (para->getBool("Pick Dual Point") && global_paraMgr.drawer.getBool("Doing Pick"))
   {
     glDrawer.draw(GLDrawer::DOT, dataMgr.getCurrentDualSamples());
   }

	if(para->getBool("Show Samples"))  
	{
		if (para->getBool("Show Samples Quad"))
			glDrawer.draw(GLDrawer::QUADE, dataMgr.getCurrentSamples());
		if (para->getBool("Show Samples Dot"))
		{
			lightOnOff(false);
			glDrawer.draw(GLDrawer::DOT, dataMgr.getCurrentSamples());
		}
		if (para->getBool("Show Samples Circle"))
			glDrawer.draw(GLDrawer::CIRCLE, dataMgr.getCurrentSamples());
		if (para->getBool("Show Samples Sphere"))
			glDrawer.draw(GLDrawer::SPHERE, dataMgr.getCurrentSamples());
	}

	if (show_something)
	{
		cout << "show something: 2 " << something++ << endl;
	}

  if(para->getBool("Show Dual Samples"))
  {
		if (para->getBool("Show Dual Samples Quad"))
      glDrawer.draw(GLDrawer::QUADE, dataMgr.getCurrentDualSamples());
		if (para->getBool("Show Dual Samples Dot"))
    {
      lightOnOff(false);
      glDrawer.draw(GLDrawer::DOT, dataMgr.getCurrentDualSamples());
    }
		if (para->getBool("Show Dual Samples Circle"))
      glDrawer.draw(GLDrawer::CIRCLE, dataMgr.getCurrentDualSamples());	
		if (para->getBool("Show Dual Samples Sphere"))
      glDrawer.draw(GLDrawer::SPHERE, dataMgr.getCurrentDualSamples());	
  }

	if (para->getBool("Show Skeltal Points"))
	{
		if (para->getBool("Show Dual Samples Quad"))
			glDrawer.draw(GLDrawer::QUADE, dataMgr.getCurrentSkelPoints());
		if (para->getBool("Show Dual Samples Dot"))
		{
			lightOnOff(false);
			glDrawer.draw(GLDrawer::DOT, dataMgr.getCurrentSkelPoints());
		}
		if (para->getBool("Show Dual Samples Circle"))
			glDrawer.draw(GLDrawer::CIRCLE, dataMgr.getCurrentSkelPoints());
		if (para->getBool("Show Dual Samples Sphere"))
			glDrawer.draw(GLDrawer::SPHERE, dataMgr.getCurrentSkelPoints());
	}

	if (show_something)
	{
		cout << "show something: 3 " << something++ << endl;
	}

	if (para->getBool("Show Target Samples"))
	{
		if (para->getBool("Show Samples Quad"))
			glDrawer.draw(GLDrawer::QUADE, dataMgr.getCurrentTargetSamples());
		if (para->getBool("Show Samples Dot"))
		{
			lightOnOff(false);
			glDrawer.draw(GLDrawer::DOT, dataMgr.getCurrentTargetSamples());
		}
		if (para->getBool("Show Samples Circle"))
			glDrawer.draw(GLDrawer::CIRCLE, dataMgr.getCurrentTargetSamples());
		if (para->getBool("Show Samples Sphere"))
			glDrawer.draw(GLDrawer::SPHERE, dataMgr.getCurrentTargetSamples());
	}

	if (para->getBool("Show Target Dual Samples"))
	{
		if (para->getBool("Show Dual Samples Quad"))
			glDrawer.draw(GLDrawer::QUADE, dataMgr.getCurrentTargetDualSamples());
		if (para->getBool("Show Dual Samples Dot"))
		{
			lightOnOff(false);
			glDrawer.draw(GLDrawer::DOT, dataMgr.getCurrentTargetDualSamples());
		}
		if (para->getBool("Show Dual Samples Circle"))
			glDrawer.draw(GLDrawer::CIRCLE, dataMgr.getCurrentTargetDualSamples());
		if (para->getBool("Show Dual Samples Sphere"))
			glDrawer.draw(GLDrawer::SPHERE, dataMgr.getCurrentTargetDualSamples());
	}

  lightOnOff(para->getBool("Light On or Off"));


	if (para->getBool("Show Dual Connection"))
	{
		//glDrawer.drawDualSampleRelations(dataMgr.getCurrentSamples(), dataMgr.getCurrentDualSamples());
		glDrawer.drawDualSampleRelations(dataMgr.getCurrentSamples(), dataMgr.getCurrentSkelPoints());

		//glDrawer.drawDualSampleRelations(dataMgr.getCurrentTargetSamples(), dataMgr.getCurrentTargetDualSamples());
	}

	if (show_something)
	{
		cout << "show something: 4 " << something++ << endl;
	}

	if (para->getBool("Show Eigen Directions"))
	{
// 		//glDrawer.drawEigenDirections(dataMgr.getCurrentDualSamples());
// 
// 		CVertex v;
// 		if (!pickList.empty() && pickList[0] >= 0)
// 		{
// 			int id = pickList[0];
// 			if (id >= 0 && id < dataMgr.getCurrentDualSamples()->vert.size())
// 			{
// 				v = dataMgr.getCurrentDualSamples()->vert[id];
// 			}
// 			else
// 			{
// 				v = dataMgr.getCurrentDualSamples()->vert[0];
// 			}
// 		}
// 		else
// 		{
// 			v = dataMgr.getCurrentDualSamples()->vert[0];
// 		}
// 		glDrawer.drawEigenDirectionsOfone(v);
// 
// //  		glw.m = dataMgr.getCurrentEllipsoid();
// //  		glw.Draw(GLW::DMSmooth, GLW::CMPerMesh, GLW::TMNone);

	}

	if (show_something)
	{
		cout << "show something: 5 " << something++ << endl;
	}


	if (para->getBool("Show Normal")) 
	{
		if (para->getBool("Show Samples"))
		{
			glDrawer.draw(GLDrawer::NORMAL, dataMgr.getCurrentSamples());
		}
    else if (para->getBool("Show Dual Samples"))
    {
      glDrawer.draw(GLDrawer::NORMAL, dataMgr.getCurrentDualSamples());
    }
		else
		{
			if(!dataMgr.isOriginalEmpty())
			  glDrawer.draw(GLDrawer::NORMAL, dataMgr.getCurrentOriginal());
		}
	}


	if (show_something)
	{
		cout << "show something: " << something++ << endl;
	}

	if (para->getBool("Show Skeleton"))
	{
		CMesh* dual_samples = dataMgr.getCurrentDualSamples();

		for (int i = 1; i < dual_samples->vert.size(); i++)
		{
			Point3f previous_position = dual_samples->vert[i-1].P();
			Point3f current_position = dual_samples->vert[i].P();

			glDrawer.glDrawLine(previous_position, current_position, cBlue, 1);

		}
	}

 	if(para->getBool("Show Original"))
 	{
 		if(!dataMgr.isOriginalEmpty())
 		{
 			if(para->getBool("Show Original Quad"))
 				glDrawer.draw(GLDrawer::QUADE, dataMgr.getCurrentOriginal());
 			if(para->getBool("Show Original Dot"))
 				glDrawer.draw(GLDrawer::DOT, dataMgr.getCurrentOriginal());
 			if(para->getBool("Show Original Circle"))
 				glDrawer.draw(GLDrawer::CIRCLE, dataMgr.getCurrentOriginal());
			if (para->getBool("Show Original Sphere"))
				glDrawer.draw(GLDrawer::SPHERE, dataMgr.getCurrentOriginal());	
 		}		
 	}

	if ("Show Correspondences")
	{
		drawCorrespondences();
	}
	if (para->getBool("Show Skeleton"))
	{
		glDrawer.drawCurveSkeleton(*dataMgr.getCurrentSkeleton());
	}

//   	if ((para->getBool("Show Radius") || para->getBool("Show Bounding Box") ))
//   	{
//   		Box3f box = dataMgr.getCurrentSamples()->bbox;
//   		glBoxWire(box);
//   
//   		CoordinateFrame(box.Diag() / 2.0).Render(this, NULL);
//   	}


		if (show_something )
		{
			cout << "show something: Pick Dual Point 6" << something++ << endl;
		}


	if (!(takeSnapTile && para->getBool("No Snap Radius")))
	{
    if (para->getBool("Pick Dual Point"))
    {
      glDrawer.drawPickPoint(dataMgr.getCurrentSamples(), pickList, para->getBool("Show Samples Dot"));
      glDrawer.drawPickPoint(dataMgr.getCurrentDualSamples(), pickList, para->getBool("Show Samples Dot"));
    }
    else
    {
      glDrawer.drawPickPoint(dataMgr.getCurrentSamples(), pickList, para->getBool("Show Samples Dot"));
      //glDrawer.drawPickPoint(dataMgr.getCurrentDualSamples(), pickList, para->getBool("Show Samples Dot"));
    }


    //glDrawer.drawPickPoint(dataMgr.getCurrentSamples(), pickList, false);
    //glDrawer.drawPickPoint(dataMgr.getCurrentDualSamples(), pickList, false);


		if (global_paraMgr.drawer.getBool("Draw Picked Point Neighbor") && para->getBool("Show Radius"))
     {
      glDrawer.drawPickedPointNeighbor(dataMgr.getCurrentSamples(), pickList);
      glDrawer.drawPickedPointNeighbor(dataMgr.getCurrentDualSamples(), pickList);
       
			//glDrawer.drawPickedPointOriginalNeighbor(dataMgr.getCurrentDualSamples(), dataMgr.getCurrentOriginal(), pickList);
			
			
			
// 			if (para->getBool("Show Samples"))
// 			{
// 				glDrawer.drawPickedPointOriginalNeighbor(dataMgr.getCurrentDualSamples(), dataMgr.getCurrentSamples(), pickList);
// 			}

			//glDrawer.drawPickedPointOriginalNeighbor(dataMgr.getCurrentSamples(), dataMgr.getCurrentOriginal(), pickList);

     }


//     if (!pickList.empty())
//     {
//       glDrawer.drawPickedDisk(dataMgr.getCurrentDualSamples(), &picked_disk);
//     }

	}


	if (isDragging && para->getBool("Multiply Pick Point"))
	{
		drawPickRect();
	}
  lightOnOff(para->getBool("Light On or Off"));
  /* The following are semitransparent, careful for the rendering order*/
  glDepthMask(GL_FALSE);

  if (para->getBool("Show Radius")&& !(takeSnapTile && para->getBool("No Snap Radius"))) 
  {
		if (!para->getBool("Pick Dual Point"))
		{
			drawNeighborhoodRadius();
		}
    
  }

	if (show_something)
	{
		cout << "show something: final 7 " << something++ << endl;

	}


  glDepthMask(GL_TRUE);
  lightOnOff(para->getBool("Light On or Off"));

	if (show_something)
	{
		cout << "show something: final" << something++ << endl;
		show_something = false;
		something = 0;
	}

	glPopMatrix();
	glPopMatrix();

	if (takeSnapTile)
  {
		cout << "snap shot!" << endl;
		pasteTile();

		double normal_width = global_paraMgr.drawer.getDouble("Normal Line Width");
		double dot_size = global_paraMgr.drawer.getDouble("Sample Dot Size");
		double original_dot_size = global_paraMgr.drawer.getDouble("Original Dot Size");
		double dual_dot_size = global_paraMgr.drawer.getDouble("Dual Sample Dot Size");

		global_paraMgr.drawer.setValue("Normal Line Width", DoubleValue(normal_width / (SnapResolutionScale * SnapResolutionScale * snapDrawScal)));
		global_paraMgr.drawer.setValue("Sample Dot Size", DoubleValue(dot_size / (SnapResolutionScale * SnapResolutionScale * snapDrawScal)));
		global_paraMgr.drawer.setValue("Original Dot Size", DoubleValue(original_dot_size / (SnapResolutionScale * SnapResolutionScale * snapDrawScal)));

		global_paraMgr.drawer.setValue("Dual Sample Dot Size", DoubleValue(dual_dot_size / (SnapResolutionScale * SnapResolutionScale * snapDrawScal)));

	}

	}



PAINT_RETURN:
	paintMutex.unlock();
}


void GLArea::lightOnOff(bool _val)
{
	if(_val)
	{
		glEnable(GL_LIGHTING);
	}
	else
	{
		glDisable(GL_LIGHTING);
	}
}

void GLArea::initAfterOpenFile()
{
	//dataMgr.downSamplesByNum();
	dataMgr.recomputeQuad();
	//initView();
  dataMgr.getInitRadiuse();
	initSetting();
	//emit needUpdateStatus();
}

void GLArea::initSetting()
{
  dataMgr.recomputeQuad();
  initView();
	wlop.setFirstIterate();
	skeletonization.setFirstIterate();
	emit needUpdateStatus();
}

void GLArea::initView()
{
	dataMgr.recomputeBox();
	if (!dataMgr.isOriginalEmpty())
	{
		gl_box = dataMgr.getCurrentOriginal()->bbox;
	}
	else if(!dataMgr.isSamplesEmpty())
	{
		gl_box = dataMgr.getCurrentSamples()->bbox;
	}
}

void GLArea::runPointCloudAlgorithm(PointCloudAlgorithm& algorithm)
{
	paintMutex.lock();

	QString name = algorithm.getParameterSet()->getString("Algorithm Name");
	cout << "*********************************** Start  " << name.toStdString() << "  ***********************************" << endl;
	int starttime, stoptime, timeused;
	starttime = clock();

	algorithm.setInput(&dataMgr);
	algorithm.run();
	algorithm.clear();

	stoptime = clock();
	timeused = stoptime - starttime;

	int currentUsedTime = timeused/double(CLOCKS_PER_SEC);

	cout << "time used:  " << timeused/double(CLOCKS_PER_SEC) << " seconds." << endl;
	cout << "*********************************** End  " << name.toStdString() << "  ***********************************" << endl;
	cout << endl << endl << endl;

	paintMutex.unlock();
}


void GLArea::openByDrop(QString fileName)
{
	if(fileName.endsWith("ply"))
	{
		if (fileName.contains("original"))
		{
			dataMgr.loadPlyToOriginal(fileName);
		}
		else if (fileName.contains("dual"))
		{
			dataMgr.loadPlyToDualSample(fileName);
		}
		else
		{
			dataMgr.loadPlyToSample(fileName);
		}
	}
	if(fileName.endsWith("View"))
	{
		loadView(fileName);
	}
	if (fileName.endsWith("VPoint"))
	{
		loadVPoint(fileName);
	}
	if (fileName.endsWith("paras"))
	{
		ifstream inpara(fileName.toStdString().c_str());
		global_paraMgr.inputAllParameters(inpara);
		inpara.close();
	}
	if (fileName.endsWith("skel"))
	{
		if (fileName.contains("target"))
		{
			dataMgr.loadTargetSkeletonFromSkel(fileName);
		}
		else
		{
			dataMgr.loadSkeletonFromSkel(fileName);
		}
	}
	else if(fileName.endsWith("RGBN"))
	{
		readRGBNormal(fileName);
	}
	if (fileName.endsWith("jpg"))
	{
		dataMgr.loadImage(fileName);
	}
  if(fileName.endsWith("xyz"))
  {
    dataMgr.loadXYZN(fileName);
  }

	show_something = true;
	cout << "finish open " << endl; 
	initAfterOpenFile();
	updateGL();
	emit needUpdateStatus();
}

void GLArea::loadDefaultModel()
{
	dataMgr.loadDefaultSphere();
	dataMgr.loadPlyToSample("default.ply");
	dataMgr.loadPlyToOriginal("default_original.ply");
	initAfterOpenFile();
	updateGL();
}

void GLArea::drawPickRect()
{
	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0,width(),height(),0,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_COLOR_LOGIC_OP);
	glLogicOp(GL_XOR);
	glColor3f(1,1,1);

	glBegin(GL_LINE_LOOP);
	glVertex2f(x1,viewport[3] - y1);
	glVertex2f(x2,viewport[3] - y1);
	glVertex2f(x2,viewport[3] - y2);
	glVertex2f(x1,viewport[3] - y2);
	glEnd();
	glDisable(GL_LOGIC_OP);

	// Closing 2D
	glPopAttrib();
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}


void GLArea::drawNeighborhoodRadius()
{
  if (dataMgr.getCurrentSamples()->vert.empty()) return;

  Point3f p;
  if(!pickList.empty() && pickList[0] >= 0)
  {
    int id = pickList[0];
    if (id >= 0 && id < dataMgr.getCurrentSamples()->vert.size())
    {
      p = dataMgr.getCurrentSamples()->vert[id].P();
    }
    else
    {
      p = dataMgr.getCurrentSamples()->vert[0].P();
    }
  }
  else
  {
    p = dataMgr.getCurrentSamples()->vert[0].P();
  }

  double h_Gaussian_para = global_paraMgr.wLop.getDouble("H Gaussian Para");
  double grid_radius = global_paraMgr.data.getDouble("CGrid Radius");

  if (!takeSnapTile && para->getBool("Show Red Radius Line"))
  {
    glColor3f(1, 0, 0);
    glLineWidth(3);

    glBegin(GL_LINES);

    glVertex3f(p[0], p[1], p[2]);
    glVertex3f(p[0], p[1] + grid_radius, p[2]);
    glEnd();
  }

  glMatrixMode(GL_MODELVIEW_MATRIX);

  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  double trans_value = para->getDouble("Radius Ball Transparency");
  glColor4f(0,0,1,trans_value);
  glShadeModel(GL_SMOOTH);


  CMesh* samples = dataMgr.getCurrentSamples();

	CMesh* ellipsoid = dataMgr.getCurrentEllipsoid();
// 	if (para->getBool("Show Eigen Directions") && !ellipsoid->vert.empty())
// 	{
// 	}
// 	else 
		if (para->getBool("Show All Radius") && samples->vn < 1000)
  {
    for(int i = 0; i < samples->vert.size(); i++)
    {
      CVertex& v = samples->vert[i];
      glPushMatrix();
      glTranslatef(v.P()[0], v.P()[1], v.P()[2]);
      glutSolidSphere(grid_radius  / sqrt(h_Gaussian_para), 40, 40);
      glPopMatrix();
    }
  }
  else
  {
    glPushMatrix();
    glTranslatef(p[0], p[1], p[2]);
    glutSolidSphere(grid_radius  / sqrt(h_Gaussian_para), 40, 40);

    if (para->getBool("Show Red Radius Line") 
      && para->getBool("Show Skeleton")
      && !dataMgr.getCurrentSkeleton()->isEmpty())
    {
      double branch_merge_radius = global_paraMgr.skeleton.getDouble("Branches Merge Max Dist");
      glColor4f(0,1,0.5,0.4);
      glutSolidSphere(branch_merge_radius, 40, 40);

    }
    glPopMatrix();
  }

  if (para->getBool("Show Samples Dot"))
  {
    glDepthMask(GL_TRUE);
  }

  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);
  glDisable(GL_BLEND);
  glDisable(GL_CULL_FACE); 
}

void GLArea::drawLightBall()
{
	// ============== LIGHT TRACKBALL ==============
	// Apply the trackball for the light direction

	glPushMatrix();
	trackball_light.GetView();
	trackball_light.Apply(!(isDefaultTrackBall()));

	static float lightPosF[]={0.0,0.0,1.0,0.0};
	glLightfv(GL_LIGHT0,GL_POSITION,lightPosF);
	static float lightPosB[]={0.0,0.0,-1.0,0.0};
	glLightfv(GL_LIGHT1,GL_POSITION,lightPosB);

	if (!(isDefaultTrackBall()))
	{
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
		glColor3f(1,1,0);
		glDisable(GL_LIGHTING);
		const unsigned int lineNum=3;
		glBegin(GL_LINES);
		for(unsigned int i=0;i<=lineNum;++i)
			for(unsigned int j=0;j<=lineNum;++j) {
				glVertex3f(-1.0f+i*2.0/lineNum,-1.0f+j*2.0/lineNum,-2);
				glVertex3f(-1.0f+i*2.0/lineNum,-1.0f+j*2.0/lineNum, 2);
			}
			glEnd();
			glPopAttrib();
	}
	glPopMatrix();
}

void GLArea::changeColor(QString paraName)
{
	QColor qcolor;
  if (global_paraMgr.drawer.hasParameter(paraName))
  {
    qcolor = global_paraMgr.drawer.getColor(paraName);
  }
  else if (global_paraMgr.glarea.hasParameter(paraName))
  {
    qcolor = global_paraMgr.glarea.getColor(paraName);
  }
  else
  {
    return;
  }
	

	qcolor = QColorDialog::getColor(qcolor);

	if(qcolor.isValid()){

		if(paraName.contains("Light"))
		{
			GLColor color(qcolor);
			float light_col[4] = {color.r, color.g, color.b, 1.0};

			if(paraName.contains("Ambient"))
			{
				glLightfv(GL_LIGHT0, GL_AMBIENT, light_col);
				glLightfv(GL_LIGHT1, GL_AMBIENT, light_col);

			}
			else if(paraName.contains("Diffuse"))
			{
				glLightfv(GL_LIGHT0, GL_DIFFUSE, light_col);
				glLightfv(GL_LIGHT1, GL_DIFFUSE, light_col);

			}
			else if(paraName.contains("Specular"))
			{
				glLightfv(GL_LIGHT0, GL_SPECULAR, light_col);
				glLightfv(GL_LIGHT1, GL_SPECULAR, light_col);
			}
			para->setValue(paraName, ColorValue(qcolor));
		}


    if (global_paraMgr.drawer.hasParameter(paraName))
    {
      global_paraMgr.drawer.setValue(paraName, ColorValue(qcolor));
    }
    else if (global_paraMgr.glarea.hasParameter(paraName))
    {
      global_paraMgr.glarea.setValue(paraName, ColorValue(qcolor));
    }
		

	} 
}

int GLArea::pickPoint(int x, int y, vector<int> &result, int width, int height, bool only_one)
{
	if(dataMgr.isSamplesEmpty())
		return -1;

	result.clear();

	if(width==0 ||height==0) return 0; 
	long hits;	
	CMesh* samples = dataMgr.getCurrentSamples();

	if (global_paraMgr.drawer.getBool("Use Pick Original"))
	{
		samples = dataMgr.getCurrentOriginal();
	}

	int sz= samples->vert.size();

	//if (global_paraMgr.drawer.getBool("Use Pick Skeleton"))
	//{
	//	sz = dataMgr.getCurrentSkeleton()->size;
	//}
	GLuint *selectBuf =new GLuint[sz*5];

	//  static unsigned int selectBuf[16384];
	glSelectBuffer(sz * 5, selectBuf);
	glRenderMode(GL_SELECT);
	glInitNames();

	/* Because LoadName() won't work with no names on the stack */
	glPushName(-1);

	double mp[16];

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
	glMatrixMode(GL_PROJECTION);
	glGetDoublev(GL_PROJECTION_MATRIX ,mp);
	glPushMatrix();
	glLoadIdentity();
	gluPickMatrix(x, y, width, height, viewport);
	glMultMatrixd(mp);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	//doPick = true;
	global_paraMgr.drawer.setValue("Doing Pick", BoolValue(true));


	//if (global_paraMgr.drawer.getBool("Use Pick Skeleton"))
	//{
	//	just_draw_skel_for_pick = true;
	//}
	paintGL();
	//just_draw_skel_for_pick = false;


	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	hits = glRenderMode(GL_RENDER);

	//xstring buf;
	if (hits <= 0)     return 0;

	vector<int> H;

	if (hits > 1 && global_paraMgr.drawer.getBool("Use Pick Mode2"))
	{
		hits--;
	}

	for(long ii=0;ii<hits;ii++)
	{
		//H.push_back( std::pair<double,unsigned int>(selectBuf[ii*4+1]/4294967295.0,selectBuf[ii*4+3]));
		H.push_back(selectBuf[ii*4+3]);
	}
	sort(H.begin(),H.end());


	//Only One Pick
	for(long i = 0 ;i < H.size();i++)
	{
		if(H[i] >= 0 && H[i] < sz)
		{
			result.push_back(H[i]);
			if(only_one)
				break;
		}
	}

	delete [] selectBuf;

	global_paraMgr.drawer.setValue("Doing Pick", BoolValue(false));

	if(only_one)
	{
		if(result.empty())
			return -1;
		else
		{
			return result[0];
		}
	}

	pickList.pop_back();

	return result.size();
}


void GLArea::setView()
{
	glViewport(0,0, this->width(),this->height());
	curSiz.setWidth(this->width());
	curSiz.setHeight(this->height());

	GLfloat fAspect = (GLfloat)curSiz.width()/ curSiz.height();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// This parameter is the one that controls:
	// HOW LARGE IS THE TRACKBALL ICON ON THE SCREEN.
	float viewRatio = 1.75f;
	float cameraDist = viewRatio / tanf(math::ToRad(fov*.5f));

	if(fov==5)
		cameraDist = 1000; // small hack for orthographic projection where camera distance is rather meaningless...
	nearPlane = cameraDist - 2.f*clipRatioNear;
	farPlane =  cameraDist + 10.f*clipRatioFar;
	if(nearPlane<=cameraDist*.1f) nearPlane=cameraDist*.1f;

	if (!takeSnapTile)
	{
		if(fov==5)	glOrtho( -viewRatio*fAspect, viewRatio*fAspect, -viewRatio, viewRatio, cameraDist - 2.f*clipRatioNear, cameraDist+2.f*clipRatioFar);
		else    		gluPerspective(fov, fAspect, nearPlane, farPlane);
	}
	else	setTiledView(fov, viewRatio, fAspect, nearPlane, farPlane, cameraDist);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, -cameraDist,0, 0, 0, 0, 1, 0);

}

void GLArea::pasteTile()
{
	glPushAttrib(GL_ENABLE_BIT);
	QImage tileBuffer=grabFrameBuffer(true).mirrored(false,true);

	if (snapBuffer.isNull())
		snapBuffer = QImage(tileBuffer.width() * ss.resolution, tileBuffer.height() * ss.resolution, tileBuffer.format());

	uchar *snapPtr = snapBuffer.bits() + (tileBuffer.bytesPerLine() * tileCol) + ((totalCols * tileRow) * tileBuffer.byteCount());
	uchar *tilePtr = tileBuffer.bits();

	for (int y=0; y < tileBuffer.height(); y++)
	{
		memcpy((void*) snapPtr, (void*) tilePtr, tileBuffer.bytesPerLine());
		snapPtr+=tileBuffer.bytesPerLine() * totalCols;
		tilePtr+=tileBuffer.bytesPerLine();
	}

	tileCol++;

	if (tileCol >= totalCols)
	{
		tileCol=0;
		tileRow++;

		if (tileRow >= totalRows)
		{
			QString outfile=QString("%1/%2%3.png")
				.arg(ss.outdir)
				.arg(ss.basename)
				.arg("");
			//.arg(ss.counter++,2,10,QChar('0'));
			bool ret = (snapBuffer.mirrored(false,true)).save(outfile,"PNG");

			takeSnapTile=false;
			recoverView();

			snapBuffer=QImage();
		}
	}
	update();
	glPopAttrib();
}

void GLArea::recoverView()
{
	int w = width();
	int h = height();

	glViewport(0, 0, (GLint)w, (GLint)h);  
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();  

	float r = w/(float)h;
	gluPerspective(60, r, 0.1, 10);
	glMatrixMode(GL_MODELVIEW);  
}

void GLArea::setTiledView(GLdouble fovY, float viewRatio, float fAspect, GLdouble zNear, GLdouble zFar,  float cameraDist)
{
	if(fovY<=5)
	{
		GLdouble fLeft   = -viewRatio*fAspect;
		GLdouble fRight  =  viewRatio*fAspect;
		GLdouble fBottom = -viewRatio;
		GLdouble fTop    =  viewRatio;

		GLdouble tDimX = fabs(fRight-fLeft) / totalCols;
		GLdouble tDimY = fabs(fTop-fBottom) / totalRows;


		glOrtho(fLeft   + tDimX * tileCol, fLeft   + tDimX * (tileCol+1),     /* left, right */
			fBottom + tDimY * tileRow, fBottom + tDimY * (tileRow+1),     /* bottom, top */
			cameraDist - 2.f*clipRatioNear, cameraDist+2.f*clipRatioFar);
	}
	else
	{
		GLdouble fTop    = zNear * tan(math::ToRad(fovY/2.0));
		GLdouble fRight  = fTop * fAspect;
		GLdouble fBottom = -fTop;
		GLdouble fLeft   = -fRight;

		// tile Dimension
		GLdouble tDimX = fabs(fRight-fLeft) / totalCols;
		GLdouble tDimY = fabs(fTop-fBottom) / totalRows;

		glFrustum(fLeft   + tDimX * tileCol, fLeft   + tDimX * (tileCol+1),
			fBottom + tDimY * tileRow, fBottom + tDimY * (tileRow+1), zNear, zFar);
	}
}

void GLArea::saveSnapshot()
{
	is_paintGL_locked = true;

	double SnapResolutionScale = para->getDouble("Snapshot Resolution");
	ss.resolution = SnapResolutionScale;
	snapDrawScal = 1. / SnapResolutionScale;
	
	totalCols = totalRows = ss.resolution;
	tileRow = tileCol = 0;
	ss.setOutDir(current_snap_path);
	ss.GetSysTime();

	if (para->getBool("SnapShot Each Iteration"))
	{
		int snape_idx = para->getDouble("Snapshot Index");

		ss.basename = QString::number(snape_idx++,10); 
		double snap_id = para->getDouble("Snapshot Index");
		global_paraMgr.setGlobalParameter("Snapshot Index", DoubleValue(snap_id));
		emit needUpdateStatus();
	
		para->setValue("Snapshot Index", DoubleValue(snape_idx));

		if (para->getBool("Need Snap Files"))
		{
			QString skel_file_name = QString(QString(".\\snapfile\\")) + ss.basename + QString(".skel");
			dataMgr.saveSkeletonAsSkel(skel_file_name);

// 			skel_file_name.replace(".skel", ".View");
// 			saveView(skel_file_name);

			skel_file_name.replace(".skel", ".paras");
			ifstream inpara(skel_file_name.toStdString().c_str());
			global_paraMgr.inputAllParameters(inpara);
			inpara.close();
		}

	}

	takeSnapTile = true;
	//updateGL();
	if (SnapResolutionScale > 1)
	{
		for (int i = 0; i < 4; i++)
		{
			is_paintGL_locked = false;
			updateGL();
			is_paintGL_locked = true;
		}
	}
	else
	{
		is_paintGL_locked = false;
		updateGL();
		is_paintGL_locked = true;
	}

	is_paintGL_locked = false;
}

void GLArea::runWlop()
{
	if (dataMgr.isOriginalEmpty())
	{
		dataMgr.subSamples();
	}

  global_paraMgr.glarea.setValue("GLarea Busying", BoolValue(true));

	cout << "wloppppppppppppppppppppppppppppp 1" << endl;

  bool is_break = false;
	for (int i = 0; i < global_paraMgr.wLop.getDouble("Num Of Iterate Time"); i++)
	{
    if (global_paraMgr.glarea.getBool("Algorithom Stop") || is_break)
    {
      is_break = true;
      break;
    }


		
		runPointCloudAlgorithm(wlop);



// 		int sleep_time = 100;
// 		QTime dieTime = QTime::currentTime().addMSecs(sleep_time);
// 		while (QTime::currentTime() < dieTime)
// 			QCoreApplication::processEvents(QEventLoop::AllEvents, sleep_time);

		updateUI();
		update();
		updateGL();
    emit needUpdateStatus();



// 		int knn = global_paraMgr.norSmooth.getInt("PCA KNN");
// 		CMesh* samples;
// 		samples = dataMgr.getCurrentSamples();
// 		vector<Point3f> remember_normal(samples->vert.size());
// 		for (int i = 0; i < samples->vert.size(); i++)
// 		{
// 			remember_normal[i] = samples->vert[i].N();
// 		}
// 		vcg::tri::PointCloudNormal<CMesh>::Param pca_para;
// 		pca_para.fittingAdjNum = knn;
// 		vcg::tri::PointCloudNormal<CMesh>::Compute(*samples, pca_para, NULL);
// 		for (int i = 0; i < samples->vert.size(); i++)
// 		{
// 			CVertex& v = samples->vert[i];
// 			if (v.N() * remember_normal[i] < 0)
// 			{
// 				v.N() *= -1;
// 			}
// 		}
	}


	if (para->getBool("SnapShot Each Iteration"))
	{
		if (!global_paraMgr.wLop.getBool("Run Regularize Normals")
			&& !global_paraMgr.wLop.getBool("Run Compute Confidence")
			&& !global_paraMgr.wLop.getBool("Run Move Backward"))
		{
			saveSnapshot();

		}

	}

// 	CMesh* samples = dataMgr.getCurrentSamples();
// 
// 	for (int i = 0; i < samples->vert.size(); i++)
// 	{
// 		CVertex& v = samples->vert[i];
// 
// 		if (i < 50)
// 		{
// 			cout << "e e: " << v.eigen_confidence << endl;
// 		}
// 	}

	cout << "wloppppppppppppppppppppppppppppp 2" << endl;

  if (is_break)
  {
    global_paraMgr.skeleton.setValue("The Skeletonlization Process Should Stop", BoolValue(false));
    global_paraMgr.glarea.setValue("Algorithom Stop", BoolValue(false));
    global_paraMgr.glarea.setValue("GLarea Busying", BoolValue(false));
    return;
  }

	para->setValue("Running Algorithm Name",
		StringValue(wlop.getParameterSet()->getString("Algorithm Name")));

	//if (para->getBool("SnapShot Each Iteration"))
	//{
	//	saveSnapshot();
	//}

	global_paraMgr.wLop.setValue("Run Anisotropic LOP", BoolValue(false));
	global_paraMgr.wLop.setValue("Run Skel WLOP", BoolValue(false));
	

	global_paraMgr.wLop.setValue("Dual Samples Represent Inner Points", BoolValue(false));
	global_paraMgr.wLop.setValue("Dual Samples Represent Skeltal Points", BoolValue(false));

	cout << "wloppppppppppppppppppppppppppppp 3" << endl;

	//show_something = true;
}

void GLArea::runSkeletonization_linear()
{
  if (dataMgr.isSamplesEmpty())
  {
    return;
  }

  runPointCloudAlgorithm(skeletonization);

  para->setValue("Running Algorithm Name",
    StringValue(skeletonization.getParameterSet()->getString("Algorithm Name")));

  emit needUpdateStatus();
}

void GLArea::runSkeletonization_paralleled()
{
	if (dataMgr.isOriginalEmpty() || dataMgr.isSamplesEmpty())
	{
		return;
	}
  global_paraMgr.glarea.setValue("GLarea Busying", BoolValue(true));

  global_paraMgr.skeleton.setValue("The Skeletonlization Process Should Stop", BoolValue(false));
  double current_radius = global_paraMgr.skeleton.getDouble("CGrid Radius");
  global_paraMgr.skeleton.setValue("Initial Radius", DoubleValue(current_radius));


  bool is_break = false;
  int MAX_SKELETON_ITERATE = 500;
  //if (global_paraMgr.skeleton.getBool("Run Auto Wlop One Step"))
  //{
  //  MAX_SKELETON_ITERATE = 1;
  //}

  for (int i = 0; i < MAX_SKELETON_ITERATE; i++)
  {
    if (global_paraMgr.glarea.getBool("Algorithom Stop") || is_break)
    {
      is_break = true;
      break;
    }
    //paintMutex.lock();
    global_paraMgr.skeleton.setValue("Run Auto Wlop One Step", BoolValue(true));
    runPointCloudAlgorithm(skeletonization);
    global_paraMgr.skeleton.setValue("Run Auto Wlop One Step", BoolValue(true));
   // paintMutex.unlock();

    if (global_paraMgr.skeleton.getBool("The Skeletonlization Process Should Stop"))
    {
      break;
    }
    emit needUpdateStatus();
  }
  global_paraMgr.glarea.setValue("GLarea Busying", BoolValue(false));
  
  if (is_break)
  {
    global_paraMgr.skeleton.setValue("The Skeletonlization Process Should Stop", BoolValue(false));
    global_paraMgr.glarea.setValue("Algorithom Stop", BoolValue(false));
    global_paraMgr.glarea.setValue("GLarea Busying", BoolValue(false));
    return;
  }

  if (global_paraMgr.skeleton.getBool("Run Auto Wlop One Stage"))
  {
    global_paraMgr.skeleton.setValue("Run Auto Wlop One Stage", BoolValue(false));
    emit needUpdateStatus();
    return;
  }

	para->setValue("Running Algorithm Name",
		StringValue(skeletonization.getParameterSet()->getString("Algorithm Name")));
	
	
}


void GLArea::runUpsampling()
{
	if (dataMgr.isSamplesEmpty())
	{
		return;
	}

	runPointCloudAlgorithm(upsampler);

	para->setValue("Running Algorithm Name",
		StringValue(upsampler.getParameterSet()->getString("Algorithm Name")));

	emit needUpdateStatus();
}


void GLArea::runNormalSmoothing()
{
	if (dataMgr.isSamplesEmpty())
	{
		return;
	}

	for (int i = 0; i < global_paraMgr.norSmooth.getInt("Number Of Iterate"); i++)
	{
		runPointCloudAlgorithm(norSmoother);
	}

	emit needUpdateStatus();
}

void GLArea::outputColor(ostream& out, QColor& color)
{
	out << color.red() << "	" << color.green() << "	" << color.blue() << endl;
}

QColor GLArea::inputColor(istream& in)
{
	QColor color;
	int r, g, b;
	in >> r >> g >> b;
	color.setRgb(r, g, b);
	return color;
}

void GLArea::saveVPoint(QString fileName)
{
	QString file = fileName;
	ofstream outfile(file.toStdString().c_str());

	char viewStr[100];
	trackball.ToAscii(viewStr);
	cout << "saveView" << viewStr << endl;
	outfile << viewStr << endl;

	trackball_light.ToAscii(viewStr);
	outfile << viewStr << endl;

	outfile << rotate_pos[0] << " " << rotate_pos[1] << " " << rotate_pos[2] << " " << endl;
	outfile << rotate_normal[0] << " " << rotate_normal[1] << " " << rotate_normal[2] << " " << endl;

	outfile.close();
}

void GLArea::saveView(QString fileName)
{
	QString file = fileName;
	ofstream outfile(file.toStdString().c_str());

	char viewStr[100];
	trackball.ToAscii(viewStr);
	cout << "saveView" << viewStr << endl;
	outfile << viewStr << endl;

	CMesh* samples = dataMgr.getCurrentSamples();
	//Point3d cut_pos = samples->vert[0].P();
	

	outfile << global_paraMgr.wLop.getDouble("CGrid Radius") << endl;
	cout << "save radius:" << global_paraMgr.wLop.getDouble("CGrid Radius") << endl;

	//outfile << cut_pos[0] << " " << cut_pos[1] << " " << cut_pos[2] << endl;	
	outfile << -1 << " " << -1 << " " << -1 << endl;	


	outfile << global_paraMgr.drawer.getDouble("Sample Dot Size") << endl;
	outfile << global_paraMgr.drawer.getDouble("Sample Draw Width") << endl;
	outfile << global_paraMgr.drawer.getDouble("Normal Line Width") << endl;
	outfile << global_paraMgr.drawer.getDouble("Normal Line Length") << endl;

 	outfile << rotate_pos[0] << " " << rotate_pos[1] << " " << rotate_pos[2] << " " << endl;
 	outfile << rotate_normal[0] << " " << rotate_normal[1] << " " << rotate_normal[2] << " " << endl;
//	outfile << 0.0 << " " << 0.0 << " " << rotate_pos[2] << " " << endl;
//	outfile << 0.0 << " " << 0.0 << " " << 1.0  << " " << endl;

	outfile << global_paraMgr.drawer.getDouble("Original Dot Size") << endl;

	outfile << global_paraMgr.drawer.getDouble("Skeleton Bone Width") << endl;
	outfile << global_paraMgr.drawer.getDouble("Skeleton Node Size") << endl;
	outfile << global_paraMgr.drawer.getDouble("Skeleton Branch Size") << endl;

	outputColor(outfile, global_paraMgr.drawer.getColor("Sample Point Color"));
	outputColor(outfile, global_paraMgr.drawer.getColor("Original Point Color"));
	outputColor(outfile, global_paraMgr.glarea.getColor("Light Ambient Color"));
	outputColor(outfile, global_paraMgr.glarea.getColor("Light Diffuse Color"));
	outputColor(outfile, global_paraMgr.glarea.getColor("Light Specular Color"));
	outputColor(outfile, global_paraMgr.drawer.getColor("Skeleton Bone Color"));
	outputColor(outfile, global_paraMgr.drawer.getColor("Skeleton Node Color"));
	outputColor(outfile, global_paraMgr.drawer.getColor("Skeleton Branch Color"));
	outputColor(outfile, global_paraMgr.drawer.getColor("Normal Line Color"));


	outfile << global_paraMgr.wLop.getDouble("CGrid Radius") << endl;
	//outfile << global_paraMgr.wLop.getDouble("Local Density Radius") << endl;
	//outfile << global_paraMgr.skeleton.getDouble("Snake Search Max Dist Blue") << endl;
	//outfile << global_paraMgr.skeleton.getDouble("Branch Search Max Dist Yellow") << endl;
	//outfile << global_paraMgr.skeleton.getDouble("Branches Merge Max Dist Orange") << endl;

	outfile << -1<< "	" << -1 << "	"<< -1 << "	"<< -1 << "	"<< endl;

	outfile << global_paraMgr.wLop.getDouble("Repulsion Mu") << endl;
	outfile << global_paraMgr.wLop.getDouble("Repulsion Mu2") << endl;


	trackball_light.ToAscii(viewStr);
	outfile << viewStr << endl;


	double init_radius = global_paraMgr.skeleton.getDouble("Initial Radius") ;
	if (init_radius > 0)
	{
		outfile << global_paraMgr.skeleton.getDouble("Initial Radius") << endl;
	}
	else
	{
		outfile << global_paraMgr.skeleton.getDouble("CGrid Radius") << endl;
	}

	outfile << global_paraMgr.glarea.getDouble("Radius Ball Transparency") << endl;


	//new parameters 2015-1-11
	outputColor(outfile, global_paraMgr.drawer.getColor("DLink Color"));
	outputColor(outfile, global_paraMgr.drawer.getColor("Backface Color"));

	outfile << global_paraMgr.drawer.getDouble("Original Draw Width") << endl;
	outfile << global_paraMgr.drawer.getDouble("Dual Sample Draw Width") << endl;
	outfile << global_paraMgr.drawer.getDouble("Dual Sample Dot Size") << endl;
	
	
	outfile << global_paraMgr.wLop.getDouble("Original Averaging KNN") << endl;
	outfile << global_paraMgr.wLop.getDouble("Increasing Step Size") << endl;
	outfile << global_paraMgr.wLop.getDouble("Local Neighbor Size For Inner Points") << endl;
	outfile << global_paraMgr.wLop.getDouble("Local Neighbor Size For Surface Points") << endl;
	outfile << global_paraMgr.wLop.getDouble("Inner Points Cooling Parameter") << endl;
	outfile << global_paraMgr.wLop.getDouble("Local Angle Threshold") << endl;
	outfile << global_paraMgr.wLop.getDouble("Eigen Neighborhood Para1") << endl;
	outfile << global_paraMgr.wLop.getDouble("Eigen Neighborhood Para2") << endl;
	
	outfile << global_paraMgr.wLop.getDouble("Average Closest Dist") << endl;
	outfile << global_paraMgr.wLop.getDouble("Average Dist To Input Threshold") << endl;
	
	outfile << global_paraMgr.wLop.getDouble("Original Confidence KNN") << endl;
	outfile << global_paraMgr.wLop.getDouble("Similarity Term Neighbor Para") << endl;
	outfile << global_paraMgr.wLop.getDouble("Density Confidence Segment Threshold") << endl;
	outfile << global_paraMgr.wLop.getDouble("Eigen Directional Threshold") << endl;
	outfile << global_paraMgr.wLop.getDouble("Save Move Dist Along Normal Para") << endl;
	outfile << global_paraMgr.wLop.getDouble("Big Repulsion Power") << endl;
	outfile << global_paraMgr.wLop.getDouble("Protect Small Tubular Para") << endl;
	outfile << global_paraMgr.wLop.getDouble("Protect High Confidence Para") << endl;

	outfile << global_paraMgr.norSmooth.getInt("PCA KNN") << endl;
	outfile << global_paraMgr.norSmooth.getDouble("Sharpe Feature Bandwidth Sigma") << endl;
	outfile << global_paraMgr.upsampling.getDouble("Feature Sigma") << endl;
	outfile << global_paraMgr.upsampling.getDouble("Dist Threshold") << endl;
	outfile << global_paraMgr.skeleton.getDouble("Sigma KNN") << endl;
	outfile << global_paraMgr.skeleton.getDouble("Eigen Feature Identification Threshold") << endl;

	outfile << global_paraMgr.wLop.getDouble("Data Outweigh Similarity Para") << endl;
	outfile << global_paraMgr.wLop.getDouble("Confidence Power") << endl;
	outfile << global_paraMgr.wLop.getDouble("KNN For Similarity") << endl;

	outfile << global_paraMgr.glarea.getDouble("Sample Confidence Color Scale") << endl;
	outfile << global_paraMgr.glarea.getDouble("Point ISO Value Shift") << endl;

	outfile.close();
}


void GLArea::loadVPoint(QString fileName)
{
	ifstream infile(fileName.toStdString().c_str());
	double temp;

	string viewStr;
	infile >> viewStr;

	trackball.SetFromAscii(viewStr.c_str());

	infile >> viewStr;
	trackball_light.SetFromAscii(viewStr.c_str());

	initLight();

	if (!infile.eof())
	{
		Point3f r_pos, r_normal;
		infile >> r_pos[0] >> r_pos[1] >> r_pos[2];
		infile >> r_normal[0] >> r_normal[1] >> r_normal[2];

		rotate_pos = r_pos;
		rotate_normal = r_normal;
	}



	infile.close();
	emit needUpdateStatus();
}

void GLArea::loadView(QString fileName)
{
	ifstream infile(fileName.toStdString().c_str());
	double temp;

	string viewStr;
	infile >> viewStr;

	trackball.SetFromAscii(viewStr.c_str());

	double radius;
	infile >> radius;
	global_paraMgr.setGlobalParameter("CGrid Radius", DoubleValue(radius));
	cout << "load radius: " << radius << endl;

	infile >> temp >> temp >> temp;

	double dot_size;
	double quade_width;
	double normal_width;
	double normal_length;

	infile >> dot_size >> quade_width >> normal_width >> normal_length;
	global_paraMgr.drawer.setValue("Sample Dot Size", DoubleValue(dot_size));
	global_paraMgr.drawer.setValue("Sample Draw Width", DoubleValue(quade_width));
	global_paraMgr.drawer.setValue("Normal Line Width", DoubleValue(normal_width));
	global_paraMgr.drawer.setValue("Normal Line Length", DoubleValue(normal_length));

	Point3f r_pos, r_normal;
	infile >> r_pos[0] >> r_pos[1] >> r_pos[2];
	infile >> r_normal[0] >> r_normal[1] >> r_normal[2];

	rotate_pos = r_pos;
	rotate_normal = r_normal;

	//double temp;
	QColor c;
	if (!infile.eof())
	{
		infile >> dot_size;
		global_paraMgr.drawer.setValue("Original Dot Size", DoubleValue(dot_size));

		infile >> temp;
		global_paraMgr.drawer.setValue("Skeleton Bone Width", DoubleValue(temp));

		infile >> temp;
		global_paraMgr.drawer.setValue("Skeleton Node Size", DoubleValue(temp));

		infile >> temp;
		global_paraMgr.drawer.setValue("Skeleton Branch Size", DoubleValue(temp));

		global_paraMgr.drawer.setValue("Sample Point Color", ColorValue(inputColor(infile)));
		global_paraMgr.drawer.setValue("Original Point Color", ColorValue(inputColor(infile)));
		global_paraMgr.glarea.setValue("Light Ambient Color", ColorValue(inputColor(infile)));
		global_paraMgr.glarea.setValue("Light Diffuse Color", ColorValue(inputColor(infile)));
		global_paraMgr.glarea.setValue("Light Specular Color", ColorValue(inputColor(infile)));
		global_paraMgr.drawer.setValue("Skeleton Bone Color", ColorValue(inputColor(infile)));
		global_paraMgr.drawer.setValue("Skeleton Node Color", ColorValue(inputColor(infile)));
		global_paraMgr.drawer.setValue("Skeleton Branch Color", ColorValue(inputColor(infile)));
		global_paraMgr.drawer.setValue("Normal Line Color", ColorValue(inputColor(infile)));
	}

	if (!infile.eof())
	{
		infile >> temp;
		//global_paraMgr.setGlobalParameter("CGrid Radius", DoubleValue(temp));

		infile >> temp;
		//global_paraMgr.setGlobalParameter("Local Density Radius", DoubleValue(temp));

	}

	if (!infile.eof())
	{
		infile >> temp;
		//global_paraMgr.skeleton.setValue("Snake Search Max Dist Blue", DoubleValue(temp));

		infile >> temp;
		//global_paraMgr.skeleton.setValue("Branch Search Max Dist Yellow", DoubleValue(temp));

		infile >> temp;
		//global_paraMgr.skeleton.setValue("Branches Merge Max Dist Orange", DoubleValue(temp));
	}

	if (!infile.eof())
	{
		infile >> temp;
		global_paraMgr.wLop.setValue("Repulsion Mu", DoubleValue(temp));

		infile >> temp;
		global_paraMgr.wLop.setValue("Repulsion Mu2", DoubleValue(temp));

		infile >> viewStr;
		trackball_light.SetFromAscii(viewStr.c_str());
	}
	initLight();

	if (!infile.eof())
	{
		double init_radius;
		infile >> init_radius;
// 		if (init_radius > 0)
// 		{
// 			global_paraMgr.setGlobalParameter("Initial Radius", DoubleValue(init_radius));
// 		}
// 		else
// 		{
// 			global_paraMgr.setGlobalParameter("Initial Radius", DoubleValue(global_paraMgr.wLop.getDouble("CGrid Radius")));
// 		}
	}

	infile >> temp;
	global_paraMgr.glarea.setValue("Radius Ball Transparency", DoubleValue(temp));

	if (!infile.eof())
	{
		//new parameters 2015-1-11
		global_paraMgr.drawer.setValue("DLink Color", ColorValue(inputColor(infile)));
		global_paraMgr.drawer.setValue("Backface Color", ColorValue(inputColor(infile)));


		infile >> temp; global_paraMgr.drawer.setValue("Original Draw Width", DoubleValue(temp));
		infile >> temp; global_paraMgr.drawer.setValue("Dual Sample Draw Width", DoubleValue(temp));
		infile >> temp; global_paraMgr.drawer.setValue("Dual Sample Dot Size", DoubleValue(temp));

		infile >> temp; global_paraMgr.wLop.setValue("Original Averaging KNN", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Increasing Step Size", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Local Neighbor Size For Inner Points", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Local Neighbor Size For Surface Points", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Inner Points Cooling Parameter", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Local Angle Threshold", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Eigen Neighborhood Para1", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Eigen Neighborhood Para2", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Average Closest Dist", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Average Dist To Input Threshold", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Original Confidence KNN", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Similarity Term Neighbor Para", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Density Confidence Segment Threshold", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Eigen Directional Threshold", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Save Move Dist Along Normal Para", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Big Repulsion Power", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Protect Small Tubular Para", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Protect High Confidence Para", DoubleValue(temp));

		infile >> temp; global_paraMgr.norSmooth.setValue("PCA KNN", IntValue(temp));
		infile >> temp; global_paraMgr.norSmooth.setValue("Sharpe Feature Bandwidth Sigma", DoubleValue(temp));
		infile >> temp; global_paraMgr.upsampling.setValue("Feature Sigma", DoubleValue(temp));
		infile >> temp; global_paraMgr.upsampling.setValue("Dist Threshold", DoubleValue(temp));
		infile >> temp; global_paraMgr.skeleton.setValue("Sigma KNN", DoubleValue(temp));
		infile >> temp; global_paraMgr.skeleton.setValue("Eigen Feature Identification Threshold", DoubleValue(temp));


		infile >> temp; global_paraMgr.wLop.setValue("Data Outweigh Similarity Para", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("Confidence Power", DoubleValue(temp));
		infile >> temp; global_paraMgr.wLop.setValue("KNN For Similarity", DoubleValue(temp));

	}


	if (!infile.eof())
	{
		infile >> temp; global_paraMgr.glarea.setValue("Sample Confidence Color Scale", DoubleValue(temp));
		infile >> temp; global_paraMgr.glarea.setValue("Point ISO Value Shift", DoubleValue(temp));

	}

	infile.close();
	emit needUpdateStatus();
}



void GLArea::wheelEvent(QWheelEvent *e) 
{
	const int WHEEL_STEP = 120;
	double change_rate = 0.1;
	double change = (e->delta() < 0) ? (1 + change_rate) : (1 - change_rate);

	double change_rate2 = 0.05;
	double change2 = (e->delta() < 0) ? (1 + change_rate2) : (1 - change_rate2);

	double size_temp = 0.0;

	
	if( (e->modifiers() & Qt::AltModifier) && (e->modifiers() & Qt::ControlModifier) )
	{
		if (para->getBool("Show Normal"))
		{
			size_temp = global_paraMgr.drawer.getDouble("Normal Line Length");
			global_paraMgr.drawer.setValue("Normal Line Length", DoubleValue(size_temp * change));
		}
		else if (para->getBool("Show Samples"))
		{
			size_temp = global_paraMgr.glarea.getDouble("Point ISO Value Shift");
			if (e->delta() < 0)
			{
				size_temp += 0.02;
			}
			else
			{
				size_temp -= 0.02;
			}
			global_paraMgr.glarea.setValue("Point ISO Value Shift", DoubleValue(size_temp));
			cout << "Point ISO Value Shift" << size_temp << endl;
		}

	}
	else if( (e->modifiers() & Qt::ShiftModifier) && (e->modifiers() & Qt::ControlModifier) )
	{
    if (para->getBool("Show Skeleton") /*&& !dataMgr.isSkeletonEmpty()*/)
    {
      size_temp = global_paraMgr.skeleton.getDouble("Branches Merge Max Dist");
      global_paraMgr.skeleton.setValue("Branches Merge Max Dist", DoubleValue(size_temp * change));
    }
		else if (global_paraMgr.drawer.getBool("Show Confidence Color"))
		{
			if (para->getBool("Show Samples"))
			{
				size_temp = global_paraMgr.glarea.getDouble("Sample Confidence Color Scale");
				global_paraMgr.glarea.setValue("Sample Confidence Color Scale", DoubleValue(size_temp * change));
				cout << "Sample Confidence Color Scale = " << size_temp * change << endl;
			}
		}
    else
    {
      size_temp = global_paraMgr.drawer.getDouble("Normal Line Width");
      global_paraMgr.drawer.setValue("Normal Line Width", DoubleValue(size_temp * change));
    }


	}
	else if((e->modifiers() & Qt::ShiftModifier) && (e->modifiers() & Qt::AltModifier))
	{
		if (para->getBool("Show Skeleton"))
		{
			size_temp = global_paraMgr.glarea.getDouble("Radius Ball Transparency") * change;
			global_paraMgr.glarea.setValue("Radius Ball Transparency", DoubleValue(size_temp));
			cout << "trans: " << size_temp << endl;
			if (size_temp < 0)
			{
				size_temp = 0;
			}
		}

		if ("Show Correspondences")
		{
			size_temp = global_paraMgr.glarea.getDouble("Show Confidence Percentage") * change2;
			global_paraMgr.glarea.setValue("Show Confidence Percentage", DoubleValue(size_temp));
			cout << "percentage: " << size_temp << endl;
			if (size_temp < 0)
			{
				size_temp = 0;
			}

		}

	}
	else
	{
		switch(e->modifiers())
		{
		case Qt::ControlModifier:

      if (para->getBool("Show Skeleton") && !dataMgr.isSkeletonEmpty())
      {
        size_temp = global_paraMgr.drawer.getDouble("Skeleton Node Size");
        global_paraMgr.drawer.setValue("Skeleton Node Size", DoubleValue(size_temp * change));
      }
			else if(para->getBool("Show Original") && para->getBool("Show Original Sphere") )
			{
				size_temp = global_paraMgr.drawer.getDouble("Original Draw Width");
				global_paraMgr.drawer.setValue("Original Draw Width", DoubleValue(size_temp * change));
			}
			else if(para->getBool("Show Samples") &&  para->getBool("Show Samples Dot")
				&&para->getBool("Show Original") && para->getBool("Show Original Dot") )
			{
				size_temp = global_paraMgr.drawer.getDouble("Sample Dot Size") * change;
				if(size_temp < 1)
				{
					size_temp = 1;
				}
				global_paraMgr.drawer.setValue("Sample Dot Size", DoubleValue(size_temp));

				size_temp = global_paraMgr.drawer.getDouble("Original Dot Size") * change;
				if(size_temp < 1)
				{
					size_temp = 1;
				}

				global_paraMgr.drawer.setValue("Original Dot Size", DoubleValue(size_temp));
			}
			else if (!para->getBool("Show Samples") && para->getBool("Show Dual Samples") &&  para->getBool("Show Dual Samples Dot"))
			{
				size_temp = global_paraMgr.drawer.getDouble("Dual Sample Dot Size") * change;
				if (size_temp < 1)
				{
					size_temp = 1;
				}
				global_paraMgr.drawer.setValue("Dual Sample Dot Size", DoubleValue(size_temp));
			}
			else if(para->getBool("Show Samples") &&  para->getBool("Show Samples Dot") )
			{
				size_temp = global_paraMgr.drawer.getDouble("Sample Dot Size") * change;
				if(size_temp < 1)
				{
					size_temp = 1;
				}
				global_paraMgr.drawer.setValue("Sample Dot Size", DoubleValue(size_temp));
			}
			else if (para->getBool("Show Samples") && (para->getBool("Show Samples Quad") || para->getBool("Show Samples Circle")))
			{
				size_temp = global_paraMgr.drawer.getDouble("Sample Draw Width") * change;
				if(size_temp < 0)
				{
					size_temp = 0.001;
				}
				cout << "draw width: " <<size_temp << endl;
				global_paraMgr.drawer.setValue("Sample Draw Width", DoubleValue(size_temp));
			}
			else if (para->getBool("Show Dual Samples") && !para->getBool("Show Dual Samples Dot"))
			{
				size_temp = global_paraMgr.drawer.getDouble("Dual Sample Draw Width") * change;
				if (size_temp < 0)
				{
					size_temp = 0.001;
				}
				//cout << "draw width: " <<size_temp << endl;
				global_paraMgr.drawer.setValue("Dual Sample Draw Width", DoubleValue(size_temp));
			}
			else if(para->getBool("Show Original") && para->getBool("Show Original Dot") )
			{
				size_temp = global_paraMgr.drawer.getDouble("Original Dot Size") * change;
				if(size_temp < 1)
				{
					size_temp = 1;
				}
				global_paraMgr.drawer.setValue("Original Dot Size", DoubleValue(size_temp));
			}
			else if (para->getBool("Show Original") && para->getBool("Show Original Circle") && !para->getBool("Show Dual Samples"))
			{
				size_temp = global_paraMgr.drawer.getDouble("Original Draw Width");
				global_paraMgr.drawer.setValue("Original Draw Width", DoubleValue(size_temp * change));
			}
			else
			{
				size_temp = global_paraMgr.drawer.getDouble("Sample Draw Width");
				global_paraMgr.drawer.setValue("Sample Draw Width", DoubleValue(size_temp * change));

				size_temp = global_paraMgr.drawer.getDouble("Original Draw Width");
				global_paraMgr.drawer.setValue("Original Draw Width", DoubleValue(size_temp * change));
			}
			emit needUpdateStatus();
			break;

		case Qt::ShiftModifier:

      if (para->getBool("Show Skeleton") && !dataMgr.isSkeletonEmpty())
      {
        size_temp = global_paraMgr.drawer.getDouble("Skeleton Bone Width");
        global_paraMgr.drawer.setValue("Skeleton Bone Width", DoubleValue(size_temp * change));
      }
      else
      {
        size_temp = global_paraMgr.data.getDouble("Down Sample Num");
        global_paraMgr.setGlobalParameter("Down Sample Num", DoubleValue(size_temp * change));
        emit needUpdateStatus();
      }

			break;

		case  Qt::AltModifier:
			size_temp = global_paraMgr.data.getDouble("CGrid Radius");
			if (size_temp < 1e-20)
			{
				size_temp = 0.05;
			}
			global_paraMgr.setGlobalParameter("CGrid Radius", DoubleValue(size_temp * change));
      global_paraMgr.setGlobalParameter("Initial Radius", DoubleValue(size_temp * change));

			//nitSetting();
			initView();
			wlop.setFirstIterate();
			skeletonization.setFirstIterate();
			emit needUpdateStatus();

			break;

		default:
			trackball.MouseWheel( e->delta()/ float(WHEEL_STEP));
			break;
		}
	}

	updateGL();
}

void GLArea::mouseMoveEvent(QMouseEvent *e)
{
	if (isRightPressed)
	{
		isDragging = true;
	}

	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);

	x2 = e->x();
	y2 = viewport[3] - e->y();

	if (isDefaultTrackBall())
	{
		trackball.MouseMove(e->x(),height()-e->y());
	}
	else 
	{
		trackball_light.MouseMove(e->x(),height()-e->y());
	}

	updateGL();
}

vcg::Trackball::Button QT2VCG(Qt::MouseButton qtbt,  Qt::KeyboardModifiers modifiers)
{
	int vcgbt= vcg::Trackball::BUTTON_NONE;
	if(qtbt & Qt::LeftButton		) vcgbt |= vcg::Trackball::BUTTON_LEFT;
	if(qtbt & Qt::RightButton		) vcgbt |= vcg::Trackball::BUTTON_RIGHT;
	if(qtbt & Qt::MidButton			) vcgbt |= vcg::Trackball::BUTTON_MIDDLE;
	if(modifiers & Qt::ShiftModifier		)	vcgbt |= vcg::Trackball::KEY_SHIFT;
	if(modifiers & Qt::ControlModifier ) vcgbt |= vcg::Trackball::KEY_CTRL;
	if(modifiers & Qt::AltModifier     ) vcgbt |= vcg::Trackball::KEY_ALT;
	return vcg::Trackball::Button(vcgbt);
}

void GLArea::mousePressEvent(QMouseEvent *e)
{


	if ((e->modifiers() & Qt::ShiftModifier) && (e->modifiers() & Qt::ControlModifier) &&
		(e->button()==Qt::LeftButton) )
		activeDefaultTrackball=false;
	else activeDefaultTrackball=true;

	if (isDefaultTrackBall())
	{
		if(e->button() == Qt::LeftButton)
			trackball.MouseDown(e->x(), height() - e->y(), QT2VCG(e->button(), e->modifiers() ) );     

		if(e->button() == Qt::RightButton) 
		{
			GLint viewport[4];
			glGetIntegerv (GL_VIEWPORT, viewport);
			x1 = e->x();
			y1 = viewport[3] - e->y();

			isRightPressed = true;
		}
	}
	else trackball_light.MouseDown(e->x(),height()-e->y(), QT2VCG(e->button(), Qt::NoModifier ) );

	//isDragging = true;
	update();
	updateGL();
}

void GLArea::mouseReleaseEvent(QMouseEvent *e)
{
	isDragging = false;

	if(e->button() == Qt::LeftButton)
		trackball.MouseUp(e->x(),height()-e->y(), QT2VCG(e->button(), e->modifiers() ) );

	if(e->button() == Qt::RightButton ) 
	{
		isRightPressed = false;
		GLint viewport[4];
		glGetIntegerv (GL_VIEWPORT, viewport);

		x2 = e->x();
		y2 = viewport[3] - e->y();

		double x3 = (x1 + x2)/2;
		double y3 = (y1 + y2)/2;

		if((e->modifiers() & Qt::ControlModifier))
		{
			pickPoint(x1,y1, fatherPickList, (x2-x1), (y1-y2), true);
			addPointByPick();
		}
		else if((e->modifiers() & Qt::AltModifier))
		{
			pickPoint(x1,y1, friendPickList, (x2-x1), (y1-y2), true);
			changePointByPick();
		}
		else if((e->modifiers() & Qt::ShiftModifier))
		{
			vector<int> no_use;
			addRBGPick( pickPoint(x1,y1, no_use, (x2-x1), (y1-y2)) );
		}
		else
		{
			pickPoint(x3,y3, pickList, abs(x2-x1), abs(y1-y2), !para->getBool("Multiply Pick Point"));
			friendPickList.clear();
			fatherPickList.clear();
			RGBPickList.assign(3, -1);
			RGB_counter = 0;

      if (pickList.empty())
      {
        para->setValue("Picked Index", DoubleValue(-1));
      }
      else
      {
        para->setValue("Picked Index", DoubleValue(pickList[0]));
      }


      //if (!pickList.empty())
      //{
      //  CMesh* dual_samples = dataMgr.getCurrentDualSamples();
      //  CVertex picked = dual_samples->vert[pickList[0]];
      //  picked_disk = NeighborDisk(picked.P(), picked.N());
      //  for (int i = 0; i < picked.neighbors.size(); i++)
      //  {
      //    picked_disk.add_point(dual_samples->vert[picked.neighbors[i]]);
      //  }
      //  picked_disk.printSlots();
      //  cout << picked_disk.getOccupyPercentage() << endl << endl;
      //}
		}
	}

	updateGL();
}



void GLArea::keyReleaseEvent ( QKeyEvent * e )
{
	if(e->key()==Qt::Key_Control) trackball.MouseUp(0,0, QT2VCG(Qt::NoButton, Qt::ControlModifier ) );
	if(e->key()==Qt::Key_Shift) trackball.MouseUp(0,0, QT2VCG(Qt::NoButton, Qt::ShiftModifier ) );
	if(e->key()==Qt::Key_Alt) trackball.MouseUp(0,0, QT2VCG(Qt::NoButton, Qt::AltModifier ) );
}


void GLArea::reorientPick()
{
//   CMesh* samples = dataMgr.getCurrentSamples();
// 
//   for (int i = 0; i < pickList.size(); i++)
//   {
//     samples->vert[pickList[i]].N() *= -1;
//   }

	if (global_paraMgr.wLop.getBool("Use Kite Points"))
	{
		CMesh* samples = dataMgr.getCurrentSamples();

		for (int i = 0; i < pickList.size(); i++)
		{
			samples->vert[pickList[i]].is_boundary = true;
		}

		return;
	}


  CMesh* samples = dataMgr.getCurrentSamples();

  for (int i = 0; i < pickList.size(); i++)
  {
    samples->vert[pickList[i]].is_fixed_sample = true;
		samples->vert[pickList[i]].is_skel_branch = true;
  }
}

void GLArea::cleanPick()
{
  CMesh* samples = dataMgr.getCurrentSamples();

	if (global_paraMgr.wLop.getBool("Use Kite Points"))
	{
		CMesh* samples = dataMgr.getCurrentSamples();

		for (int i = 0; i < pickList.size(); i++)
		{
			samples->vert[pickList[i]].is_boundary = false;
		}

		return;
	}


  for (int i = 0; i < pickList.size(); i++)
  {
    samples->vert[pickList[i]].is_fixed_sample = false;
		samples->vert[pickList[i]].is_skel_branch = false;
  }
}


void GLArea::removePickPoint()
{
	CMesh* samples = dataMgr.getCurrentSamples();
  CMesh* dual_samples = dataMgr.getCurrentDualSamples();
	bool use_random_erase = para->getBool("Random Erase");

	CMesh::VertexIterator vi;
	int j = 0;
	for(vi = samples->vert.begin(); vi != samples->vert.end(); ++vi, ++j)
	{
		vi->m_index = j;
		vi->neighbors.clear();
	}

	int step_size = 2;
	if (!use_random_erase)
	{
		step_size = 1;
	}
	if (!para->getBool("Multiply Pick Point"))
	{
		for (int i = 0; i < pickList.size(); i += step_size)
		{
			if(pickList[i] < 0 || pickList[i] >= samples->vert.size())
				continue;

			CVertex v = samples->vert[pickList[i]]; 
			int idx = v.m_index;
			samples->vert.erase(samples->vert.begin() + idx);
			dual_samples->vert.erase(dual_samples->vert.begin() + idx);
		}
	}
	else
	{		
		for (int i = 0; i < pickList.size(); i += step_size)
		{
			samples->vert[pickList[i]].is_skel_ignore = true;
      dual_samples->vert[pickList[i]].is_skel_ignore = true;

		}

		vector<CVertex> save_sample_vert;
		vector<CVertex> save_dual_sample_vert;

		for (int i = 0; i < samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];
			CVertex& dual_v = dual_samples->vert[i];

			if (!v.is_skel_ignore)
			{
				save_sample_vert.push_back(v);
				save_dual_sample_vert.push_back(dual_v);
			}
		}

		samples->vert.clear();
		dual_samples->vert.clear();

		for (int i = 0; i < save_sample_vert.size(); i++)
		{
			samples->vert.push_back(save_sample_vert[i]);
			dual_samples->vert.push_back(save_dual_sample_vert[i]);
		}
	}

	samples->vn = samples->vert.size();
	dual_samples->vn = dual_samples->vert.size();

	j = 0;
	for(vi = samples->vert.begin(); vi != samples->vert.end(); ++vi, ++j)
	{
		vi->m_index = j;
		dual_samples->vert[j].m_index = j;
	}

	cleanPickPoints();
}

void GLArea::sprayErasePick()
{
	cout << "spray erase" << endl;

	vector<CVertex> pick_union;
	CMesh* samples = dataMgr.getCurrentSamples();

	CMesh::VertexIterator vi;
	int j = 0;
	for (vi = samples->vert.begin(); vi != samples->vert.end(); ++vi, ++j)
	{
		vi->m_index = j;
		vi->neighbors.clear();
	}

	Point3f union_center(0., 0., 0.);
	int picked_size = 0;
	for (int i = 0; i < pickList.size(); i++)
	{
		if (pickList[i] < 0 || pickList[i] >= samples->vert.size())
			continue;

		CVertex &v = samples->vert[pickList[i]];
		
		pick_union.push_back(v);
		union_center += v.P();
		picked_size++;
	}
	union_center /= picked_size;
	cout << "picked size: " << picked_size;

	double max_dist2 = 0;
	for (int i = 0; i < picked_size; i++)
	{
		CVertex& v = pick_union[i];
		double dist2 = GlobalFun::computeEulerDistSquare(v.P(), union_center);
		if (dist2 > max_dist2)
		{
			max_dist2 = dist2;
		}
	}
	double threshold_dist2 = max_dist2 * 0.9 * 0.9;

	double iradius16 = -4 / max_dist2;
	for (int i = 0; i < picked_size; i++)
	{
		CVertex& v = pick_union[i];
		double dist2 = GlobalFun::computeEulerDistSquare(v.P(), union_center);
		if (dist2 > threshold_dist2)
		{
			continue;
		}
		double probability = exp(dist2 * iradius16);

		double r = (rand() % 1000) * 0.001;
		if (r < probability)
		{
			v.is_skel_ignore = true;
		}
	}

	for (int i = 0; i < picked_size; i++)
	{
		CVertex& v = pick_union[i];
		if (v.is_skel_ignore)
		{
			samples->vert[v.m_index].is_skel_ignore = true;
		}
	}


	vector<CVertex> save_sample_vert;
	for (int i = 0; i < samples->vert.size(); i++)
	{
		CVertex& v = samples->vert[i];
		if (!v.is_skel_ignore)
		{
			save_sample_vert.push_back(v);
		}
	}

	samples->vert.clear();
	for (int i = 0; i < save_sample_vert.size(); i++)
	{
		samples->vert.push_back(save_sample_vert[i]);
	}
	samples->vn = samples->vert.size();

	j = 0;
	for (vi = samples->vert.begin(); vi != samples->vert.end(); ++vi, ++j)
	{
		vi->m_index = j;
	}
	cleanPickPoints();

}


void GLArea::addPointByPick()
{
	if (dataMgr.isSamplesEmpty())
		return;

	if(pickList.empty() || fatherPickList.empty())
		return;

	CMesh* mesh = dataMgr.getCurrentSamples();
	CVertex newv;

	for(int ii = 0; ii < pickList.size(); ii++) 
	{
		int i = pickList[ii];

		if(i < 0 )
			return;

		CVertex &v = mesh->vert[i];   

		for(int jj = 0; jj < fatherPickList.size(); jj++)
		{
			int j = fatherPickList[jj];

			if(j < 0)
				return;

			CVertex &t = mesh->vert[j];

			if(v.P() == t.P())
			{
				cout << "overlap choose point!!" << endl;
				return;
			}

			//if(!v.neighbors.empty() && !t.neighbors.empty())
			//{
			//	dupsampler.setInput(&dataMgr);
			//	dupsampler.insertOnePointByHand(v.m_index, t.m_index);
			//	dupsampler.clear();
			//}
			//else
			//{
			newv = v;
			newv.P() = (v.P() + t.P()) / 2.0;  
			newv.m_index = mesh->vert.size();
			mesh->vert.push_back(newv);
			//}
		}
	}

	mesh->vn = mesh->vert.size();

	updateGL();
}

void GLArea::cleanPickPoints()
{
	pickList.clear();
	friendPickList.clear();
	fatherPickList.clear();
	RGBPickList.clear();
	glDrawer.cleanPickPoint();
	//glDrawer.weight_color_original.clear();
}

void GLArea::changePointByPick()
{
	if (dataMgr.isSamplesEmpty())
		return;

	CMesh* mesh = dataMgr.getCurrentSamples();
	CVertex newv;

	for(int ii = 0; ii < pickList.size(); ii++) 
	{
		int i = pickList[ii];

		if(i < 0 )
			continue;

		CVertex &v = mesh->vert[i];
		Point3f &p = v.P();     

		for(int jj = 0; jj < friendPickList.size(); jj++)
		{
			int j = friendPickList[jj];

			if(j < 0)
				continue;

			CVertex &t = mesh->vert[j];
			Point3f &q = t.P();
			v.N() = t.N();

			break;

		}
	}

	updateGL();
}


void GLArea::addRBGPick(int pick_index)
{
	vector<Point3f> RGB_normals;

	if(dataMgr.isSamplesEmpty())
		return;

	CMesh* mesh = dataMgr.getCurrentSamples();
	if(pick_index < 0 || pick_index >= mesh->vert.size())
	{
		cout << "pick GRB wrong!" << endl;
		return;
	}

	rotate_normal = mesh->vert[pick_index].N();
	rotate_pos = mesh->vert[pick_index].P();

	cout << "PICK!!" << endl;


// 	RGBPickList[RGB_counter++] = pick_index;
// 
// 	if(RGB_counter == 3)
// 	{
// 		cout << "constructRGBNormals" << endl;
// 
// 		vector<Point3f> normals;
// 		for(int i = 0; i < 3; i++)
// 		{
// 			normals.push_back(mesh->vert[RGBPickList[i]].N());
// 		}
// 
// 		RGB_normals.assign(3, Point3f(0.0, 0.0, 0.0));
// 		RGB_normals[0] = normals[0];
// 		RGB_normals[1] = normals[0] ^ normals[2];
// 
// 		if(RGB_normals[1] * normals[1] < 0)
// 			RGB_normals[1] = -RGB_normals[1];
// 
// 		RGB_normals[2] = normals[0]	^ normals[1];
// 		if(RGB_normals[2] * normals[2] < 0)
// 			RGB_normals[2] = -RGB_normals[2];
// 
// 
// 		glDrawer.setRGBNormals(RGB_normals);
// 		RGB_counter = 0;
// 	}

	emit needUpdateStatus();
}

void GLArea::readRGBNormal(QString fileName)
{
	ifstream infile(fileName.toStdString().c_str());

	vector<Point3f> rgb_normal(3, Point3f(0, 0, 0));
	for(int i = 0; i < rgb_normal.size(); i++)
	{
		infile >> rgb_normal[i].X() >> rgb_normal[i].Y() >> rgb_normal[i].Z();
	}

	glDrawer.setRGBNormals(rgb_normal);

	infile.clear();
}

void GLArea::drawCorrespondences()
{
	CMesh* samples = dataMgr.getCurrentSamples();
	CMesh* dual_samples = dataMgr.getCurrentDualSamples();
	CMesh* target_samples = dataMgr.getCurrentTargetSamples();
	CMesh* target_dual_samples = dataMgr.getCurrentTargetDualSamples();

	double normal_width = global_paraMgr.drawer.getDouble("Normal Line Width");
	double show_percentage = para->getDouble("Show Confidence Percentage");

	if (para->getBool("Show Dual Samples"))
	{
		for (int i = 0; i < samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];
			CVertex& dual_v = dual_samples->vert[i];

			if (v.eigen_confidence < show_percentage)
			{
				continue;
			}
			int target_index = v.target_index;
			if (target_index < 0)
			{
				continue;
			}
			CVertex& target_v = target_dual_samples->vert[target_index];
			glDrawer.glDrawLine(dual_v.P(), target_v.P(), cBlue, normal_width);
		}
	}
	else
	{
		for (int i = 0; i < samples->vert.size(); i++)
		{
			CVertex& v = samples->vert[i];
			if (v.eigen_confidence < show_percentage)
			{
				continue;
			}
			int target_index = v.target_index;
			if (target_index < 0)
			{
				continue;
			}
			CVertex& target_v = target_samples->vert[target_index];
			glDrawer.glDrawLine(v.P(), target_v.P(), cYellow, normal_width);
		}

	}
}


void GLArea::rotatingAnimation()
{
	if (-rotate_angle > 361)
		return;

	need_rotate = true;
	//rotate_angle -= 1.7f;
	rotate_angle -= 7.f;
	updateGL();


	cout << rotate_angle << endl;

	double total_time = 6500;
	double sleep_time = 45;
	double delta = 360 / (total_time / sleep_time);
	do
	{
		rotate_angle -= delta;
		updateGL();
		Sleep(sleep_time);
	} while (-rotate_angle < 361);

	rotate_angle += 360;
}