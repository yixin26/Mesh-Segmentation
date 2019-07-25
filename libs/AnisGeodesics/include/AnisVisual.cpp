#include "AnisVisual.h"
#include "AnisGeodesic.h"
#include <algorithm>
#include <GL/glu.h>
#include <iostream>
using namespace std;
#define _USE_MATH_DEFINES 
#include <math.h>

unsigned AnisVisual::count;
unsigned AnisVisual::lineWidth=1;
unsigned AnisVisual::shapeSize=30;
unsigned AnisVisual::tensorShapeType=0;
bool AnisVisual::isDoSrcToAll=false;
bool AnisVisual::showDistanceField;
bool AnisVisual::showIsoline;
bool AnisVisual::showAllPaths;
bool AnisVisual::showInvalidTriangle=false;

void AnisVisual::drawBadFaces()
{
	if(anisGeodesy->anisMesh==NULL)return;

	glColor3f(.1,.1,.9);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_TRIANGLES);
	unsigned preID=0;
	auto f_it=anisGeodesy->anisMesh->getFaces().begin();
	for(unsigned i=0;i<anisGeodesy->badFaces.size();i++)
	{
		unsigned faceID=anisGeodesy->badFaces[i];
		std::advance(f_it,faceID-preID);
		preID=faceID;
		MyMesh::VertexIter v[]={f_it->vertex_iter(0),f_it->vertex_iter(1),f_it->vertex_iter(2)};
		std::vector<Vec3> triangle;
		for(unsigned j=0;j<3;j++)
		{
			triangle.push_back(v[j]->coordinate());
		}
		Vec3 norm = triangle[1]-triangle[0];
		norm = norm.cross(triangle[2]-triangle[1]);
		norm.normalize();

		glNormal3f(norm.x,norm.y,norm.z);
		for(unsigned j=0;j<3;j++)
		{
			glVertex3f(triangle[j].x,triangle[j].y,triangle[j].z);
		}
	}
	glEnd();
}

void AnisVisual::drawMesh()
{
	MyMesh* iMesh;
	if (anisGeodesy->anisMetric==EUCLIDEAN) iMesh=anisGeodesy->orgMesh;
	else iMesh=anisGeodesy->anisMesh;

	if(iMesh==NULL) return;
	glBegin(GL_TRIANGLES);
	for(auto f_it=iMesh->getFaces().begin();f_it!=iMesh->getFaces().end();f_it++)
	{
		MyMesh::VertexIter v[]={f_it->vertex_iter(0),f_it->vertex_iter(1),f_it->vertex_iter(2)};
		for(unsigned j=0;j<3;j++)
		{
			glNormal3f(v[j]->normal().x,v[j]->normal().y,v[j]->normal().z);
			glVertex3f(v[j]->coordinate().x,v[j]->coordinate().y,v[j]->coordinate().z);
		}
	}
	glEnd();
}
void AnisVisual::solidMeshList()
{
	glNewList(solidMesh,GL_COMPILE);
	drawSolidMesh();
	glEndList();
}
void AnisVisual::drawSolidMesh()
{
// 	glColor3d(.5, 0.539, 0.527);
	glColor3d(.5, 0.539, 0.527);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	drawMesh(); 	
	if(showInvalidTriangle) drawBadFaces();
}
void AnisVisual::wireMeshList()
{
	glNewList(wireMesh,GL_COMPILE);
	drawWireMesh();
	glEndList();
}
void AnisVisual::drawWireMesh(double w,Vec3 c)
{
	glColor3d(c.x,c.y,c.z);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glLineWidth(w);
	drawMesh();
}
void AnisVisual::tensorLineList(unsigned dir)
{
	if(dir==0)
		glNewList(tensorLineMax,GL_COMPILE);
	else
		glNewList(tensorLineMin,GL_COMPILE);	

	drawDirection(dir);

	glEndList();
}
void AnisVisual::drawDirection(unsigned dir)
{
	if(anisGeodesy->anisMesh==NULL)return;
	glLineWidth(1.0);
	if(dir==0)glColor3f(0.9,0.1,0.1);
	else glColor3f(0.1,0.1,0.9);

	unsigned maxNum = anisGeodesy->orgMesh->getVertices().size();
	MetricType mt = anisGeodesy->anisMetric;
	glBegin(GL_LINES);
	for(MyMesh::VertexIter v_it=anisGeodesy->anisMesh->getVertices().begin();v_it!=anisGeodesy->anisMesh->getVertices().end();v_it++)
	{
		if(maxNum<=v_it->id()) break;

		double mag[]={sqrt(v_it->magnitude(0)),sqrt(v_it->magnitude(1))};
		if (mt==EUCLIDEAN ||mt==POTTMANN || mt== CAMPEN_MAX || mt==CAMPEN_MIN)
		{
			if(mag[0]>mag[1]) std::swap(mag[0],mag[1]);
			mag[0] = mag[0]/mag[1];
			mag[1] = 1.;
		}
		else if (mt==KOVACS_MAX || mt== KOVACS_MIN)
		{
			if(anisGeodesy->maxAnis != 0 )
			{
				mag[0] /= anisGeodesy->maxAnis;
				mag[1] /= anisGeodesy->maxAnis;
			}
		}

		Vec3 &p=v_it->coordinate();
		Vec3 c1=v_it->direction(dir);
		c1 *= 0.0002*shapeSize*mag[dir];
		Vec3 c2=-c1;
		c1+=p;c2+=p;
		for(unsigned j=0;j<2;j++)
		{
			glVertex3f(c1.x,c1.y,c1.z);
			glVertex3f(c2.x,c2.y,c2.z);
		}
	}
	glEnd();
}
void AnisVisual::ellipse(const MyMesh::VertexIter& v_it,std::vector<double>& sin360,std::vector<double>& cos360)
{
	double mag1=sqrt(v_it->magnitude(0));
	double mag2=sqrt(v_it->magnitude(1));

	MetricType mt = anisGeodesy->anisMetric;
	if(mt==EUCLIDEAN || mt==POTTMANN)
	{
		mag1=1;
		mag2=1;
	}
	else if (mt== CAMPEN_MAX || mt==CAMPEN_MIN)
	{
		if(mag1>mag2) std::swap(mag1,mag2);
		mag1 = mag1/mag2;
		mag2 = 1.;
	}
	else if (mt==KOVACS_MAX || mt== KOVACS_MIN)
	{
		if(anisGeodesy->maxAnis != 0 )
		{
			mag1 /= anisGeodesy->maxAnis;
			mag2 /= anisGeodesy->maxAnis;
		}
	}

	Vec3 axis_x=v_it->direction(0);
	axis_x *= 0.0002*shapeSize*mag1;
	Vec3 axis_y=v_it->direction(1);
	axis_y *= 0.0002*shapeSize*mag2;

	Vec3 normal=v_it->normal();
	normal *= 0.0002;
	Vec3 origin=v_it->coordinate()+normal;

	Vec2 dir;
	std::vector<Vec3> ellipse(36);
	
	for (unsigned i=0; i<36; i++)
		ellipse[i]=origin+cos360[i]*axis_x+sin360[i]*axis_y;
	
	Vec3 n = v_it->normal();

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_POLYGON);
	for (int i=35; i>=0; i--) 
	{
		Vec3& p=ellipse[i];
		glNormal3d(n.x,n.y,n.z);
		glVertex3d(p.x,p.y,p.z);
	}
	glEnd();
}
void AnisVisual::tensorShape(const MyMesh::VertexIter& v_it,std::vector<double>& sin360,std::vector<double>& cos360)
{
	Vec3 axis_x=v_it->direction(0);
	axis_x *= 0.00015*shapeSize;
	Vec3 axis_y=v_it->direction(1);
	axis_y *= 0.00015*shapeSize;

	Vec3 normal=v_it->normal();
	normal *= 0.0002;
	Vec3 origin=v_it->coordinate()+normal;

	Vec2 dir;
	std::vector<Vec3> ellipse(36);
	
	double mag1=sqrt(v_it->magnitude(0));
	double mag2=sqrt(v_it->magnitude(1));
	if(mag1>mag2) std::swap(mag1,mag2);

	MetricType mt = anisGeodesy->anisMetric;
	if (mt==EUCLIDEAN ||mt==POTTMANN)
	{
		mag1 *= mag1;mag2 *= mag2;
	}

	for (unsigned i=0; i<36; i++) 
	{
		Vec3& p=ellipse[i];
		//Vec3& n=normals[i];

		double edgeCost= sqrt(pow(cos360[i],2)*mag1+pow(sin360[i],2)*mag2);
		p=origin+edgeCost*(cos360[i]*axis_x+sin360[i]*axis_y);
	}
	
	Vec3 n = v_it->normal();

	ellipse.push_back(ellipse.front());
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_POLYGON);
		for (unsigned i=0; i<19; i++) 
		{
			Vec3& p=ellipse[i];
			glNormal3d(n.x,n.y,n.z);
			glVertex3d(p.x,p.y,p.z);
		}
	glEnd();
	glBegin(GL_POLYGON);
		for (unsigned i=18; i<=36; i++) 
		{
			Vec3& p=ellipse[i];
			glNormal3d(n.x,n.y,n.z);
			glVertex3d(p.x,p.y,p.z);
		}
	glEnd();
}
void AnisVisual::tensorEllipseList()
{
	glNewList(tensorEllipse,GL_COMPILE);
	drawEllipse();
	glEndList();
}
void AnisVisual::drawEllipse()
{
	if(anisGeodesy->anisMesh==NULL)return;
	glColor3f(0,0,1);

	std::vector<double> sin360,cos360;
	cos360.resize(36);sin360.resize(36);
	for(unsigned i=0;i<36;i++)
	{
		cos360[i]=cos(double(i)*M_PI/18.);
		sin360[i]=sin(double(i)*M_PI/18.);
	}

	unsigned maxNum = anisGeodesy->orgMesh->getVertices().size();
	if(tensorShapeType==0)
	{
		for(MyMesh::VertexIter v_it=anisGeodesy->anisMesh->getVertices().begin();v_it!=anisGeodesy->anisMesh->getVertices().end();v_it++){
			if(maxNum<=v_it->id()) break;
			ellipse(v_it,sin360,cos360);
		}
	}
	else{
		for(MyMesh::VertexIter v_it=anisGeodesy->anisMesh->getVertices().begin();v_it!=anisGeodesy->anisMesh->getVertices().end();v_it++)
		{
			if(maxNum<=v_it->id()) break;
			tensorShape(v_it,sin360,cos360);
		}
	}

}
void AnisVisual::genDistanceFieldColor()
{
	if(!srcToAllPathsLength.empty())
	{
		double maxL = *std::max_element(srcToAllPathsLength.begin(),srcToAllPathsLength.end());
		double minL = *std::min_element(srcToAllPathsLength.begin(),srcToAllPathsLength.end());
		for(unsigned i=0;i<srcToAllPathsLength.size();i++)
		{
			if(srcToAllPathsLength[i]==-1.) srcToAllPathsLength[i]=maxL;
		}
		maxL-=minL; maxL= 0.66/maxL;

		std::vector<double> mags;
		for(unsigned i=0;i<srcToAllPathsLength.size();i++)
		{
			mags.push_back(0.66-(srcToAllPathsLength[i]-minL)*maxL);//reverse color
		}
		distanceFieldColor.clear();
		ColorEngine::HslToRgb(mags,distanceFieldColor);
	}
}
void AnisVisual::distanceFieldList()
{
	glNewList(distanceField,GL_COMPILE);
	drawDistanceField();
	glEndList();
}
bool AnisVisual::drawDistanceField()
{
	if(anisGeodesy->anisMesh==NULL) return false;

	if(!srcToAllPathsLength.empty())
	{
		glDisable(GL_LIGHTING);

		if(distanceFieldColor.empty()) genDistanceFieldColor();

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glBegin(GL_TRIANGLES);
		MyMesh* inMesh= anisGeodesy->orgMesh;
		for(auto f_it=inMesh->getFaces().begin();f_it!=inMesh->getFaces().end();f_it++)
		{
			MyMesh::VertexIter v[]={f_it->vertex_iter(0),f_it->vertex_iter(1),f_it->vertex_iter(2)};
			for(unsigned j=0;j<3;j++)
			{
				glColor3f(distanceFieldColor[v[j]->id()].r,distanceFieldColor[v[j]->id()].g,distanceFieldColor[v[j]->id()].b);
				glNormal3f(v[j]->normal().x,v[j]->normal().y,v[j]->normal().z);
				glVertex3f(v[j]->coordinate().x,v[j]->coordinate().y,v[j]->coordinate().z);
			}
		}
		glEnd();
		return true;
	}
	return false;
}
void AnisVisual::isoLineList()
{
	glNewList(isoline,GL_COMPILE);
	drawIsoline();	
	glEndList();
}
void AnisVisual::genIsoLine()
{
	if(!srcToAllPaths.empty())
	{
		isoLines.clear();
		unsigned num=30;
		isoLines.resize(num+1);

		double maxL = *std::max_element(srcToAllPathsLength.begin(),srcToAllPathsLength.end());
		double minL = *std::min_element(srcToAllPathsLength.begin(),srcToAllPathsLength.end());
		std::vector<unsigned> verTag(srcToAllPathsLength.size());
		std::vector<double> verVal(srcToAllPathsLength.size());
		for(unsigned i=0;i<srcToAllPathsLength.size();i++)
		{
			verVal[i]=double(num)*(srcToAllPathsLength[i]-minL)/maxL;
			verTag[i]=unsigned(verVal[i]);
		}

		MyMesh* inMesh= anisGeodesy->orgMesh;
		for(auto f_it=inMesh->getFaces().begin();f_it!=inMesh->getFaces().end();f_it++)
		{
			MyMesh::VertexIter v[]={f_it->vertex_iter(0),f_it->vertex_iter(1),f_it->vertex_iter(2)};
			for(unsigned i1=0;i1<3;i1++)
			{
				unsigned i2=(i1+1)%3;
				unsigned i3=(i1+2)%3;
				if(verTag[v[i1]->id()]!=verTag[v[i2]->id()] && verTag[v[i1]->id()]!=verTag[v[i3]->id()])
				{
					if(verTag[v[i2]->id()]==verTag[v[i3]->id()])
					{
						unsigned ind=min(verTag[v[i1]->id()],verTag[v[i2]->id()]);
						unsigned maxT=max(verTag[v[i1]->id()],verTag[v[i2]->id()]);
						for(unsigned j=0;j<maxT-ind;j++)
						{
							double jid=ind+j;
							double jmax=ind+j+1;
							double rat1= fabs(verVal[v[i1]->id()]-double(jmax))+1e-10;
							double rat2= fabs(verVal[v[i2]->id()]-double(jmax))+1e-10;
							double rsum = rat1+rat2;
							rat1 /= rsum; rat2/=rsum;
							Vec3 mid1= rat2*v[i1]->coordinate()+rat1*v[i2]->coordinate();

							rat1= fabs(verVal[v[i1]->id()]-double(jmax))+1e-10;
							rat2= fabs(verVal[v[i3]->id()]-double(jmax))+1e-10;
							rsum = rat1+rat2;
							rat1 /= rsum; rat2/=rsum;
							Vec3 mid2= rat2*v[i1]->coordinate()+rat1*v[i3]->coordinate();

							isoLines[jid].push_back(mid1);isoLines[jid].push_back(mid2);							
						}
					}
					else
					{
						std::vector<double> val_(3); val_[0]=verVal[v[i1]->id()];val_[1]=verVal[v[i2]->id()];val_[2]=verVal[v[i3]->id()];
						std::vector<unsigned> tag_(3); tag_[0]=verTag[v[i1]->id()];tag_[1]=verTag[v[i2]->id()];tag_[2]=verTag[v[i3]->id()];
						std::vector<unsigned> exactInd(3); exactInd[0]=i1;exactInd[1]=i2;exactInd[2]=i3;
						std::vector<unsigned> ind_(3);
						if(tag_[0]<tag_[1] && tag_[0]<tag_[2]) ind_[0]=0;
						else if(tag_[0]>tag_[1] && tag_[0]>tag_[2]) ind_[2]=0;
						else ind_[1]=0;
						if(tag_[1]<tag_[0] && tag_[1]<tag_[2]) ind_[0]=1;
						else if(tag_[1]>tag_[0] && tag_[1]>tag_[2]) ind_[2]=1;
						else ind_[1]=1;
						if(tag_[2]<tag_[0] && tag_[2]<tag_[1]) ind_[0]=2;						
						else if(tag_[2]>tag_[0] && tag_[2]>tag_[1]) ind_[2]=2;
						else ind_[1]=2;

						sort(val_.begin(),val_.end());sort(tag_.begin(),tag_.end());

						unsigned ind=tag_[0];
						unsigned maxT=tag_[1];
						for(unsigned j=0;j<maxT-ind;j++)
						{
							double jid=ind+j;
							double jmax=ind+j+1;

							double rat1= fabs(val_[0]-double(jmax))+1e-10;
							double rat2= fabs(val_[1]-double(jmax))+1e-10;
							double rsum = rat1+rat2;
							rat1 /= rsum; rat2/=rsum;
							Vec3 mid1= rat2*v[exactInd[ind_[0]]]->coordinate()+rat1*v[exactInd[ind_[1]]]->coordinate();

							rat1= fabs(val_[0]-double(jmax))+1e-10;
							rat2= fabs(val_[2]-double(jmax))+1e-10;
							rsum = rat1+rat2;
							rat1 /= rsum; rat2/=rsum;
							Vec3 mid2= rat2*v[exactInd[ind_[0]]]->coordinate()+rat1*v[exactInd[ind_[2]]]->coordinate();

							isoLines[jid].push_back(mid1);isoLines[jid].push_back(mid2);
						}


						ind=tag_[1];
						maxT=tag_[2];
						for(unsigned j=0;j<maxT-ind;j++)
						{
							double jid=ind+j;
							double jmax=ind+j+1;

							double rat1= fabs(val_[2]-double(jmax))+1e-10;
							double rat2= fabs(double(jmax) -val_[0])+1e-10;
							double rsum = rat1+rat2;
							rat1 /= rsum; rat2/=rsum;
							Vec3 mid1= rat2*v[exactInd[ind_[2]]]->coordinate()+rat1*v[exactInd[ind_[0]]]->coordinate();

							rat1= fabs(val_[2]-double(jmax))+1e-10;
							rat2= fabs(val_[1]-double(jmax))+1e-10;
							rsum = rat1+rat2;
							rat1 /= rsum; rat2/=rsum;
							Vec3 mid2= rat2*v[exactInd[ind_[2]]]->coordinate()+rat1*v[exactInd[ind_[1]]]->coordinate();

							isoLines[jid].push_back(mid1);isoLines[jid].push_back(mid2);
						}
					}
					break;
				}
			}
		}
	}
}

void AnisVisual::drawIsoline()
{
	if(anisGeodesy->anisMesh==NULL)return;
	if(!srcToAllPaths.empty())
	{
		if(isoLines.empty()) genIsoLine();
		unsigned num = isoLines.size()-1;
		std::vector<double> mags(num);
		double nor=0.66/double(num);
		for(unsigned i=0;i<num;i++)
		{
			mags[i]=double(i+1)*nor;
		}
		std::vector<ColorEngine::color> colorsList;
		ColorEngine::HslToRgb(mags,colorsList);
		reverse(colorsList.begin(),colorsList.end());

		glLineWidth(lineWidth);
		for(unsigned i=0;i<num;i++)
		{
			glColor3f(colorsList[i].r,colorsList[i].g,colorsList[i].b);
			glBegin(GL_LINES);
			for(unsigned j=0;j<isoLines[i].size();j+=2)
			{
				glVertex3f(isoLines[i][j].x,isoLines[i][j].y,isoLines[i][j].z);
				glVertex3f(isoLines[i][j+1].x,isoLines[i][j+1].y,isoLines[i][j+1].z);
			}
			glEnd();
		}
	}
}
void AnisVisual::allPathsList()
{
	glNewList(allPaths,GL_COMPILE);
	drawAllPaths();
	glEndList();	
}
std::vector<unsigned > ordering(const std::vector<double>&datas)
{
	std::vector<unsigned> ranks(datas.size(),0);
	std::vector<std::vector<unsigned> > rank; 
	rank.push_back(ranks);
	std::vector<std::vector<double>> data;
	data.push_back(datas);
	for(unsigned i=0;i<data.size();i++)
	{
		std::vector<double> tdat = data[i];
		for(unsigned j=0;j<tdat.size();j++)
		{
			tdat[j]/=2.;
		}
		for(unsigned j=0;j<tdat.size();j++)
		{
			unsigned minInd = j;
			for(unsigned k=0;k<tdat.size();k++)
			{
				if(tdat[minInd]>tdat[k])
					minInd = k;
			}
			rank[i][minInd]=j;
			tdat[minInd] = FLT_MAX;
		}
	}
	return rank.front();
}
void AnisVisual::genDisFieldRank()
{
	if(!srcToAllPathsLength.empty())
	{
		unsigned nodeSize=srcToAllPathsLength.size();
		std::vector<unsigned> randList = ordering(srcToAllPathsLength);
		disFieldMagnitudeRank.clear();
		disFieldMagnitudeRank.resize(nodeSize);
		for(unsigned i=0;i<nodeSize;i++)
			disFieldMagnitudeRank[randList[i]]=nodeSize-i-1;
	}
}
void AnisVisual::drawAllPaths()
{
	if(anisGeodesy->anisMesh==NULL)return;
	if(!srcToAllPaths.empty())
	{
		if(distanceFieldColor.empty()) genDistanceFieldColor();
		if(disFieldMagnitudeRank.empty()) genDisFieldRank();

		glLineWidth(lineWidth);
		for(unsigned i=0;i<srcToAllPaths.size();i++)
		{
			glColor3f(distanceFieldColor[disFieldMagnitudeRank[i]].r,
				distanceFieldColor[disFieldMagnitudeRank[i]].g,distanceFieldColor[disFieldMagnitudeRank[i]].b);
			glBegin(GL_LINE_STRIP);
			for(unsigned j=0;j<srcToAllPaths[disFieldMagnitudeRank[i]].size();j++)
			{
				glVertex3f(srcToAllPaths[disFieldMagnitudeRank[i]][j].x,srcToAllPaths[disFieldMagnitudeRank[i]][j].y,srcToAllPaths[disFieldMagnitudeRank[i]][j].z);
			}
			glEnd();
		}
	}
}
