/***************************************************************
 * Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
 **************************************************************/

#ifndef MYMESH
#define MYMESH
#include <vector>
#include <list>
#include <algorithm>
#include <math.h>
#include "amlVec.h"

typedef AML::double2 Vec2;
typedef AML::double3 Vec3;

class MyMesh
{
public:
	MyMesh(){};
	MyMesh(MyMesh*inMesh);
	~MyMesh(){};

	struct BBox{
		BBox(bool val = false) :isInit(val){}
		bool isInit;
		Vec3 bBoxMin;
		Vec3 bBoxMax;
		Vec3 centroid;
	};

	class Vertex;
	class Edge;
	class Face;

	typedef std::list<Vertex>::iterator VertexIter;
	typedef std::list<Edge>::iterator EdgeIter;
	typedef std::list<Face>::iterator FaceIter;

	class Vertex
	{
	public:
		Vertex(){
			m_id= 0;
			m_coords= Vec3(0,0,0);
			m_normal= Vec3(0,0,0);
			m_color= Vec3(0,0,0);
			m_dir[0]=Vec3(0,0,0);m_dir[1]= Vec3(0,0,0);
			m_mag[0]= 0;m_mag[1]= 0;
		};
		Vertex(const Vertex& in)
		{
			m_id= in.m_id;
			m_coords= in.m_coords;
			m_normal= in.m_normal;
			m_color= in.m_color;
			m_dir[0]=in.m_dir[0];m_dir[1]= in.m_dir[1];
			m_mag[0]= in.m_mag[0];m_mag[1]= in.m_mag[1];
			//m_face_iter=in.m_face_iter;
			m_edge_iter=in.m_edge_iter;
			m_vertex_iter=in.m_vertex_iter;
		};
		~Vertex(){};
		Vertex& operator=(const Vertex& in)
		{
			m_id= in.m_id;
			m_coords= in.m_coords;
			m_normal= in.m_normal;
			m_color= in.m_color;
			m_dir[0]=in.m_dir[0];m_dir[1]= in.m_dir[1];
			m_mag[0]= in.m_mag[0];m_mag[1]= in.m_mag[1];
			//m_face_iter=in.m_face_iter;
			m_edge_iter=in.m_edge_iter;
			m_vertex_iter=in.m_vertex_iter;
			return *this;
		}
		unsigned& id(){return m_id;}
		Vec3& coordinate(){return m_coords;}
		Vec3& color(){return m_color;}
		Vec3& normal(){return m_normal;}
		Vec3& direction(unsigned ind){return m_dir[ind];}
		double& magnitude(unsigned ind){return m_mag[ind];}

		std::vector<VertexIter>& vertex_iter(){return m_vertex_iter;}
		std::vector<EdgeIter>& edge_iter(){return m_edge_iter;}
		//std::vector<FaceIter>& face_iter(){return m_face_iter;}

	private:		
		unsigned m_id;
		Vec3 m_coords,m_normal,m_color;
		Vec3 m_dir[2];
		double m_mag[2];

		std::vector<VertexIter> m_vertex_iter;
		std::vector<EdgeIter> m_edge_iter;
		//std::vector<FaceIter> m_face_iter;
	};

	class Face
	{
	public:
		Face(){
			m_id=0;
			m_centre= Vec3(0,0,0);
			m_valid= true;
			m_cost=0;
		};
		~Face(){};
		Face(const Face& in)
		{
			m_id= in.m_id;
			m_centre= in.m_centre;
			m_vertex_iter[0] = in.m_vertex_iter[0];m_vertex_iter[1] = in.m_vertex_iter[1];m_vertex_iter[2] = in.m_vertex_iter[2];
			m_edge_iter[0]= in.m_edge_iter[0];m_edge_iter[1]= in.m_edge_iter[1];m_edge_iter[2]= in.m_edge_iter[2];
			m_valid= in.m_valid;
			m_cost=in.m_cost;
		};
		Face& operator=(const Face& in)
		{
			m_id= in.m_id;
			m_centre= in.m_centre;
			m_vertex_iter[0] = in.m_vertex_iter[0];m_vertex_iter[1] = in.m_vertex_iter[1];m_vertex_iter[2] = in.m_vertex_iter[2];
			m_edge_iter[0]= in.m_edge_iter[0];m_edge_iter[1]= in.m_edge_iter[1];m_edge_iter[2]= in.m_edge_iter[2];
			m_valid= in.m_valid;
			m_cost=in.m_cost;
			return *this;
		}

		unsigned& id(){return m_id;}
		Vec3& center_point(){return m_centre;}
		VertexIter& vertex_iter(unsigned i){return m_vertex_iter[i];}
		EdgeIter& edge_iter(unsigned i){return m_edge_iter[i];}
		bool& valid_triangle(){return m_valid;}
		double& cost(){return m_cost;}

		EdgeIter opposite_edge(VertexIter v)
		{
			for(unsigned i=0; i<3; ++i){
				EdgeIter e = m_edge_iter[i];
				if(e->vertex_iter(0)->id()!=v->id() && e->vertex_iter(1)->id()!=v->id()){
					return e;
				}
			}
		}
		VertexIter opposite_vertex(EdgeIter e)
		{
			for(unsigned i=0; i<3; ++i){
				VertexIter v = m_vertex_iter[i];				
				if(e->vertex_iter(0)->id()!=v->id() && e->vertex_iter(1)->id()!=v->id()){
					return v;
				}
			}
		}

		VertexIter next_vertex(VertexIter v)
		{
			for(unsigned i=0; i<3; ++i){
				if(v->id() == m_vertex_iter[i]->id()){
					return m_vertex_iter[(i+1)%3];
				}			
			}
		}

		bool triangleInequality()
		{
			double& l0 = m_edge_iter[0]->cost();
			double& l1 = m_edge_iter[1]->cost();
			double& l2 = m_edge_iter[2]->cost();

			if((l0+l1)>=l2 && (fabs(l0-l1)<=l2)){
				m_valid=true; return true;}
			else{
				/*triangleCost(1); */m_valid=false; return false;}
		}

		double triangleCost(unsigned k)
		{
			if(k==0){ //by perimeter of cost
				double l0 = m_edge_iter[0]->cost();
				double l1 = m_edge_iter[1]->cost();
				double l2 = m_edge_iter[2]->cost();
				m_cost= l0+l1+l2;			
			}
			if(k==1 || k==2){
				double a = m_edge_iter[0]->length();
				double b = m_edge_iter[1]->length();
				double c = m_edge_iter[2]->length();
				if(k==1) //by perimeter of length
					m_cost= a+b+c;
				if(k==2){//by area of length
					double p=(a+b+c)/2.;
					m_cost= /*p*2+*/sqrt(p*(p-a)*(p-b)*(p-c));
				}
			}
			return m_cost;
		}

	private:
		unsigned m_id;
		Vec3 m_centre;
		VertexIter m_vertex_iter[3];
		EdgeIter m_edge_iter[3];

		double m_cost;
		bool m_valid;
	};

	class Edge
	{
	public:
		Edge(){
			m_id= 0;
			m_length= 0;
			m_cost= 0;
			m_manifold= true;
		};
		~Edge(){};
		Edge(const Edge& in)
		{
			m_id= in.m_id;
			m_length= in.m_length;
			m_cost= in.m_cost;
			m_vertex_iter[0] = in.m_vertex_iter[0];m_vertex_iter[1] = in.m_vertex_iter[1];
			m_adjface_iter[0]= in.m_adjface_iter[0];m_adjface_iter[1]= in.m_adjface_iter[1];
			m_manifold= in.m_manifold;
		};
		Edge& operator=(const Edge& in)
		{
			m_id= in.m_id;
			m_length= in.m_length;
			m_cost= in.m_cost;
			m_vertex_iter[0] = in.m_vertex_iter[0];m_vertex_iter[1] = in.m_vertex_iter[1];
			m_adjface_iter[0]= in.m_adjface_iter[0];m_adjface_iter[1]= in.m_adjface_iter[1];
			m_manifold= in.m_manifold;
			return *this;
		}

		unsigned& id(){return m_id;}
		double& length(){return m_length;}
		double& cost(){return m_cost;}
		VertexIter& vertex_iter(unsigned i){return m_vertex_iter[i];}
		FaceIter& face_iter(unsigned i){return m_adjface_iter[i];}
		bool& manifold(){return m_manifold;}

		FaceIter opposite_face(FaceIter f)
		{
			if(!m_manifold)
				return f;
			if(m_adjface_iter[0]==f)
				return m_adjface_iter[1];
			else
				return m_adjface_iter[0];
		}

	private:
		unsigned m_id;
		double m_length;
		double m_cost;
		VertexIter m_vertex_iter[2];
		FaceIter m_adjface_iter[2];
		bool m_manifold;
	};

	std::list<Vertex>& getVertices(){return m_vertices;}
	std::list<Edge>& getEdges(){return m_edges;}
	std::list<Face>& getFaces(){return m_faces;}

	const std::vector<VertexIter>& getVIter(){return v_iters;}
	const std::vector<EdgeIter>& getEIter(){return e_iters;}
	const std::vector<FaceIter>& getFIter(){return f_iters;}

	void initialize(std::vector<Vec3> &inpoints, std::vector<std::vector<unsigned> > &infaces);

private:
	std::list<Vertex> m_vertices;
	std::list<Edge> m_edges;
	std::list<Face> m_faces;

	std::vector<VertexIter> v_iters;
	std::vector<EdgeIter> e_iters;
	std::vector<FaceIter> f_iters;

	void build_adjacencies(std::vector<VertexIter> &m_vertices_iter,std::vector<FaceIter> &m_faces_iter);
};//end of mesh
inline MyMesh::MyMesh(MyMesh* inMesh)
{
	if(inMesh==NULL) return;

	std::vector<VertexIter>& m_vertices_iter = v_iters;
	std::vector<EdgeIter>& m_edges_iter=e_iters;
	std::vector<FaceIter>& m_faces_iter=f_iters;
	m_vertices_iter.clear();
	m_edges_iter.clear();
	m_faces_iter.clear();

	unsigned v_num=inMesh->getVertices().size();
	m_vertices.resize(v_num);
	m_vertices_iter.resize(v_num);
	VertexIter v_it=m_vertices.begin();
	VertexIter iv_it=inMesh->getVertices().begin();
	for(unsigned i=0;i<v_num;i++){
		*v_it=*iv_it;
		m_vertices_iter[i]=v_it;
		v_it++;iv_it++;
	}
	unsigned e_num=inMesh->getEdges().size();
	m_edges.resize(e_num);
	m_edges_iter.resize(e_num);
	EdgeIter e_it=m_edges.begin();
	EdgeIter ie_it=inMesh->getEdges().begin();
	for(unsigned i=0;i<e_num;i++){
		*e_it=*ie_it;
		m_edges_iter[i]=e_it;
		e_it++;ie_it++;
	}
	unsigned f_num=inMesh->getFaces().size();
	m_faces.resize(f_num);
	m_faces_iter.resize(f_num);
	FaceIter f_it=m_faces.begin();
	FaceIter if_it=inMesh->getFaces().begin();
	for(unsigned i=0;i<f_num;i++){
		*f_it=*if_it;
		m_faces_iter[i]=f_it;
		f_it++;if_it++;
	}
	//update pointers;
	for(unsigned i=0;i<e_num;i++){
		m_edges_iter[i]->vertex_iter(0) = m_vertices_iter[m_edges_iter[i]->vertex_iter(0)->id()];
		m_edges_iter[i]->vertex_iter(1) = m_vertices_iter[m_edges_iter[i]->vertex_iter(1)->id()];
		m_edges_iter[i]->face_iter(0) = m_faces_iter[m_edges_iter[i]->face_iter(0)->id()];
		if(m_edges_iter[i]->manifold())
			m_edges_iter[i]->face_iter(1) = m_faces_iter[m_edges_iter[i]->face_iter(1)->id()];
	}
	for(unsigned i=0;i<f_num;i++){
		for(unsigned j=0;j<3;j++){
			m_faces_iter[i]->vertex_iter(j) = m_vertices_iter[m_faces_iter[i]->vertex_iter(j)->id()];
			m_faces_iter[i]->edge_iter(j) = m_edges_iter[m_faces_iter[i]->edge_iter(j)->id()];
		}
	}
	for(unsigned i=0;i<v_num;i++){
		for(unsigned j=0;j<m_vertices_iter[i]->vertex_iter().size();j++)
			m_vertices_iter[i]->vertex_iter()[j] = m_vertices_iter[m_vertices_iter[i]->vertex_iter()[j]->id()];
		for(unsigned j=0;j<m_vertices_iter[i]->edge_iter().size();j++)
			m_vertices_iter[i]->edge_iter()[j] = m_edges_iter[m_vertices_iter[i]->edge_iter()[j]->id()];
	}
}

inline void MyMesh::initialize(std::vector<Vec3> &inpoints, std::vector<std::vector<unsigned> > &infaces)
{
	std::vector<VertexIter>& m_vertices_iter = v_iters;
	std::vector<FaceIter>& m_faces_iter=f_iters;
	m_vertices_iter.clear();
	m_faces_iter.clear();

	m_vertices.resize(inpoints.size());
	unsigned i=0;
	for(MyMesh::VertexIter v_it=m_vertices.begin();v_it!=m_vertices.end();v_it++){
		m_vertices_iter.push_back(v_it);
		v_it->id() = i;
		v_it->coordinate()=inpoints[i];
		i++;
	}

	m_faces.resize(infaces.size());
	i=0;
	for(MyMesh::FaceIter f_it=m_faces.begin();f_it!=m_faces.end();f_it++){
		m_faces_iter.push_back(f_it);
		f_it->id() = i;
		f_it->center_point()=Vec3(0,0,0);
		for(unsigned j=0;j<3;j++){
			VertexIter &v_it = m_vertices_iter[infaces[i][j]];
			f_it->vertex_iter(j)= v_it;
			f_it->center_point()+=v_it->coordinate();
		}
		f_it->center_point()/=3.;
		i++;
	}
	build_adjacencies(m_vertices_iter,m_faces_iter);
}
inline void MyMesh::build_adjacencies(std::vector<VertexIter> &m_vertices_iter,std::vector<FaceIter> &m_faces_iter)
{
	//find all edges
	//i.e. find all half-edges, sort and combine them into edges
	std::vector<std::vector<unsigned> > half_edges(m_faces.size()*3);
	unsigned k = 0;
	unsigned i=0;
	for(MyMesh::FaceIter f_it=m_faces.begin();f_it!=m_faces.end();f_it++){
		for(unsigned j=0; j<3; ++j){
			half_edges[k].resize(3);
			half_edges[k][2] = i;
			unsigned vertex_id_1 = f_it->vertex_iter(j)->id();
			unsigned vertex_id_2 = f_it->vertex_iter((j+1)%3)->id();
			half_edges[k][0] = vertex_id_1<=vertex_id_2?vertex_id_1:vertex_id_2;
			half_edges[k][1] = vertex_id_1>vertex_id_2?vertex_id_1:vertex_id_2;

			k++;	
		}
		i++;
	}
	std::sort(half_edges.begin(), half_edges.end());

	unsigned number_of_edges = 1;
	for(unsigned i=1; i<half_edges.size(); ++i)
	{
		if(half_edges[i][0] != half_edges[i-1][0] || half_edges[i][1] != half_edges[i-1][1])
		{
			++number_of_edges;
		}
	}

	//		Edges->adjacent Vertices and Faces
	m_edges.resize(number_of_edges);
	EdgeIter e_it = m_edges.begin();
	unsigned edge_id = 0;
	for(unsigned i=0; i<half_edges.size();)
	{
		e_it->id() = edge_id++;
		e_it->vertex_iter(0) = m_vertices_iter[half_edges[i][0]];
		e_it->vertex_iter(1) = m_vertices_iter[half_edges[i][1]];
		e_it->length()=(e_it->vertex_iter(0)->coordinate()-e_it->vertex_iter(1)->coordinate()).length();
		if(i != half_edges.size()-1 && half_edges[i][0] == half_edges[i+1][0]&&half_edges[i][1] == half_edges[i+1][1])	//double edge
		{
			e_it->face_iter(0) = m_faces_iter[half_edges[i][2]];
			e_it->face_iter(1) = m_faces_iter[half_edges[i+1][2]];
			e_it->manifold()=true;
			i += 2;
		}
		else			//single edge
		{
			e_it->face_iter(0) = m_faces_iter[half_edges[i][2]];
			e_it->manifold()=false;
			i += 1;
		}
		e_it++;
	}

	//	Faces->adjacent Edges
	std::vector<unsigned> inds(m_faces.size(),0);
	for(MyMesh::EdgeIter e_it=m_edges.begin();e_it!=m_edges.end();e_it++){
		FaceIter f_it = e_it->face_iter(0);
		f_it->edge_iter(inds[f_it->id()])=e_it;
		inds[f_it->id()]++;
		if(e_it->manifold()){
			f_it = e_it->face_iter(1);
			f_it->edge_iter(inds[f_it->id()])=e_it;
			inds[f_it->id()]++;
		}
	}

	// Vertices->adjacent edges and vertices;
	for(MyMesh::EdgeIter e_it=m_edges.begin();e_it!=m_edges.end();e_it++){
		VertexIter v[] = {e_it->vertex_iter(0),e_it->vertex_iter(1)};
		v[0]->vertex_iter().push_back(v[1]);v[0]->edge_iter().push_back(e_it);
		v[1]->vertex_iter().push_back(v[0]);v[1]->edge_iter().push_back(e_it);
	}

	v_iters.clear();
	e_iters.clear();
	f_iters.clear();
	for(VertexIter v_it=m_vertices.begin();v_it!=m_vertices.end();v_it++)
		v_iters.push_back(v_it);
	for(EdgeIter e_it=m_edges.begin();e_it!=m_edges.end();e_it++)
		e_iters.push_back(e_it);
	for(FaceIter f_it=m_faces.begin();f_it!=m_faces.end();f_it++)
		f_iters.push_back(f_it);

}

#endif