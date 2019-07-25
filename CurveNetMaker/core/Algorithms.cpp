/***************************************************************
 * Author:    Yixin Zhuang (yixin.zhuang@gmail.com)
 **************************************************************/


#include "core/Segmentation.h"

#include <fstream>
#include <queue>

//i use boost_1_67_0
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
using namespace boost;
//boost graph
typedef adjacency_list < vecS, vecS, undirectedS,
	no_property, property < edge_weight_t, double > > graph_t;
typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
typedef std::pair<int, int> Edge;

using namespace GeoProperty;

#include "andres/graph/graph.hxx"
#include "andres/graph/complete-graph.hxx"
#include "andres/graph/multicut-lifted/greedy-additive.hxx"
#include "andres/graph/multicut-lifted/kernighan-lin.hxx"

void Mitani_WaterShed::distanceToFeature_ShortestPath(std::list<MyMesh::Vertex>& vers, std::list<MyMesh::Edge>& edges, std::set<int>& sources, std::vector<double>& res)
{
	//build dual graph;

	int superNode = vers.size();
	int vNum = vers.size() + 1; //will add a super node
	int eNum = edges.size() + sources.size();//will add edges connecting super node and sources

	// init graph
	graph_t g(vNum);
	property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
	for (auto i = edges.begin(); i != edges.end(); i++)
	{		
		int nodeID[] = { i->vertex_iter(0)->id(), i->vertex_iter(1)->id() };
		graph_traits < graph_t >::edge_descriptor e;
		bool inserted;
		boost::tie(e, inserted) = add_edge(nodeID[0], nodeID[1], g);
		if (!inserted) cout << "insert edge failed" << endl;
		weightmap[e] = i->length();
	}
	for (auto i = sources.begin(); i != sources.end(); i++)
	{
		graph_traits < graph_t >::edge_descriptor e;
		bool inserted;
		boost::tie(e, inserted) = add_edge(*i, superNode, g);
		if (!inserted) cout << "insert edge failed" << endl;
		weightmap[e] = 0;
	}

	//shortest path;
	std::vector<vertex_descriptor> p(vNum);
	std::vector<double> d(vNum);
	vertex_descriptor s = vertex(superNode, g);
	dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));

	d.pop_back();
	double maxDist = *std::max_element(d.begin(), d.end());
	res.resize(d.size());
	for (int i = 0; i < d.size(); i++)
	{
		res[i] = maxDist - d[i];
	}
}
void Mitani_WaterShed::distanceToFeature_ShortestPath(MGraph& pg, std::set<int>& sources, std::vector<double>& res)
{
	//build dual graph;

	int superNode = pg.vNum;
	int vNum = pg.vNum + 1; //will add a super node
	int eNum = pg.eNum + sources.size();//will add edges connecting super node and sources

	// init graph
	graph_t g(vNum);
	property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
	for (auto i = pg.es.begin(); i != pg.es.end(); i++)
	{
		int nodeID[] = { i->n1, i->n2 };
		graph_traits < graph_t >::edge_descriptor e;
		bool inserted;
		boost::tie(e, inserted) = add_edge(nodeID[0], nodeID[1], g);
		if (!inserted) cout << "insert edge failed" << endl;
		weightmap[e] = i->w;
	}
	for (auto i = sources.begin(); i != sources.end(); i++)
	{
		graph_traits < graph_t >::edge_descriptor e;
		bool inserted;
		boost::tie(e, inserted) = add_edge(*i, superNode, g);
		if (!inserted) cout << "insert edge failed" << endl;
		weightmap[e] = 0;
	}

	//shortest path;
	std::vector<vertex_descriptor> p(vNum);
	std::vector<double> d(vNum);
	vertex_descriptor s = vertex(superNode, g);
	dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));

	d.pop_back();
	double maxDist = *std::max_element(d.begin(), d.end());
	res.resize(d.size());
	for (int i = 0; i < d.size(); i++)
	{
		res[i] = maxDist - d[i];
	}
}

int Mitani_WaterShed::grow(const std::vector<std::list<MyMesh::Vertex>::iterator>& vers, const std::vector<double>& distanceToFeature, std::vector<bool>& isExtrema, std::vector<unsigned>& vertexLabels)
{
	//init indexed items, and the flood fronts;
	int vNum = distanceToFeature.size();
	std::vector<bool> unIndexed(vNum, true);

	vertexLabels.clear(); vertexLabels.resize(vNum, UINT_MAX);
	int labelId = 0;
	std::priority_queue<HighFunction> floodFront;
	for (int i = 0; i < vNum; i++)
	{
		if (isExtrema[i])
		{
			unIndexed[i] = false;
			vertexLabels[i] = labelId++;
			floodFront.push(HighFunction(i, distanceToFeature[i]));
		}
	}

	const auto vts = vers;
	while (!floodFront.empty())
	{
		HighFunction v = floodFront.top(); floodFront.pop();

		for (int i = 0; i < vts[v.m_id]->vertex_iter().size(); i++)
		{
			int id = vts[v.m_id]->vertex_iter()[i]->id();
			if (unIndexed[id]) //unindexed neighbore;
			{
				HighFunction newV(id, distanceToFeature[id]);
				floodFront.push(newV);

				unIndexed[id] = false;
				vertexLabels[id] = vertexLabels[v.m_id];
			}
		}
	}

	cout << "cluster number:" << labelId << endl;
	return labelId;
}
//////////////////////////////////////////////////////////////////////////
//grow2 -- Implement of Paper: Analysis and Comparison of Algorithms for Morse Decompositions on Triangulated Terrains
//by Maria Vitali, Leila De Floriani, Paola Magillo
//detailed in section 5:watershed algorithms, page 9.
//the output contains not just the cluster labels of vertices, but also the boundary vertices(between basins). 
int Mitani_WaterShed::grow2(MGraph& pg, const std::vector<double>& distanceToFeature, const std::vector<std::set<unsigned> >& neighbours, std::vector<unsigned>& vertexLabels)
{
	//init indexed items, and the flood fronts;
	int vNum = distanceToFeature.size();
	std::vector<bool> unIndexed(vNum, true);
	std::priority_queue<HighFunction> floodFront;
	for (int i = 0; i < vNum; i++)
	{
		floodFront.push(HighFunction(i, distanceToFeature[i]));
	}
	
	int labelId = 0;
	std::vector<int> labels(vNum, -1);
	while (!floodFront.empty())
	{
		HighFunction v = floodFront.top(); floodFront.pop();
		std::set<int> ls;
		for (auto n = neighbours[v.m_id].begin(); n != neighbours[v.m_id].end(); n++)
		{
			if (!unIndexed[*n] && labels[*n] != -2) //neighbours assigned as basin;
			{
				ls.insert(labels[*n]);
			}
		}

		if (ls.empty())
		{
			labels[v.m_id] = labelId; labelId++;//new cluster begins...
		}
		else if (ls.size() == 1)
		{
			labels[v.m_id] = *ls.begin();//basin node
		}
		else
		{
			labels[v.m_id] = -2;//watershed node
		}

		unIndexed[v.m_id] = false;
	}

	vertexLabels.resize(vNum, 0);
	for (int i = 0; i < vNum; i++)
	{
		if (labels[i] == -2)
		{
			vertexLabels[i] = UINT_MAX; //watershed node
		}
	}

	return labelId;
}
void Mitani_WaterShed::thresholdClusters(MyMesh* mesh, std::vector<unsigned>& faceLabels, unsigned& segNum, std::vector<bool>& isFeature, double thres)
{
	// In our case, we use area as measurement of cluster size, for merging. Other criteria, like perimeter, would only require few changes of the code.
	double averageArea = 0;
	std::vector<double> faceAreas(mesh->getFaces().size());
	for (auto i = mesh->getFaces().begin(); i != mesh->getFaces().end(); i++)
	{
		faceAreas[i->id()] = i->triangleCost(2);
		averageArea += faceAreas[i->id()];
	}

	thres = averageArea * thres*0.01; //thres percent (thres%) of total area;

	//initialize area of all clusters;
	std::vector<double> scores(segNum, 0.0);
	for (auto i = mesh->getFaces().begin(); i != mesh->getFaces().end(); i++)
	{
		scores[faceLabels[i->id()]] += faceAreas[i->id()];
	}
	
	//shared boundary length between clusters;
	std::vector<double> sharesRow(segNum, 0.0);
	std::vector<std::vector<double> > sharesMatrix(segNum, sharesRow);
	for (auto i = mesh->getEdges().begin(); i != mesh->getEdges().end(); i++)
	{
		if (i->manifold())
		{
			int l1 = faceLabels[i->face_iter(0)->id()];
			int l2 = faceLabels[i->face_iter(1)->id()];

			if (l1 != l2 && !isFeature[i->face_iter(0)->id()] && !isFeature[i->face_iter(1)->id()])
			{
				sharesMatrix[l1][l2] = sharesMatrix[l2][l1] += i->length(); //in dual graph case, use i->dualLength;
			}
		}
	}

	//build adjacency and priority queue;
	std::multiset<ClusterArea> clusterQueue;
	for (int i = 0; i < scores.size(); i++)
	{
		clusterQueue.insert(ClusterArea(i, scores[i]));
	}
	if (clusterQueue.size() != scores.size())
	{
		cout << "numerical issue" << endl;
	}
	//adjcents between clusters;
	std::list<std::multiset<DecreaseOrderExt> > clusterAdjcencyList(segNum);
	std::vector<std::list<std::multiset<DecreaseOrderExt> >::iterator> clusterAdjcency;
	for (auto i = clusterAdjcencyList.begin(); i != clusterAdjcencyList.end(); i++)
	{
		clusterAdjcency.push_back(i);
	}

	for (auto i = mesh->getEdges().begin(); i != mesh->getEdges().end(); i++)
	{
		if (i->manifold())
		{
			int l1 = faceLabels[i->face_iter(0)->id()];
			int l2 = faceLabels[i->face_iter(1)->id()];

			if (l1 != l2)
			{
				if (clusterAdjcency[l1]->find(DecreaseOrderExt(l2, scores[l2], sharesMatrix[l1][l2])) == clusterAdjcency[l1]->end())
					clusterAdjcency[l1]->insert(DecreaseOrderExt(l2, scores[l2], sharesMatrix[l1][l2]));

				if (clusterAdjcency[l2]->find(DecreaseOrderExt(l1, scores[l1], sharesMatrix[l1][l2])) == clusterAdjcency[l2]->end())
					clusterAdjcency[l2]->insert(DecreaseOrderExt(l1, scores[l1], sharesMatrix[l1][l2]));
			}
		}
	}

	//clusters
	std::list<std::set<int>> clusters(segNum);
	std::vector<std::list<std::set<int>>::iterator> iclusters;
	for (auto i = clusters.begin(); i != clusters.end(); i++)
	{
		iclusters.push_back(i);
	}
	for (int i = 0; i < faceLabels.size(); i++)
	{
		iclusters[faceLabels[i]]->insert(i);
	}

// 	cout << "thres:" << thres << endl;
	while (!clusterQueue.empty())
	{
		ClusterArea ci = *clusterQueue.begin(); clusterQueue.erase(clusterQueue.begin());
		if (ci.m_val > thres)
		{
			break;
		}
		if (clusterAdjcency[ci.m_id]->empty())
		{
			continue;
		}

		DecreaseOrderExt cj_ = *clusterAdjcency[ci.m_id]->begin();
		ClusterArea cj(cj_.m_id, cj_.m_val);

		clusterAdjcency[ci.m_id]->erase(clusterAdjcency[ci.m_id]->begin()); //remove cj from ci's

		auto tp = clusterAdjcency[cj.m_id]->begin();
		for (; tp != clusterAdjcency[cj.m_id]->end(); tp++)
		{
			if (tp->m_id == ci.m_id )
				break;
		}
		if (tp == clusterAdjcency[cj.m_id]->end())
		{
			cout << "invalid merge" << endl;
		}
		else
		{
			clusterAdjcency[cj.m_id]->erase(tp);
		}

		auto i = std::find(clusterQueue.begin(), clusterQueue.end(), cj);
		if (i == clusterQueue.end())
		{
			cout << "invalid reference" << endl;
		}
		else
		{
			clusterQueue.erase(i);
		}

		//0. merge ci to cj.
		//1. cj is still in clusterQueue, need to update it, the new area is not just sum of ci and cj.
		//2. clear ci's clusterAdjacency, and add adj to cj, if not existed in cj. add cj to ci's adjacency, if not existed.
		
		//0. merge ci to cj.
		for (auto i = iclusters[ci.m_id]->begin(); i != iclusters[ci.m_id]->end(); i++)
		{
			faceLabels[*i] = cj.m_id;
		}
		iclusters[cj.m_id]->insert(iclusters[ci.m_id]->begin(), iclusters[ci.m_id]->end());
		iclusters[ci.m_id]->clear();


		//1. cj is still in clusterQueue, need to update it, the new area is  just sum of ci and cj.
		//2. clear ci's clusterAdjacency, and add it's adjacency to cj, if not existed in cj, and cj to it's adjacency, if not existed in...
		// ci's = {ci_a, ci_...}
		double newScore = ci.m_val + cj.m_val;
		clusterQueue.insert(ClusterArea(cj.m_id,newScore));

		for (auto i = clusterAdjcency[ci.m_id]->begin(); i != clusterAdjcency[ci.m_id]->end(); i++)
		{
			int ci_a = i->m_id; // for each cluster

			//if ci_a not in cj, then add to cj, if exist, then update their shared boundary, m_val_ext;
			//if cj not in ci_a, then add to ci_a, if exist, then update their shared boundary, m_val_ext;
			auto ip = clusterAdjcency[ci_a]->begin();
			for (; ip != clusterAdjcency[ci_a]->end(); ip++)
			{
				if (ip->m_id == cj.m_id && ip->m_val == cj.m_val)
					break;
			}
			if (ip == clusterAdjcency[ci_a]->end())
			{
				clusterAdjcency[ci_a]->insert(DecreaseOrderExt(cj.m_id, newScore, i->m_val_ext));
			}
			else
			{
				double temp = ip->m_val_ext;
				clusterAdjcency[ci_a]->erase(ip);
				clusterAdjcency[ci_a]->insert(DecreaseOrderExt(cj.m_id, newScore, i->m_val_ext + temp));
			}

			ip = clusterAdjcency[cj.m_id]->begin();
			for (; ip != clusterAdjcency[cj.m_id]->end(); ip++)
			{
				if (ip->m_id == i->m_id && ip->m_val == i->m_val)
					break;
			}
			if (ip == clusterAdjcency[cj.m_id]->end())
			{
				clusterAdjcency[cj.m_id]->insert(DecreaseOrderExt(i->m_id, i->m_val, i->m_val_ext));
			}
			else
			{
				double temp = ip->m_val_ext;
				clusterAdjcency[cj.m_id]->erase(ip);
				clusterAdjcency[cj.m_id]->insert(DecreaseOrderExt(i->m_id, i->m_val, i->m_val_ext + temp));
			}


			//remove ci from ci_a
			ip = clusterAdjcency[ci_a]->begin();
			for (; ip != clusterAdjcency[ci_a]->end(); ip++)
			{
				if (ip->m_id == ci.m_id && ip->m_val == ci.m_val)
					break;
			}
			if (ip == clusterAdjcency[ci_a]->end())
			{
				cout << endl << ci_a << " error " << ci.m_id << " " << ci.m_val << endl;
			}
			else
			{
				clusterAdjcency[ci_a]->erase(ip);
			}
		}

		clusterAdjcency[ci.m_id]->clear();

		//for all cj_a in cj, even though not appeared in ci, need to update cj_a's element cj...
		for (auto i = clusterAdjcency[cj.m_id]->begin(); i != clusterAdjcency[cj.m_id]->end();i++)
		{
			auto ip = clusterAdjcency[i->m_id]->begin();
			for (; ip != clusterAdjcency[i->m_id]->end(); ip++)
			{
				if (ip->m_id == cj.m_id)
					break;
			}
			if (ip == clusterAdjcency[i->m_id]->end())
			{
				cout << "bad reference" << endl;
			}
			else
			{
				clusterAdjcency[i->m_id]->erase(ip);
				clusterAdjcency[i->m_id]->insert(DecreaseOrderExt(cj.m_id, newScore,i->m_val_ext));
			}			
		}

	}

	//update label;
	int ind = 0;
	for (auto i = clusters.begin(); i != clusters.end(); i++)
	{
		if (i->empty())
			continue;

		for (auto f = i->begin(); f != i->end(); f++)
		{
			faceLabels[*f] = ind;
		}

		ind++;
	}
	segNum = ind;
}

void MeshSegment::convertFaceToVertexLabelling(std::vector<unsigned>& vertexLabel,std::vector<unsigned>& faceLabel)
{
	unsigned vNum = myMesh->getVertices().size();
	std::vector<std::vector<unsigned> > vertexMultiLabels(vNum);
	for (auto i = myMesh->getFaces().begin(); i != myMesh->getFaces().end(); i++)
	{
		for (unsigned j = 0; j < 3; j++)
		{
			vertexMultiLabels[i->vertex_iter(j)->id()].push_back(faceLabel[i->id()]);
		}
	}

	vertexLabel.clear(); vertexLabel.resize(vNum);
	for (int i = 0; i < vNum; i++)
	{
		std::map<unsigned, unsigned> st;
		std::map<unsigned, unsigned>::iterator it;
		for (int j = 0; j < vertexMultiLabels[i].size(); j++)
		{
			it = st.find(vertexMultiLabels[i][j]);
			if (it == st.end())
			{
				st[vertexMultiLabels[i][j]] = 1;
			}
			else
			{
				st[vertexMultiLabels[i][j]] = it->second + 1;
			}
		}

		std::pair<unsigned, unsigned> majorPair(vertexMultiLabels[i][0],1);
		for (it = st.begin(); it != st.end(); it++)
		{
			if (majorPair.second < it->second)
			{
				majorPair.first = it->first;
				majorPair.second = it->second;
			}
		}
		vertexLabel[i] = majorPair.first;
	}

	std::vector<double> weights;
	getWeights(edgeWeightParameter,weights);

	double cutCost = 0;
	for (auto e = myMesh->getEdges().begin(); e != myMesh->getEdges().end(); e++)
	{
		if (vertexLabel[e->vertex_iter(0)->id()] != vertexLabel[e->vertex_iter(1)->id()])
		{
			cutCost += weights[e->id()];
			graphFeature.ef[e->id()].isCut = true;
		}
		else
		{
			graphFeature.ef[e->id()].isCut = false;
		}
	}
	cout << "Cost:" << cutCost << endl;

	unsigned eNum = myMesh->getEdges().size();
	vertexLabel.clear(); vertexLabel.resize(vNum);
	segNumber = 0;
	std::vector<bool> visitedVer(vNum, false);
	std::vector<bool> visitedEdge(eNum, false);
	for (int i = 0; i < vNum; i++)
	{
		if (visitedVer[i] == true) continue;

		std::vector<int> frontVer(1, i);
		//propagates vertex with same label, and group them into a subgraph
		while (!frontVer.empty())
		{
			int fv = frontVer.back(); frontVer.pop_back();
			vertexLabel[fv] = segNumber;

			auto v = mg.vs[fv].vadjs.begin();
			auto e = mg.vs[fv].eadjs.begin();
			for (; v != mg.vs[fv].vadjs.end(); v++, e++)
			{
				if (graphFeature.ef[*e].isCut || visitedEdge[*e] == true)
					continue;

				visitedEdge[*e] = true;

				if (visitedVer[*v] == true)
					continue;

				visitedVer[*v] = true;

				frontVer.push_back(*v);
			}
		}
		segNumber++;
	}

	cout << endl << "cluster number :" << segNumber << endl;
}
void MeshSegment::convertVertexToFaceLabelling(std::vector<unsigned>& vertexLabel,std::vector<unsigned>& faceLabel)
{
	unsigned fNum = myMesh->getFaces().size();
	std::vector<std::vector<unsigned> > faceMultiLabels(fNum);
	for (auto i = myMesh->getFaces().begin(); i != myMesh->getFaces().end(); i++)
	{
		for (unsigned j = 0; j < 3; j++)
		{
			faceMultiLabels[i->id()].push_back(vertexLabel[i->vertex_iter(j)->id()]);
		}
	}
	
	faceLabel.clear();	faceLabel.resize(fNum);
	for (int i = 0; i < fNum; i++)
	{
		std::map<unsigned, unsigned> st;
		std::map<unsigned, unsigned>::iterator it;
		for (int j = 0; j < faceMultiLabels[i].size(); j++)
		{
			it = st.find(faceMultiLabels[i][j]);
			if (it == st.end())
			{
				st[faceMultiLabels[i][j]] = 1;
			}
			else
			{
				st[faceMultiLabels[i][j]] = it->second + 1;
			}
		}

		std::pair<unsigned, unsigned> majorPair(faceMultiLabels[i][0], 1);
		for (it = st.begin(); it != st.end(); it++)
		{
			if (majorPair.second < it->second)
			{
				majorPair.first = it->first;
				majorPair.second = it->second;
			}
		}

		faceLabel[i] = majorPair.first;
	}
}

///////////////////OVER SEGMENTATION///////////////////////////////
void MeshSegment::Mitani_Watershed()
{
	clock_t tstr = clock();

	//find sources;
	std::set<int> sources;
	for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
	{
		if (graphFeature.ef[e_it->id()].lab != -1)
		{
			sources.insert(e_it->vertex_iter(0)->id());
			sources.insert(e_it->vertex_iter(1)->id());
		}
	}

	ws.distanceToFeature_ShortestPath(myMesh->getVertices(), myMesh->getEdges(), sources, ws.distanceToFeature);

	//local extrema, and growing
	auto& vers = myMesh->getVIter();
	std::vector<bool> isExtrema(vers.size(), true);
	int ringSize = mitani_Watershed_Ringsize;
	for (int i = 0; i < vers.size(); i++)
	{
		int v1 = i;
		if (isExtrema[v1])
		{
			int times = ringSize;
			std::set<int> fronts; fronts.insert(v1);
			std::set<int> visited = fronts;
			std::set<int> nextFronts = fronts;
			while (times > 0)
			{
				for (auto j = fronts.begin(); j != fronts.end(); j++)
				{
					auto vj = vers[*j];
					for (int k = 0; k < vj->vertex_iter().size(); k++)
					{
						int v2 = vj->vertex_iter()[k]->id();
						if (visited.find(v2)!=visited.end())continue;
						nextFronts.insert(v2);
						visited.insert(v2);
					}
				}
				fronts = nextFronts;
				times--;
			}

			visited.erase(v1);
			for (auto j = visited.begin(); j != visited.end(); j++)
			{
				int v2 = *j;
				ws.distanceToFeature[v1] < ws.distanceToFeature[v2] ? isExtrema[v2] = false : isExtrema[v1] = false;
			}
		}
	}

	auto& multiCutLabels = vertexLabel;
	auto& multiCutNum = segNumber;
	multiCutNum = ws.grow(myMesh->getVIter(), ws.distanceToFeature, isExtrema, multiCutLabels);

	clock_t totalTime = clock() - tstr;
	cout << "Watershed Segmentation Use " << totalTime / 1000 << "sec" << totalTime % 1000 << "mm" << endl;

	std::vector<double> weights;
	getWeights(edgeWeightParameter, weights);

	vertexLabelBeforeMerge = vertexLabel;

	double cutCost = 0;
	for (auto e = myMesh->getEdges().begin(); e != myMesh->getEdges().end(); e++)
	{
		if (multiCutLabels[e->vertex_iter(0)->id()] != multiCutLabels[e->vertex_iter(1)->id()])
		{
			cutCost += weights[e->id()];
			graphFeature.ef[e->id()].isCut = true;
		}
		else
		{
			graphFeature.ef[e->id()].isCut = false;
		}
	}
	cout << "Cost:" << cutCost << endl;
}
void MeshSegment::Mitani_Watershed_Dual()
{
	clock_t tstr = clock();

	//build dual graph;
	std::list<MyMesh::Vertex> versArray(myMesh->getFaces().size());
	std::vector<std::list<MyMesh::Vertex>::iterator> vers(myMesh->getFaces().size());
	int ind = 0;
	for (auto v = versArray.begin(); v != versArray.end(); v++, ind++)
	{
		v->id() = ind;
		vers[ind] = v;
	}

	for (auto f = myMesh->getFaces().begin(); f != myMesh->getFaces().end(); f++)
	{
		vers[f->id()]->coordinate() = f->center_point();
	}

	std::vector<double> edgeCosts(myMesh->getEdges().size());
	std::vector<double> edgeLengths(myMesh->getEdges().size());
	{
		std::vector<Tensor>& tAnis = faceAnis;

		for (auto e_it = myMesh->getEdges().begin(); e_it != myMesh->getEdges().end(); e_it++)
		{
			if (e_it->manifold())
			{
				MyMesh::FaceIter f[] = { e_it->face_iter(0), e_it->face_iter(1) };
				unsigned fid[] = { f[0]->id(), f[1]->id() };
				Vec3 edgeVec = graphFeature.cf[f[0]->id()].rep - graphFeature.cf[f[1]->id()].rep;
				double edgeCost = 0;
				for (unsigned k = 0; k < 2; k++)
				{
					Vec2 edgeVec2D = Vec2(edgeVec.dot(tAnis[fid[k]].dir1), edgeVec.dot(tAnis[fid[k]].dir2));// two direction were swapped, since we computed dual edge's cost; and in faceAnis, mag1>mag2;
					edgeCost += sqrt(pow(edgeVec2D.x, 2)*tAnis[fid[k]].mag1 + pow(edgeVec2D.y, 2)*tAnis[fid[k]].mag2);
				}
				edgeCosts[e_it->id()] = edgeCost*0.5;
				edgeLengths[e_it->id()] = edgeVec.length();
			}
			else
			{
				MyMesh::FaceIter f[] = { e_it->face_iter(0) };
				unsigned fid[] = { f[0]->id(), f[1]->id() };
				Vec3 edgeVec = graphFeature.cf[f[0]->id()].rep - graphFeature.ef[e_it->id()].rep;
				double edgeCost = 0;
				for (unsigned k = 0; k < 1; k++)
				{
					Vec2 edgeVec2D = Vec2(edgeVec.dot(tAnis[fid[k]].dir1), edgeVec.dot(tAnis[fid[k]].dir2));// two direction were swapped, since we computed dual edge's cost; and in faceAnis, mag1>mag2;
					edgeCost += sqrt(pow(edgeVec2D.x, 2)*tAnis[fid[k]].mag1 + pow(edgeVec2D.y, 2)*tAnis[fid[k]].mag2);
				}
				edgeCosts[e_it->id()] = edgeCost + 1e-6; //make sure the boundary edge has weight larger than 0;
				edgeLengths[e_it->id()] = edgeVec.length();
			}
		}
	}

	std::list<MyMesh::Edge> edgesArray;
	for (auto e = myMesh->getEdges().begin(); e != myMesh->getEdges().end(); e++)
	{
		if (!e->manifold()) continue;

		int f1 = e->face_iter(0)->id();
		int f2 = e->face_iter(1)->id();

		vers[f1]->vertex_iter().push_back(vers[f2]);
		vers[f2]->vertex_iter().push_back(vers[f1]);

		MyMesh::Edge de;
		{//compute anis distance
			if (isAnisGeodesics_Watershed)
			{
				de.length() = edgeCosts[e->id()];
			}			
			else
			{
				de.length() = edgeLengths[e->id()];
			}
		}
		de.vertex_iter(0) = vers[f1];
		de.vertex_iter(1) = vers[f2];
		edgesArray.push_back(de);
	}
	std::vector<std::list<MyMesh::Edge>::iterator> edges(edgesArray.size());
	ind = 0;
	for (auto e = edgesArray.begin(); e != edgesArray.end(); e++, ind++)
	{
		e->id() = ind;
		edges[ind] = e;
	}

	//find sources;
	std::set<int> sources;
	for (unsigned i = 0; i<crestEdgesVisible.size(); i++)
	{
		if (!crestEdgesVisible[i])
			continue;

		unsigned faceID = crestEdges[i][2];
		if (faceID>myMesh->getFaces().size())
		{
			continue; // crest line data is not alway clean
		}
		sources.insert(faceID);
	}
	for (auto f = userSketches.begin(); f != userSketches.end(); f++)
	{
		sources.insert(f->fid);
	}

	ws.distanceToFeature_ShortestPath(versArray, edgesArray, sources, ws.distanceToFeature);

	//local extrema, and growing
	std::vector<std::set<unsigned> > neighbours(myMesh->getFaces().size());
	for (auto f = myMesh->getFaces().begin(); f != myMesh->getFaces().end(); f++)
	{
		auto& fs = neighbours[f->id()];
		for (unsigned j = 0; j < 3; j++)
		{
			for (auto e = 0; e < f->vertex_iter(j)->edge_iter().size(); e++)
			{
				fs.insert(f->vertex_iter(j)->edge_iter()[e]->face_iter(0)->id());
				if (f->vertex_iter(j)->edge_iter()[e]->manifold())
					fs.insert(f->vertex_iter(j)->edge_iter()[e]->face_iter(1)->id());
			}
		}
	}

// 	std::vector<bool> isExtrema(neighbours.size(), true);
// 	for (unsigned i = 0; i < neighbours.size(); i++)
// 	{
// 		for (auto j = neighbours[i].begin(); j != neighbours[i].end(); j++)
// 		{
// 			if(ws.distanceToFeature[i] > ws.distanceToFeature[*j])
// 			{
// 				isExtrema[i] = false;
// 			}
// 		}
// 	}

	std::vector<unsigned> watershedLabels;
	ws.grow2(mg,ws.distanceToFeature, neighbours, watershedLabels);

	clock_t totalTime = clock() - tstr;
	cout << "Watershed Segmentation Use " << totalTime / 1000 << "sec" << totalTime % 1000 << "mm" << endl;

	std::vector<double> weights;
	getWeights(edgeWeightParameter, weights);

	double cutCost = 0;
	for (auto e = myMesh->getEdges().begin(); e != myMesh->getEdges().end(); e++)
	{
		if (!e->manifold()) continue;

		if (watershedLabels[e->face_iter(0)->id()] == UINT_MAX && watershedLabels[e->face_iter(1)->id()] == UINT_MAX)
		{
			cutCost += weights[e->id()];
			graphFeature.ef[e->id()].isCut = true;
		}
		else
		{
			graphFeature.ef[e->id()].isCut = false;
		}
	}
	cout << "Cost:" << cutCost << endl;

	unsigned vNum = myMesh->getVertices().size();
	unsigned eNum = myMesh->getEdges().size();
	vertexLabel.clear(); vertexLabel.resize(vNum);
	segNumber = 0;
	std::vector<bool> visitedVer(vNum, false);
	std::vector<bool> visitedEdge(eNum, false);
	for (int i = 0; i < vNum; i++)
	{
		if (visitedVer[i] == true) continue;

		std::vector<int> frontVer(1, i);
		//propagates vertex with same label, and group them into a subgraph
		while (!frontVer.empty())
		{
			int fv = frontVer.back(); frontVer.pop_back();
			vertexLabel[fv] = segNumber;

			auto v = mg.vs[fv].vadjs.begin();
			auto e = mg.vs[fv].eadjs.begin();
			for (; v != mg.vs[fv].vadjs.end(); v++, e++)
			{
				if (graphFeature.ef[*e].isCut || visitedEdge[*e] == true)
					continue;

				visitedEdge[*e] = true;

				if (visitedVer[*v] == true)
					continue;

				visitedVer[*v] = true;

				frontVer.push_back(*v);
			}
		}
		segNumber++;
	}

	cout << endl << "cluster number:" << segNumber << endl;

	vertexLabelBeforeMerge = vertexLabel;
}
////////////////////MERGING//////////////////////////////////
void MeshSegment::LMP_Merging(const std::vector<unsigned>& vertexLabel, std::vector< std::list< std::list< MGTriple > >::iterator >&pCIter,
	unsigned labelId, std::vector<unsigned>& newLabels)
{
	unsigned numE = 0;
	unsigned numV = labelId;
	std::vector<std::pair<std::pair<unsigned, unsigned>, double> > edges;
	for (int i = 0; i < pCIter.size(); i++){
		for (auto tpc = pCIter[i]->begin(); tpc != pCIter[i]->end(); tpc++)
		{
			if (i != tpc->j)
			{
				numE++;
				edges.push_back(std::pair<std::pair<unsigned, unsigned>, double>(std::pair<unsigned, unsigned>(i, tpc->j), tpc->val));
			}
		}
	}
	cout << "number of edges:" << numE << "  number of vertices:" << numV << endl;

	andres::graph::Graph<> graph;
	graph.insertVertices(numV);
	std::vector<double> weights(numE);
	for (int i = 0; i < numE; i++)
	{
		graph.insertEdge(edges[i].first.first, edges[i].first.second);
		weights[i] = edges[i].second;
	}

	std::vector<char> edge_labels(graph.numberOfEdges(), 1);
	andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(graph, graph, weights, edge_labels);

	std::vector<char> out_labels(graph.numberOfEdges(), 1);
	andres::graph::multicut_lifted::kernighanLin(graph, graph, weights, edge_labels, out_labels);

	std::vector<unsigned> res(numV);
	std::vector<bool> visitedVer(numV, false);
	std::vector<bool> visitedEdge(numE, false);
	unsigned multiCutNum = 0;
	for (int i = 0; i < numV; i++)
	{
		if (visitedVer[i] == true) continue;

		std::vector<int> frontVer(1, i);
		//propagates vertex with same label, and group them into a subgraph
		while (!frontVer.empty())
		{
			int fv = frontVer.back(); frontVer.pop_back();
			res[fv] = multiCutNum;
			
			
			for (auto p = graph.adjacenciesToVertexBegin(fv); p != graph.adjacenciesToVertexEnd(fv); p++)
			{
				auto eid = p->edge();
				if (out_labels[eid] || visitedEdge[eid] == true)
					continue;

				visitedEdge[eid] = true;

				if (visitedVer[p->vertex()] == true)
					continue;

				visitedVer[p->vertex()] = true;

				frontVer.push_back(p->vertex());
			}
		}
		multiCutNum++;
	}
	//////////////////////////////////////////////////////////////////////////
	newLabels.clear(); newLabels.resize(vertexLabel.size());
	for (size_t i = 0; i < vertexLabel.size(); i++)
	{
		newLabels[i] = res[vertexLabel[i]];
	}
}
void MeshSegment::Mitani_Watershed_Dual_Partial(std::vector<unsigned>& fs)
{
	if (fs.empty())return;

	std::set<int> patchLabels;
	for (int i = 0; i < fs.size(); i++)
	{
		auto&f = myMesh->getFIter()[fs[i]];
		for (int j = 0; j < 3; j++)
		{
			patchLabels.insert(vertexLabel[f->vertex_iter(j)->id()]);
		}
	}

	std::vector<bool> faceInPatches(myMesh->getFaces().size(), true);
	subgraphVers.clear(); subgraphVers.resize(myMesh->getVertices().size(), true);
	unsigned numV = faceInPatches.size();
	for (auto f = myMesh->getFaces().begin(); f != myMesh->getFaces().end(); f++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (std::find(patchLabels.begin(), patchLabels.end(), vertexLabel[f->vertex_iter(j)->id()]) == patchLabels.end())
			{
				faceInPatches[f->id()] = false;
				subgraphVers[f->vertex_iter(j)->id()] = false;
			}
		}
		if (faceInPatches[f->id()] == false) numV--;
	}

	MGraph pg;//partial graph;
	pg.vNum = numV;
	pg.vs.clear(); pg.vs.resize(pg.vNum);

	std::map<unsigned, unsigned> localIndToGlobalInd;
	std::map<unsigned, unsigned> globalIndToLocalInd;
	unsigned ind = 0;
	for (auto f = myMesh->getFaces().begin(); f != myMesh->getFaces().end(); f++)
	{
		if (faceInPatches[f->id()])
		{
			pg.vs[ind].id = ind;
			pg.vs[ind].coord = f->center_point();
			localIndToGlobalInd[ind] = f->id();
			globalIndToLocalInd[f->id()] = ind;
			ind++;
		}
	}

	pg.es.clear();
	ind = 0;
	for (auto e = myMesh->getEdges().begin(); e != myMesh->getEdges().end(); e++)
	{
		if (!e->manifold()) continue;

		int f1 = e->face_iter(0)->id();
		int f2 = e->face_iter(1)->id();
		if (!faceInPatches[f1] || !faceInPatches[f2]) continue;

		f1 = globalIndToLocalInd[f1];
		f2 = globalIndToLocalInd[f2];

		pg.vs[f1].vadjs.push_back(f2);
		pg.vs[f2].vadjs.push_back(f1);

		MGraph::Edge de;
		de.w = e->cost();
		de.n1 = f1;
		de.n2 = f2;
		pg.es.push_back(de);

		pg.vs[f1].eadjs.push_back(ind);
		pg.vs[f2].eadjs.push_back(ind);
		ind++;
	}
	pg.eNum = pg.es.size();

	//find sources;
	std::set<int> sources;
	{
		for (unsigned i = 0; i < crestEdgesVisible.size(); i++)
		{
			if (!crestEdgesVisible[i])
				continue;

			unsigned faceID = crestEdges[i][2];
			if (faceID > myMesh->getFaces().size())
			{
				continue; // crest line data is not alway clean
			}
			if (!faceInPatches[faceID]) continue;

			sources.insert(globalIndToLocalInd[faceID]);
		}
		if (!graphCutLocally)
		{
			for (unsigned i = 0; i < userSketches.size(); i++)
			{
				if (!faceInPatches[userSketches[i].fid])continue;

				sources.insert(globalIndToLocalInd[userSketches[i].fid]);
			}
			for (unsigned i = 0; i < fs.size(); i++)
			{
				sources.insert(globalIndToLocalInd[fs[i]]);
			}
		}
	}

	std::vector<double> distanceToFeature;
	ws.distanceToFeature_ShortestPath(pg, sources, distanceToFeature);

	//localhood
	std::vector<std::set<unsigned> > neighbours(pg.vNum);
	for (auto f = myMesh->getFaces().begin(); f != myMesh->getFaces().end(); f++)
	{
		if (faceInPatches[f->id()] == false)continue;

		auto& fs = neighbours[globalIndToLocalInd[f->id()]];
		for (unsigned j = 0; j < 3; j++)
		{
			for (auto e = 0; e < f->vertex_iter(j)->edge_iter().size(); e++)
			{
				unsigned tid = f->vertex_iter(j)->edge_iter()[e]->face_iter(0)->id();
				if (faceInPatches[tid])
					fs.insert(globalIndToLocalInd[tid]);
				if (f->vertex_iter(j)->edge_iter()[e]->manifold())
				{
					tid = f->vertex_iter(j)->edge_iter()[e]->face_iter(1)->id();
					if (faceInPatches[tid])
						fs.insert(globalIndToLocalInd[tid]);
				}
			}
		}
	}

	std::vector<unsigned> watershedLabels;
	ws.grow2(pg, distanceToFeature, neighbours, watershedLabels);

	double cutCost = 0;
	for (auto e = myMesh->getEdges().begin(); e != myMesh->getEdges().end(); e++)
	{
		if (!e->manifold()) continue;

		if (vertexLabel[e->vertex_iter(0)->id()] != vertexLabel[e->vertex_iter(1)->id()])
			graphFeature.ef[e->id()].isCut = true;
		else
			graphFeature.ef[e->id()].isCut = false;

		int f1 = e->face_iter(0)->id();
		int f2 = e->face_iter(1)->id();

		if (!faceInPatches[f1] && !faceInPatches[f2])
			continue;
		else if (faceInPatches[f1] && !faceInPatches[f2] && watershedLabels[globalIndToLocalInd[f1]] == UINT_MAX)
		{
			graphFeature.ef[e->id()].isCut = true;
			continue;
		}
		else if (!faceInPatches[f1] && faceInPatches[f2] && watershedLabels[globalIndToLocalInd[f2]] == UINT_MAX)
		{
			graphFeature.ef[e->id()].isCut = true;
			continue;
		}

		if (watershedLabels[globalIndToLocalInd[f1]] == UINT_MAX && watershedLabels[globalIndToLocalInd[f2]] == UINT_MAX)
		{
			graphFeature.ef[e->id()].isCut = true;
		}
		else
		{
			graphFeature.ef[e->id()].isCut = false;
		}
	}

	unsigned vNum = myMesh->getVertices().size();
	unsigned eNum = myMesh->getEdges().size();
	vertexLabel.clear(); vertexLabel.resize(vNum);
	segNumber = 0;
	std::vector<bool> visitedVer(vNum, false);
	std::vector<bool> visitedEdge(eNum, false);
	for (int i = 0; i < vNum; i++)
	{
		if (visitedVer[i] == true) continue;

		std::vector<int> frontVer(1, i);
		//propagates vertex with same label, and group them into a subgraph
		while (!frontVer.empty())
		{
			int fv = frontVer.back(); frontVer.pop_back();
			vertexLabel[fv] = segNumber;

			auto v = mg.vs[fv].vadjs.begin();
			auto e = mg.vs[fv].eadjs.begin();
			for (; v != mg.vs[fv].vadjs.end(); v++, e++)
			{
				if (graphFeature.ef[*e].isCut || visitedEdge[*e] == true)
					continue;

				visitedEdge[*e] = true;

				if (visitedVer[*v] == true)
					continue;

				visitedVer[*v] = true;

				frontVer.push_back(*v);
			}
		}
		segNumber++;
	}

	cout << "cluster number:" << segNumber << endl;

	vertexLabelBeforeMerge = vertexLabel;

	//update initial cut graph;
	if (false)
	{
		std::vector<bool> isCut(myMesh->getEdges().size());
		for (auto e = myMesh->getEdges().begin(); e != myMesh->getEdges().end(); e++)
		{
			if (!e->manifold()) continue;

			if (vertexLabelInit[e->vertex_iter(0)->id()] != vertexLabelInit[e->vertex_iter(1)->id()])
				isCut[e->id()] = true;
			else
				isCut[e->id()] = false;

			int f1 = e->face_iter(0)->id();
			int f2 = e->face_iter(1)->id();

			if (!faceInPatches[f1] && !faceInPatches[f2])
				continue;
			else if (faceInPatches[f1] && !faceInPatches[f2] && watershedLabels[globalIndToLocalInd[f1]] == UINT_MAX)
			{
				isCut[e->id()] = true;
				continue;
			}
			else if (!faceInPatches[f1] && faceInPatches[f2] && watershedLabels[globalIndToLocalInd[f2]] == UINT_MAX)
			{
				isCut[e->id()] = true;
				continue;
			}

			if (watershedLabels[globalIndToLocalInd[f1]] == UINT_MAX && watershedLabels[globalIndToLocalInd[f2]] == UINT_MAX)
			{
				isCut[e->id()] = true;
			}
			else
			{
				isCut[e->id()] = false;
			}
		}
		vertexLabelInit.clear(); vertexLabelInit.resize(vNum);
		unsigned initSegNumber = 0;
		visitedVer.clear(); visitedVer.resize(vNum, false);
		visitedEdge.clear(); visitedEdge.resize(eNum, false);
		for (int i = 0; i < vNum; i++)
		{
			if (visitedVer[i] == true) continue;

			std::vector<int> frontVer(1, i);
			//propagates vertex with same label, and group them into a subgraph
			while (!frontVer.empty())
			{
				int fv = frontVer.back(); frontVer.pop_back();
				vertexLabelInit[fv] = initSegNumber;

				auto v = mg.vs[fv].vadjs.begin();
				auto e = mg.vs[fv].eadjs.begin();
				for (; v != mg.vs[fv].vadjs.end(); v++, e++)
				{
					if (isCut[*e] || visitedEdge[*e] == true)
						continue;

					visitedEdge[*e] = true;

					if (visitedVer[*v] == true)
						continue;

					visitedVer[*v] = true;

					frontVer.push_back(*v);
				}
			}
			initSegNumber++;
		}
	}

	//for debugging
	ws.distanceToFeature.clear();
	for (int i = 0; i < myMesh->getFaces().size(); i++)
	{
		if (globalIndToLocalInd.find(i) == globalIndToLocalInd.end())
			ws.distanceToFeature.push_back(0);
		else
			ws.distanceToFeature.push_back(distanceToFeature[globalIndToLocalInd[i]]);
	}
}
void MeshSegment::LMP_Partitioning()
{
	clock_t tstr = clock();

	unsigned vNum = myMesh->getVertices().size();
	unsigned eNum = myMesh->getEdges().size();

	andres::graph::Graph<> graph;
	graph.insertVertices(vNum);
	for (auto e = myMesh->getEdges().begin(); e != myMesh->getEdges().end(); e++)
	{
		graph.insertEdge(e->vertex_iter(0)->id(), e->vertex_iter(1)->id());
	}
	std::vector<double> weights;
	getWeights(edgeWeightParameter, weights);

	std::vector<char> edge_labels(graph.numberOfEdges(), 1);
	//andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(graph, graph, weights, edge_labels);

	std::vector<char> out_labels(graph.numberOfEdges(), 1);
	//andres::graph::multicut_lifted::kernighanLin(graph, graph, weights, edge_labels, out_labels);


	auto& multiCutLabels = vertexLabel;
	auto& multiCutNum = segNumber;
	multiCutLabels.clear(); multiCutLabels.resize(myMesh->getVertices().size());
	multiCutNum = 0;
	std::vector<bool> visitedVer(vNum, false);
	std::vector<bool> visitedEdge(eNum, false);
	for (int i = 0; i < vNum; i++)
	{
		if (visitedVer[i] == true) continue;

		std::vector<int> frontVer(1, i);
		//propagates vertex with same label, and group them into a subgraph
		while (!frontVer.empty())
		{
			int fv = frontVer.back(); frontVer.pop_back();
			multiCutLabels[fv] = multiCutNum;

			auto v = mg.vs[fv].vadjs.begin();
			auto e = mg.vs[fv].eadjs.begin();
			for (; v != mg.vs[fv].vadjs.end(); v++, e++)
			{
				if (out_labels[*e] || visitedEdge[*e] == true)
					continue;

				visitedEdge[*e] = true;

				if (visitedVer[*v] == true)
					continue;

				visitedVer[*v] = true;

				frontVer.push_back(*v);
			}
		}
		multiCutNum++;
	}

	clock_t tinit = clock() - tstr;
	cout << endl << "LMP takes:" << tinit / 1000 << "sec" << tinit % 1000 << "mm to generate " << multiCutNum << " patches." << endl;
	//////////////////////////////////////////////////////////////////////////

	vertexLabelBeforeMerge = vertexLabel;

	double cutCost = 0;
	for (int i = 0; i < eNum; i++)
	{
		if (out_labels[i])
		{
			cutCost += weights[i];
			graphFeature.ef[i].isCut = true;
		}
		else
		{
			graphFeature.ef[i].isCut = false;
		}
	}
	cout << "Cost:" << cutCost << endl;
}