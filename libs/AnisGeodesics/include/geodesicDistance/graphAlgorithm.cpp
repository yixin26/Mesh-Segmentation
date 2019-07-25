#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <float.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include "graphAlgorithm.h"

#include "fflas-ffpack/ffpack.h"
//#include "fflas-ffpack/modular-balanced.h"
#include "fflas-ffpack/modular-positive.h"
#include "Matio.h"
//typedef ModularBalanced<double> Field;
typedef Modular<double> Field;


graphAlgorithm::graphAlgorithm()
{
}
graphAlgorithm::~graphAlgorithm()
{
}

void graphAlgorithm::computeGeneralHamiltonGraph(const Capacity &org1,
												 const ArcWithWeight &org2,
												 const bool &org3,
												 MultipleArcs &tar)
{
	Capacity capacity=org1;
	ArcWithWeight arcsWithWeights=org2;
	bool connectRequire = org3;
	MultipleArcs usedArcsID;
	usedArcsID.clear();

    State initState;
	std::vector<State > queueState;

	initState.usedArcs.clear();
	initState.availableCapacity = capacity;
	initState.nextArcID = 0;
	initState.weight = 0;

    queueState.push_back(initState);

    while(true){
        double minWeight= FLT_MAX;
        int minID =-1;
        for(int i=0;i<queueState.size();i++){
            if(queueState[i].weight<minWeight){
                minID=i; minWeight=queueState[i].weight;
            }
        }
		if(minID==-1){
			usedArcsID.clear();
			break;}
        State currentState = queueState[minID];
        queueState.erase(queueState.begin()+minID);

		if(*(std::max_element(currentState.availableCapacity.begin(),
			currentState.availableCapacity.end()))==0){

			if(connectRequire){
				usedArcsID = currentState.usedArcs;
				std::vector<std::pair<int,int> > usedArcs;
				for(int i=0;i<usedArcsID.size();i++){
					int n1 = arcsWithWeights[usedArcsID[i].first].second.first;
					int n2 = arcsWithWeights[usedArcsID[i].first].second.second;
					usedArcs.push_back(std::pair<int,int>(n1,n2));
				}
				std::vector<int> usedNodes;
				usedNodes.push_back(usedArcs[0].first);
				usedNodes.push_back(usedArcs[0].second);
				usedArcs.erase(usedArcs.begin());
				while(!usedArcs.empty()){
					bool isChanged = false;
					int nodeSize=usedNodes.size();
					for(int j=0;j<nodeSize;j++){
						for(int k=0;k<usedArcs.size();k++){
							if(usedNodes[j]==usedArcs[k].first){
								if(find(usedNodes.begin(),usedNodes.end(),usedArcs[k].second)
									== usedNodes.end())
									usedNodes.push_back(usedArcs[k].second);
								usedArcs.erase(usedArcs.begin()+k);
								isChanged=true;
								break;
							}
							if(usedNodes[j]==usedArcs[k].second){
								if(find(usedNodes.begin(),usedNodes.end(),usedArcs[k].first)
									== usedNodes.end())
									usedNodes.push_back(usedArcs[k].first);
								usedArcs.erase(usedArcs.begin()+k);
								isChanged=true;
								break;
							}
						}
					}
					if(usedNodes.size()==capacity.size())
						break;
					if(!isChanged)break;
				}
				if(usedNodes.size()!=capacity.size())
					continue;
			}
			usedArcsID = currentState.usedArcs;
			break;
        }
        else{
		//expand current state;
			int nodeID[2];
			while(currentState.nextArcID<arcsWithWeights.size()){
				nodeID[0]=arcsWithWeights[currentState.nextArcID].second.first;
				nodeID[1]=arcsWithWeights[currentState.nextArcID].second.second;
				if(	currentState.availableCapacity[nodeID[0]]>=1
					&&currentState.availableCapacity[nodeID[1]]>=1)
					break;
				currentState.nextArcID++;
			}
			if(currentState.nextArcID<arcsWithWeights.size()){
				int degree = currentState.availableCapacity[nodeID[0]]<=
									currentState.availableCapacity[nodeID[1]]?
									currentState.availableCapacity[nodeID[0]]:
									currentState.availableCapacity[nodeID[1]];
				State saveState = currentState;
				int usedArcsID,isUsed=0;
				//currentState.weight += arcsWithWeights[currentState.nextArcID].first;
				currentState.nextArcID++;
				queueState.push_back(currentState);
				currentState=saveState;
				while(degree>0){
					currentState.availableCapacity[nodeID[0]]--;
					currentState.availableCapacity[nodeID[1]]--;
				//	int usedArcsID,isUsed=0;
					for(usedArcsID=0;usedArcsID<currentState.usedArcs.size();usedArcsID++){
						if(currentState.usedArcs[usedArcsID].first == (currentState.nextArcID)){
							isUsed=1;break;}
					}
					if(isUsed==1)
						currentState.usedArcs[usedArcsID].second++;
					else
						currentState.usedArcs.push_back(std::pair<int,int>(currentState.nextArcID,1));

					currentState.weight += arcsWithWeights[currentState.nextArcID].first;
					currentState.nextArcID++;
					queueState.push_back(currentState);
					currentState.nextArcID--;
					degree--;
				}
				currentState.nextArcID++;
			}
        }//end of else
    }//end of while

	tar = usedArcsID;
}

void graphAlgorithm::computeGeneralHamiltonGraph(const Capacity &org1,
								 const ArcWithWeight &org2, const PairArcsWithWeight & org3,
								 const double org4, const bool &org5, MultipleArcs &tar)
{
	Capacity capacity=org1;
	ArcWithWeight arcsWithWeights=org2;
	PairArcsWithWeight pairArcsWithWeights= org3;
	double balance = org4;
	bool connectRequire = org5;
	MultipleArcs usedArcsID;
	usedArcsID.clear();

	StateExtend initState;
	std::vector<StateExtend > queueState;

	initState.usedArcs.clear();
	initState.availableCapacity = capacity;
	initState.nextArcID = 0;
	initState.totalWeight = 0;
	initState.dihedralWeights.assign(capacity.size(),0);
	initState.adjacentArcs.resize(capacity.size());

	queueState.push_back(initState);

	while(true){
		double minWeight= FLT_MAX;
		int minID =-1;
		for(int i=0;i<queueState.size();i++){
			if(queueState[i].totalWeight<minWeight){
				minID=i; minWeight=queueState[i].totalWeight;
			}
		}
		if(minID==-1){
			usedArcsID.clear();	break;}

		StateExtend currentState = queueState[minID];
		queueState.erase(queueState.begin()+minID);

		if(*(std::max_element(currentState.availableCapacity.begin(),
			currentState.availableCapacity.end()))==0){

				if(connectRequire){
					usedArcsID = currentState.usedArcs;
					std::vector<std::pair<int,int> > usedArcs;
					for(int i=0;i<usedArcsID.size();i++){
						int n1 = arcsWithWeights[usedArcsID[i].first].second.first;
						int n2 = arcsWithWeights[usedArcsID[i].first].second.second;
						usedArcs.push_back(std::pair<int,int>(n1,n2));
					}
					std::vector<int> usedNodes;
					usedNodes.push_back(usedArcs[0].first);
					usedNodes.push_back(usedArcs[0].second);
					usedArcs.erase(usedArcs.begin());
					while(!usedArcs.empty()){
						bool isChanged = false;
						int nodeSize=usedNodes.size();
						for(int j=0;j<nodeSize;j++){
							for(int k=0;k<usedArcs.size();k++){
								if(usedNodes[j]==usedArcs[k].first){
									if(find(usedNodes.begin(),usedNodes.end(),usedArcs[k].second)
										== usedNodes.end())
										usedNodes.push_back(usedArcs[k].second);
									usedArcs.erase(usedArcs.begin()+k);
									isChanged=true;
									break;
								}
								if(usedNodes[j]==usedArcs[k].second){
									if(find(usedNodes.begin(),usedNodes.end(),usedArcs[k].first)
										== usedNodes.end())
										usedNodes.push_back(usedArcs[k].first);
									usedArcs.erase(usedArcs.begin()+k);
									isChanged=true;
									break;
								}
							}
						}
						if(usedNodes.size()==capacity.size())
							break;
						if(!isChanged)break;
					}
					if(usedNodes.size()!=capacity.size())
						continue;
				}
				usedArcsID = currentState.usedArcs;
				break;
		}
		else{
			//expand current state;
			int nodeID[2];
			while(currentState.nextArcID<arcsWithWeights.size()){
				nodeID[0]=arcsWithWeights[currentState.nextArcID].second.first;
				nodeID[1]=arcsWithWeights[currentState.nextArcID].second.second;
				if(	currentState.availableCapacity[nodeID[0]]>=1
					&&currentState.availableCapacity[nodeID[1]]>=1)
					break;
				currentState.nextArcID++;
			}
			if(currentState.nextArcID<arcsWithWeights.size()){
				int degree = currentState.availableCapacity[nodeID[0]]<=
							currentState.availableCapacity[nodeID[1]]?
							currentState.availableCapacity[nodeID[0]]:
							currentState.availableCapacity[nodeID[1]];
				StateExtend saveState = currentState;
				currentState.nextArcID++;
				queueState.push_back(currentState);
				currentState=saveState;
				while(degree>0){
					currentState.availableCapacity[nodeID[0]]--;
					currentState.availableCapacity[nodeID[1]]--;
					//	int usedArcsID,isUsed=0;
					int usedArcsID,isUsed=0;
					for(usedArcsID=0;usedArcsID<currentState.usedArcs.size();usedArcsID++){
						if(currentState.usedArcs[usedArcsID].first == (currentState.nextArcID)){
							isUsed=1;break;}
					}
					if(isUsed==1)
						currentState.usedArcs[usedArcsID].second++;
					else
						currentState.usedArcs.push_back(std::pair<int,int>(currentState.nextArcID,1));

					currentState.totalWeight += arcsWithWeights[currentState.nextArcID].first;

					//update weight;
					for(int n=0;n<2;n++){
						double w=0;
						for(int i=0;i<currentState.adjacentArcs[nodeID[n]].size();i++){
							if(w<pairArcsWithWeights[currentState.nextArcID][currentState.adjacentArcs[nodeID[n]][i]])
								w=pairArcsWithWeights[currentState.nextArcID][currentState.adjacentArcs[nodeID[n]][i]];
						}
						if(w>currentState.dihedralWeights[nodeID[n]]){
							currentState.totalWeight+=balance*(w-currentState.dihedralWeights[nodeID[n]]);
						}
						else{
							currentState.dihedralWeights[nodeID[n]]=w;
						}
						currentState.adjacentArcs[nodeID[n]].push_back(currentState.nextArcID);
					}

					currentState.nextArcID++;
					queueState.push_back(currentState);
					currentState.nextArcID--;
					degree--;
				}
				currentState.nextArcID++;
			}
		}//end of else
	}//end of while

	tar = usedArcsID;
}

void graphAlgorithm::computeCurveNormal(const std::vector<AML::double3> &org, AML::double3 &tar)
{
	std::vector<double> weights(org.size()-1);
	double size = double(org.size());
	for(int i=0;i<org.size()-1;i++)
		weights[i]=pow(1-double(i)/((size-1.0)==0?1:(size-1.0)),5);
	double sum=0.;
	for(int i=0;i<weights.size();i++)
		sum+=weights[i];
	for(int i=0;i<weights.size();i++)
		weights[i]/=sum;
	std::vector<AML::double3> vector(org.size()-1);
	for(int i=0;i<org.size()-1;i++){
		vector[i]=(org[i+1]-org[i]);
		vector[i].normalize();
	}
	for(int i=0;i<3;i++){
		sum=0;
		for(int j=0;j<weights.size();j++)
			sum+= weights[j]*vector[j][i];
		tar[i]=sum;
	}
	tar.normalize();
}

void graphAlgorithm::computeTransportMatrix(const std::vector<AML::double3> &org, std::vector<double> &tar)
{
	AML::double3 e1 = (org[1]-org[0]);
	e1.normalize();
	AML::double3 e2 = (org.back()-org[org.size()-2]);
	e2.normalize();;
	AML::double3 randomVector(1,0,0);
	AML::double3 coord[2];
	coord[0] = (e1.cross(randomVector));
	coord[0].normalize();
	coord[1] = e1.cross(coord[0]);
	AML::double3 coordTrans[]={coord[0],coord[1]};
	std::vector<AML::double3> vector(org.size()-1);
	for(int i=0;i<org.size()-1;i++){
		vector[i]=(org[i+1]-org[i]);
		vector[i].normalize();
	}
	AML::double3 tangent;
	for(int j=0;j<vector.size()-1;j++){
		tangent = (vector[j]+vector[j+1]);
		tangent.normalize();
		for(int i=0;i<2;i++){
			coordTrans[i] -= 2*(coordTrans[i].dot(tangent))*tangent;
		}
	}
	AML::double3 leftMatrix[3],rightMatrix[3];
	for(int i=0;i<3;i++){
		leftMatrix[i]= AML::double3(e1[i],coord[0][i],coord[1][i]);
		rightMatrix[i]=AML::double3(e2[i],coordTrans[0][i],coordTrans[1][i]);
	}
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			tar.push_back(leftMatrix[i].dot(rightMatrix[j]));
		}
	}
}

double graphAlgorithm::compute_angle(const AML::double3 &org1,const AML::double3 &org2)
{
	AML::double3 ang_str = org1;
	AML::double3 ang_end = org2;
	ang_str.normalize();
	ang_end.normalize();
	double res = ang_str.dot(ang_end);
	res = res < -1.0 ? -1.0 : res;
	res = res > 1.0 ? 1.0 : res;
	return acos(res);
}

//构造动态二维数组
double** create_Array(int n, int m)
{
	double** a = (double**) malloc(sizeof(double*)*n);
	for (int k=0; k<n; ++k)
		a[k] = (double*)malloc(sizeof(double)*m);
	return a;
}
//LU分解法回代过程，并求出方程组的解
std::vector<double> get_X(double** &L,double** &U,double* &b,int n)
{
	int i,j;

	double* y=(double*)malloc(sizeof(double)*n);
	double* x=(double*)malloc(sizeof(double)*n);

	//计算Ly=b中的y
	for(i=0;i<n;i++)    
	{
		y[i]=b[i];
		for(j=0;j<=i-1;j++)
			y[i]-=(L[i][j]*y[j]);
		y[i]/=L[i][i];
	}

	//计算Ux=y中的x
	for(i=n-1;i>=0;i--)   
	{
		x[i]=y[i];
		for(j=i+1;j<n;j++)
			x[i]-=(U[i][j]*x[j]);
		x[i]/=U[i][i];
	}

	std::vector<double> temp(n);
	for(i=0;i<n;i++){
		temp[i]=x[i];
	}
	free(y);
	free(x);
	return temp;
}

//释放动态二维数组的内存空间
void free_Array(double **a,int n)
{
	for(int m=0;m <n;m++)
		delete[]   a[m];
	delete[]   a;
}
void linearSystemSolver(const std::vector<std::vector<double> > &A,
						const std::vector<double> &ib, std::vector<double> &x)
{
	int n,m;
	int i,j,k;
	double sum;

	n=A.size();
	m=A.size();
	//建立并输入增广矩阵
	double** a=create_Array(n,m);
	for(i=0;i<n;i++){
		for(j=0;j<m;j++){
			a[i][j]=A[i][j];
		}
	}
	double *b=(double*)malloc(sizeof(double)*n);
/*	for(i=0;i<n;i++)
		b[i]=a[i][m-1];
*/
	for(i=0;i<n;i++)
		b[i]=ib[i];

	double** L=create_Array(n,n);
	double** U=create_Array(n,n);

	//初始化L矩阵
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
		{
			if(i==j)
				L[i][j]=1;
			else
				L[i][j]=0;
		}
	}
	//初始化U矩阵
	for(i=0; i<n; i++)
	{
		U[0][i] = (double)(a[0][i]/L[0][0]);
	}
	for(i=1;i<n;i++)
		for(j=0;j<n;j++)
			U[i][j]=0;

	//计算出L和U矩阵
	for(i=0; i<n-1; i++)
	{
		for(j=i+1; j<n; j++)
		{
			for(k=0,sum=0; k<n; k++)
			{
				if(k != i)
					sum += L[j][k]*U[k][i];
			}
			L[j][i] = (double)((a[j][i]-sum)/U[i][i]);
		}
		for(j=i+1; j<n; j++)
		{
			for(k=0,sum=0; k<n; k++)
			{
				if(k != i+1)
					sum += L[i+1][k]*U[k][j];
			}
			U[i+1][j] = (double)((a[i+1][j]-sum));
		}
	}
	
	//求解   
	x=get_X(L,U,b,n);
	
	//释放内存
	free_Array(a,n);
	free_Array(U,n);
	free_Array(L,n);
	free(b);
}

bool graphAlgorithm::tutteEmbedding(Graph &graph,std::vector<std::pair<int,AML::double3> > &fixedNode)
{
	std::vector<bool> isFixed(graph.nodes.size(),false);
	for(int i=0;i<fixedNode.size();i++){
		isFixed[fixedNode[i].first]=true;
	}
	std::vector<std::pair<int,AML::double3> > freeNode;
	std::vector<int> index(isFixed.size(),-1);
	for(int i=0;i<isFixed.size();i++){
		if(!isFixed[i]){
			AML::double3 temp;
			freeNode.push_back(std::pair<int,AML::double3>(i,temp));
			index[i]=freeNode.size()-1;
		}
	}	
	for(int i=0;i<fixedNode.size();i++){
		index[fixedNode[i].first]=i;
	}
	std::vector<std::vector<double> > A(freeNode.size(),std::vector<double>(freeNode.size(),0));
	for(int i=0;i<freeNode.size();i++){
		int nodeID = freeNode[i].first;
		int degree = graph.nodes[nodeID].arcID.size();
		for(int j=0;j<degree;j++){
			int arcID = graph.nodes[nodeID].arcID[j];
			int adjcentNodeID;
			if(graph.nodes[nodeID].arcDirection[j]==1){
				adjcentNodeID=graph.arcs[arcID].endNodesID.second;
			}
			else{
				adjcentNodeID=graph.arcs[arcID].endNodesID.first;
			}
			if(!isFixed[adjcentNodeID]){
				A[i][index[adjcentNodeID]]=1.0/double(degree);
			}
		}
		A[i][i]=-1.;
	}

	Field F(2.0);
	Field::Element zero;
	F.init(zero,0.0);
	Field::Element * X=NULL;
	int n=A.size();
	int p=n;
	X = new Field::Element[n*p];
	for (int i=0;i<n*p;++i)
		X[i] = zero;

	for (int i=0;i<n;++i){
		for (int j=0;j<p;++j){
			F.init(X[p*i+j],double(A[i][j]));
		}
	}

	Field::Element new_rank=0;
	new_rank = FFPACK::Det(F, n, p, X, p);
	if(new_rank==0) return false;

	AML::double3 temp(0,0,0);
	std::vector<AML::double3> vec_b(n,temp);
	for(int i=0;i<freeNode.size();i++){
		int nodeID = freeNode[i].first;
		int degree = graph.nodes[nodeID].arcID.size();
		for(int j=0;j<degree;j++){
			int arcID = graph.nodes[nodeID].arcID[j];
			int adjcentNodeID;
			if(graph.nodes[nodeID].arcDirection[j]==1){
				adjcentNodeID=graph.arcs[arcID].endNodesID.second;
			}
			else{
				adjcentNodeID=graph.arcs[arcID].endNodesID.first;
			}
			if(isFixed[adjcentNodeID]){
				vec_b[i].x -= fixedNode[index[adjcentNodeID]].second.x/double(degree);
				vec_b[i].y -= fixedNode[index[adjcentNodeID]].second.y/double(degree);
			}
		}
	}

	std::vector<double> vec_b_x(n);
	std::vector<double> vec_x_x(n);
	std::vector<double> vec_b_y(n);
	std::vector<double> vec_x_y(n);
	for (int i=0;i<n;++i){
		vec_b_x[i]=double(vec_b[i][0]);
		vec_b_y[i]=double(vec_b[i][1]);
	}

	linearSystemSolver(A,vec_b_x,vec_x_x);
	linearSystemSolver(A,vec_b_y,vec_x_y);

	for(int i=0;i<graph.nodes.size();i++){
		if(isFixed[i]){
			graph.nodes[i].pos=fixedNode[index[i]].second;
		}
		else{
			AML::double3 temp(vec_x_x[index[i]],vec_x_y[index[i]],0.);
			graph.nodes[i].pos=temp;
		}
	}
	for(int i=0;i<graph.arcs.size();i++){
		graph.arcs[i].vertexList.clear();
	}
	return true;
}
