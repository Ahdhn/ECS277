//Modified version of https://rosettacode.org/wiki/K-d_tree
//TODO: support insertion and deletion 

#ifndef _KDTREE_
#define _KDTREE_

#include <vector>
#include <algorithm>
#include <iostream>

#include "common.h"

class KdTree
{
public:
	KdTree(int dimension, int numPoints) :
		DIM(dimension), numPoints(numPoints){};

	//Call this one time to build everything 	
	void BuildTree(const double** Verts);
	void BuildTree(const std::vector<std::vector<double>>&Verts);
	
	template<typename T_d>
	void BuildTree(const scat_data_t<T_d> * m_data);

	//Call this to get the nearest node
	int FindNearest(double*point);
	int FindNearest(double xx, double yy, double zz);
	void FindNNearest(double *point, const  int N, std::vector<int>&nNearest);
	

	~KdTree(){};

private:

	void Swap(struct kd_node_t*, struct kd_node_t*);
	struct kd_node_t* FindMedian(struct kd_node_t*, struct kd_node_t*, int);
	struct kd_node_t* Construct(struct kd_node_t *t, int len, int i);

	double Dist(const struct kd_node_t *a, const struct kd_node_t *b);
	void getNearestRecursive(struct kd_node_t *root, const struct kd_node_t *testNode, int i, struct kd_node_t**best, double *best_dist);

	//bool inSkip(int id, int*skip);
	//void getNearestRecursive_skip(struct kd_node_t *root, const struct kd_node_t *testNode, int i, struct kd_node_t**best, double *best_dist, int*skip);

	void getNNearestRecursive(int N, struct kd_node_t *root, const struct kd_node_t *testNode, int i, struct kd_node_t **best, double *best_dist, int&nFound);



	const int DIM;
	const int numPoints;
	struct kd_node_t *root;



};
//*********** Construction 
void KdTree::BuildTree(const double** Verts)
{
	//Build the tree, copy data and update the root 
	/*if (dimension < 0 || numPoints <= 0){
	std::cerr << "Error in KdTree::BuildTree:: Invalid parameters" << std::endl;
	exit(1);
	}*/

	struct kd_node_t *myPoints = new kd_node_t[numPoints];
	for (int i = 0; i < numPoints; i++){
		myPoints[i].myID = i;
		myPoints[i].x = new double[DIM];
		for (int j = 0; j < DIM; j++){
			myPoints[i].x[j] = Verts[i][j];
		}
	}

	root = Construct(myPoints, numPoints, 0);

}
void KdTree::BuildTree(const std::vector<std::vector<double>>&Verts)
{
	//Build the tree, copy data and update the root 
	/*if (dimension < 0 || numPoints <= 0){
	std::cerr << "Error in KdTree::BuildTree:: Invalid parameters" << std::endl;
	exit(1);
	}*/

	struct kd_node_t *myPoints = new kd_node_t[numPoints];
	for (int i = 0; i < numPoints; i++){
		myPoints[i].myID = i;
		myPoints[i].x = new double[DIM];
		for (int j = 0; j < DIM; j++){
			myPoints[i].x[j] = Verts[i][j];
		}
	}

	root = Construct(myPoints, numPoints, 0);

}
template<typename T_d>
void KdTree::BuildTree(const scat_data_t<T_d> * m_data){

	struct kd_node_t *myPoints = new kd_node_t[numPoints];

	for (int i = 0; i < numPoints; i++){
		myPoints[i].myID = i;
		myPoints[i].x = new double[DIM];
		myPoints[i].x[0] = m_data[i].x;
		myPoints[i].x[1] = m_data[i].y;
		myPoints[i].x[2] = m_data[i].z;
	}

	root = Construct(myPoints, numPoints, 0);
}

struct kd_node_t* KdTree::Construct(struct kd_node_t *t, int len, int i)
{
	//Constrct the tree and update the root

	struct kd_node_t *n = NULL;
	if (n = FindMedian(t, t + len, i)){
		i = (i + 1) % DIM;
		n->left = Construct(t, n - t, i);
		n->right = Construct(n + 1, t + len - (n + 1), i);
	}
	return n;

}
void KdTree::Swap(struct kd_node_t*x, struct kd_node_t*y)
{
	for (int j = 0; j < DIM; j++){
		std::swap(x->x[j], y->x[j]);
	}
	std::swap(x->myID, y->myID);
}
struct kd_node_t* KdTree::FindMedian(struct kd_node_t*start, struct kd_node_t*end, int id)
{
	if (end <= start)return NULL;
	if (end == start + 1){ return start; }

	struct  kd_node_t*p, *store, *md(start + (end - start) / 2);
	double pivot;


	while (1){
		pivot = md->x[id];
		Swap(md, end - 1);



		for (store = p = start; p < end; p++){
			if (p->x[id] < pivot){
				if (p != store){
					Swap(p, store);
				}
				store++;
			}
		}

		Swap(store, end - 1);

		if (store->x[id] == md->x[id]){ return md; }

		if (store > md)end = store;
		else start = store;
	}

}

//*********** Query 
double KdTree::Dist(const struct kd_node_t *a, const struct kd_node_t *b)
{
	double t, d = 0;
	int myDim = DIM;
	while (myDim--) {
		t = a->x[myDim] - b->x[myDim];
		d += t * t;
	}
	return d;
}
int KdTree::FindNearest(double xx, double yy, double zz)
{
	//get the id of the nearest sample to 'point'
	struct kd_node_t testNode;
	testNode.x = new double[3];
	testNode.x[0] = xx;
	testNode.x[1] = yy;
	testNode.x[2] = zz;

	struct kd_node_t *best = 0;

	double best_dist;

	getNearestRecursive(root, &testNode, 0, &best, &best_dist);

	delete testNode.x;
	return best->myID;

}
int KdTree::FindNearest(double*point)
{
	//get the id of the nearest sample to 'point'
	struct kd_node_t testNode;
	testNode.x = point;

	struct kd_node_t *best = 0;

	double best_dist;

	getNearestRecursive(root, &testNode, 0, &best, &best_dist);

	return best->myID;

}
void KdTree::getNearestRecursive(struct kd_node_t *root, const struct kd_node_t *testNode, int i, struct kd_node_t**best, double *best_dist)
{
	if (!root){ return; }
	double d = Dist(root, testNode);
	double dx = root->x[i] - testNode->x[i];
	double dx2 = dx*dx;

	if (!*best || d < *best_dist){
		*best_dist = d;
		*best = root;
	}

	if (++i >= DIM){ i = 0; }

	getNearestRecursive(dx > 0 ? root->left : root->right, testNode, i, best, best_dist);
	if (dx2 >= *best_dist){ return; }
	getNearestRecursive(dx > 0 ? root->right : root->left, testNode, i, best, best_dist);

}

void KdTree::FindNNearest(double*point, const int N, std::vector<int>&nNearest)
{
	//get the id's of the N nearest samples to 'point'

	if (N >= numPoints){
		std::cout << "ERROR::KdTree::FindNNearest:: Can not find " << N << " nearest points in a " << numPoints << "-point-cloud" << std::endl;
		exit(1);

	}


	struct kd_node_t testNode;
	testNode.x = point;

	struct kd_node_t **best = new struct kd_node_t*[N];


	double best_dist = DBL_MAX;
	int nFound = 0;

	getNNearestRecursive(N, root, &testNode, 0, best, &best_dist, nFound);

	if (nFound != N){
		printf("\n Error::KdTree::FindNNearest:: can find all the N nearest points");
		exit(1);
	}

	//std::vector<int> nNearest;
	nNearest.clear();
	for (int i = 0; i < N; i++){
		nNearest.push_back(best[i]->myID);
	}

	//return nNearest;

}
void KdTree::getNNearestRecursive(int N, struct kd_node_t *root, const struct kd_node_t *testNode, int i, struct kd_node_t **best, double *best_dist, int&nFound)
{
	if (!root){ return; }
	double d = Dist(root, testNode); //sq distance 
	double dx = root->x[i] - testNode->x[i];
	double dx2 = dx*dx;

	if (d < *best_dist){

		root->dist = d;

		if (nFound < N){
			//if we have only too few
			best[nFound++] = root;
		}
		else{
			//we should update the furthest 
			int farID;
			double far_dist = -1.0;
			for (int i = 0; i < nFound; i++){
				if (far_dist < best[i]->dist){
					far_dist = best[i]->dist;
					farID = i;
				}
			}
			best[farID] = root;

			//only update the distance after you collect N number of samples
			//not effiecient but it works 
			//other wise, if you hit the closest sample is the root we won't go any further 
			*best_dist = best[0]->dist;
			for (int i = 1; i < nFound; i++){
				*best_dist = std::max(*best_dist, best[i]->dist);
			}
		}




	}

	if (++i >= DIM){ i = 0; }

	getNNearestRecursive(N, dx > 0 ? root->left : root->right, testNode, i, best, best_dist, nFound);
	if (dx2 >= *best_dist){ return; }
	getNNearestRecursive(N, dx > 0 ? root->right : root->left, testNode, i, best, best_dist, nFound);

}

#endif