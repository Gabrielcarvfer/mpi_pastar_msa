#ifndef WEIGHTEDSP_H
#define WEIGHTEDSP_H

#include "Sequences.h"
#include "Cost.h"
#include "HeuristicHPair.h"
#include <vector>

//Macros
#define MAX_SEQ_SIZE 1000
#define DASH '-'
#define BIG     999999
#define BIG_MIN_VAL 1.0E20
#define DIAG        0           /* code for traceback */
#define VERT        1           /* code for traceback */
#define HORZ        2           /* code for traceback */
#define GapCost     8

#define EfectiveGapCost 0//8


//Structs


class TreeNode
{
public:
    TreeNode(TreeNode * parent  =nullptr,
             TreeNode * leftSon =nullptr,
             TreeNode * rightSon=nullptr,
             TreeNode * brother =nullptr,
             float weight=1.0,
             float w=0.0,
             float W=0.0,
             float v=0.0,
             float V=0.0,
             int sequenceNumber=-1)
             {
             this->parent = parent;
             this->leftSon = leftSon;
             this->rightSon = rightSon;
             this->brother = brother;
             this->weight = weight;
             this->w = w;
             this->W = W;
             this->v = v;
             this->V = V;
             this->sequenceNumber = sequenceNumber;
             };


    //Pointers to build a tree
    TreeNode * parent;
    TreeNode * leftSon;
    TreeNode * rightSon;
    TreeNode * brother;

    float weight;
    float w;
    float W;
    float v;
    float V;
    int sequenceNumber; // 0 to N for a sequence
};


//Functions
void join_nodes(int min_i, int min_j, std::vector<TreeNode*> *tree,std::vector<TreeNode*> *nodesList);



float convert_path_to_cost(int I, int J,int n, int m, int **dd, int **hh, int **vv, std::vector<std::string>  Seqs);

void primer(std::vector<std::string>  seq, float ** Dij, float **scale);

/*! Recursive function to calculate the cost of path linking two tree leafs
*/
float compute_path_cost_rec(TreeNode *A, TreeNode *B, int * path_length, float *** Dij);

/*! Calculate path cost between two leafs and divide by path length linking them
*/
float compute_path_cost(int i, int j, std::vector<TreeNode*> *tree, float ***Dij);



/*! Build phylogenetic tree and calculate costs between leafs (sequences that are going to be aligned)
*/
void phylogeneticThreeNeighborJoin(std::vector<std::string>  seq, std::vector<TreeNode*> tree, std::vector<TreeNode*> nodesList, float *** Dij);

void compute_weights_from_tree(float product, float sum, TreeNode* no, TreeNode* brother, float *** weightMatrix, TreeNode ** pN);

/*! Calculate weights for weighted sum-of-pairs usings Altschul's rationale-2 method
*/
void weightAltschulsRationale2(Sequences * seqs);

int minOf3(int a, int b, int c);

#endif //WEIGHTEDSP_H
