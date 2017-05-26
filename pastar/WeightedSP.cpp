/*!
 * \authors {Gabriel Ferreira}
 * \copyright MIT License
 *  Mostly rewriting Altschul et al Weighted Sum-of-Pairs code and planning to implement Gotoh's alternative
 *  You can find the original source of MSA in https://www.ncbi.nlm.nih.gov/CBBresearch/Schaffer/msa.html
 */

#include "include/WeightedSP.hpp"

#include <iostream>
#include <string>

/*
void phylogeneticThreeUpgma(std::vector<std::string>  seq)
{

}
*/



inline int minOf3(int a, int b, int c)
{
    if (a < b)
        if (a < c)
            return a;
        else
            return c;
    else
        if (b < c)
            return b;
        else
            return c;

}

void rpt(TreeNode * A)
{
	if (A->sequenceNumber >= 0)
        printf("Leaf #%d        Distance to parent = %7.2f\n",1+A->sequenceNumber,A->weight);
	else
	{
		if (A->sequenceNumber == -2)
            printf("---------------- Tree given from ancestor ----------------\n");
		else
            printf("Internal Node  Distance to parent = %7.2f\n",A->weight);
		printf("On the left:   ");
		rpt(A->leftSon);
		printf("On the right:  ");
		rpt(A->rightSon);
	}
}


float compute_path_cost_to_leafs(TreeNode * A, float total, int * count2)
{
    if (A->sequenceNumber >= 0)
        return(total + A->weight);
    (*count2)++;
    return(compute_path_cost_to_leafs(A->leftSon,A->weight+total,count2) + compute_path_cost_to_leafs(A->rightSon,A->weight+total,count2));
}

float compute_curr_cost(int i, int j, std::vector<TreeNode*>*tree,float *** Dij)
{
	float diz=0.0,djz=0.0;
	int t, count2 = 1, nodesRemaining = tree->size();

	for (t=0; t<nodesRemaining; ++t) 
        if (t!=i && t!=j) 
        {
    		diz += compute_path_cost(i,t,tree,Dij);
    		djz += compute_path_cost(j,t,tree,Dij);
	    }
	diz = diz / (nodesRemaining - 2);
	djz = djz / (nodesRemaining - 2);
	return((compute_path_cost(i,j,tree,Dij) + diz - djz)/2 - compute_path_cost_to_leafs((*tree)[i],0.0,&count2)/count2);
}

void join_nodes(int min_i, int min_j, std::vector<TreeNode*> *tree,std::vector<TreeNode*> *nodesList,float ***Dij)
{
    TreeNode *leftNode,*rightNode,*newNode;
    leftNode = rightNode = newNode = nullptr;

    //Get pointer to left node and remove from node list
    leftNode = (*tree)[min_i];
    leftNode->weight = compute_curr_cost(min_i,min_j,tree, Dij);

    //Get pointer to right node and remove from node list
    rightNode = (*tree)[min_j];
    rightNode->weight = compute_curr_cost(min_j,min_i,tree,Dij);

    //Create new internal node
    newNode = new TreeNode(nullptr,leftNode,rightNode,nullptr,0.0,0.0, 0.0, 0.0, 0.0, -1);

    //Update previously removed nodes
    leftNode->brother = rightNode;
    rightNode->brother = leftNode;
    leftNode->parent = rightNode->parent = newNode;

    //Emplace new node to nodeslist
    (*nodesList).push_back(newNode);

    //Remove references of removed nodes and replace with new node plus latest node in tree list
    (*tree)[min_i] = newNode;
    (*tree)[min_j] = (*tree)[(*tree).size()-1];
    (*tree).erase((*tree).end()-1);
}

float convert_path_to_cost(int I, int J,int n, int m, int **dd, int **hh, int **vv, std::vector<std::string> *Seqs)
{
    int i,j,V,H,M;
    int dir=DIAG;   /* Direction of previous edge in traceback  */
    int match=0;    /* Number of matches in an optimal alignment    */
    float f,g;

    for (i=n,j=m; (i||j);)
    {
        V=vv[i][j]-(dir==VERT ? (j==m ? EfectiveGapCost:GapCost) : 0);
        H=hh[i][j]-(dir==HORZ ? (i==n ? EfectiveGapCost:GapCost) : 0);
        M=minOf3(V,H,dd[i][j]);

        if  (!j || M==V)
        {
            dir=VERT;
            --i;
        }
        else if (!i || M==H)
        {
            dir=HORZ;
             --j;
        }
        else
        {
            dir=DIAG;
            match+= (*Seqs)[I][i]==(*Seqs)[J][j];
            --i;
            --j;
        }
    }
    int ret = ((int) (0.5+1000.0*(n-match+m-match)/(n+m)));
    return ret;
}

void primer(std::vector<std::string> *seq, float *** Dij, float ***scale)
{
    int    I,J,i,j,Gi,Gj,n,m;
    int num_seq = seq->size();, max_seq_length = 1000;

    if (num_seq == 0)
        exit(-1);

    std::string seqA, seqB;

    //Allocate memory for diagonal, horizontal and vertical matrices
    int **dd = nullptr, **hh = nullptr,**vv= nullptr;

    dd = new int*[max_seq_length]();
    hh = new int*[max_seq_length]();
    vv = new int*[max_seq_length]();

    (*scale)[0][0] = 1;
    for (int k = 0; k < max_seq_length; k++)
    {
        dd[k] = new int[max_seq_length]();
        hh[k] = new int[max_seq_length]();
        vv[k] = new int[max_seq_length]();

        for (int l = 1; l < max_seq_length;l++)
        {
            dd[k][l] = dd[k][l-1]+max_seq_length;
            hh[k][l] = hh[k][l-1]+max_seq_length;
            vv[k][l] = vv[k][l-1]+max_seq_length;
            if (k < num_seq && l < num_seq)
                (*scale)[k][l] = 1;
        }
    }

    //For each sequence pair (I,J)
    for (I=0;I<num_seq-1;I++)
    {
        for (seqA=(*seq)[I],n=seqA.length(),J=I+1; J<num_seq; J++)
        {
            seqB=(*seq)[J];
            m=seqB.length();

            /* compute distance from <0,0> to <i,j> */
            dd[0][0]=0;
            hh[0][0]=vv[0][0]=EfectiveGapCost;

            //For each element in the sequence J
            for (j=1; j<=m; j++)
            {
                vv[0][j] = dd[0][j]   = BIG;
                hh[0][j] = hh[0][j-1] + Cost::cost(DASH,seqB[j]);
            }

            //For each element in sequence I
            for (i=1; i<=n; i++)
            {
                hh[i][0] = dd[i][0]   = BIG;
                vv[i][0] = vv[i-1][0] + Cost::cost(seqA[i],DASH);
            }

            //For each element in sequences I and J
            for (i=1; i<n; i++)
            {
                Gi = ( i==(n-1) ? EfectiveGapCost : GapCost);
                for (j=1; j<m; j++)
                {

                    Gj = (j==(m-1) ? EfectiveGapCost : GapCost);

                    dd[i][j] = minOf3(dd[i-1][j-1] , hh[i-1][j-1] , vv[i-1][j-1] ) + Cost::cost(seqA[i], seqB[j]);

                    hh[i][j] = minOf3(dd[i][j-1]+Gi, hh[i][j-1]   , vv[i][j-1]+Gi) + Cost::cost(DASH   , seqB[j]);

                    vv[i][j] = minOf3(dd[i-1][j]+Gj, hh[i-1][j]+Gj, vv[i-1][j]   ) + Cost::cost(seqA[i], DASH   );
                }
            }

            //Calculate path cost based on the path
            (*scale)[J][I] = convert_path_to_cost(I,J,n-1,m-1,dd,hh,vv,seq);

            //Altschul's rationale-2 needs distances >=1
            if ((*scale)[J][I]<=0)
                (*scale)[J][I]=1;

            //Copy calculated values to distance matrix
            (*Dij)[J][I] = (*Dij)[I][J] = (*scale)[I][J] = (*scale)[J][I];
        }
    }

    //Free allocated memory for direction matrices
    for (int k = 0; k < num_seq; k++)
    {
        delete dd[k];
        delete hh[k];
        delete vv[k];
    }
    delete[] dd;
    delete[] hh;
    delete[] vv;
}

/*! Recursive function to calculate the cost of path linking two tree leafs
*/
float compute_path_cost_rec(TreeNode *A, TreeNode *B, int * path_length, float *** Dij)
{
    //If left node is an internal node (sequenceNumber == -1), then continue recursion on its children
    if (A->sequenceNumber < 0)
    {
        ++(*path_length);
        return(compute_path_cost_rec(A->leftSon,B,path_length,Dij) + compute_path_cost_rec(A->rightSon,B,path_length,Dij));
    }
    //If right node is an internal node (sequenceNumber == -1), then continue recursion on its children
    else if (B->sequenceNumber < 0)
    {
        ++(*path_length);
        return(compute_path_cost_rec(A,B->leftSon,path_length,Dij) + compute_path_cost_rec(A,B->rightSon,path_length,Dij));
    }
    //If both nodes are leafs, then return weight calculated on primer
    return ((*Dij)[A->sequenceNumber][B->sequenceNumber]); 
}

/*! Calculate path cost between two leafs and divide by path length linking them, identifying leafs by node position on tree list
*/
float compute_path_cost(int i, int j, std::vector<TreeNode*> *tree, float ***Dij)
{
    int path_length = 1;
    float cost = compute_path_cost_rec((*tree)[i],(*tree)[j],&path_length,Dij);

    //Calculates the path cost that links two leafs and divide by the number of nodes linking them (length)
    return (float)(cost/path_length);
}

/*! Calculate path cost between two leafs and divide by path length linking them, identifying leafs by node pointer
*/
float compute_path_cost_n(TreeNode* A, TreeNode* B, std::vector<TreeNode*> *tree, float ***Dij)
{
    int path_length = 1;
    float cost = compute_path_cost_rec(A,B,&path_length,Dij);

    //Calculates the path cost that links two leafs and divide by the number of nodes linking them (length)
    return (float)(cost/path_length);
}

float compute_S(int i, int j, int numNodes, std::vector<TreeNode*> *tree, float ***Dij)
{
	int t, tt;
	float s1=0, s2=0;

	for (t=0;t<numNodes;t++)
		if (t!=i && t!=j)
		{
			s1 += compute_path_cost(i,t,tree,Dij)+compute_path_cost(j,t,tree,Dij);
        }
	s1 = s1 / (2 * (numNodes - 2));

	for (t=0; t<numNodes-1; t++)
		for (tt=t+1; tt<numNodes; tt++)
	    	if (t!=i && t!=j && tt!=i && tt!=j)
	    	{
                s2 += compute_path_cost(t,tt,tree,Dij);
            }
    s2 = s2 / (numNodes- 2);
    float total = (s1 + s2 + compute_path_cost(i,j,tree,Dij) / 2);
    return total;
}



#define BIG_MIN_VAL 1.0E20

/*! Build phylogenetic tree and calculate costs between leafs (sequences that are going to be aligned)
*/
void phylogeneticThreeNeighborJoin(std::vector<std::string>  seq, std::vector<TreeNode*> *tree, std::vector<TreeNode*> *nodesList, float *** Dij)
{
    int num_seq,nodesRemaining;
    float tmp,min = BIG_MIN_VAL;
    int min_i,min_j;


    min_i = min_j = 0;
    //std::cout << "neighbors join 1" << std::endl;
    //Get number of sequences being aligned
    num_seq = nodesRemaining = seq.size();

    //Create leafs and add to tree
    for (int i = 0; i < num_seq; i++)
    {
        TreeNode * node = new TreeNode(nullptr, nullptr, nullptr, nullptr, 0.0, 0.0, 0.0, 0.0, 0.0, i);
        (*tree).push_back(node);
        (*nodesList).push_back(node);
    }

    //std::cout << "neighbors join 2" << std::endl;
    //Build tree until we have only two nodes remaining (start as num_seq)
    while(nodesRemaining > 2)
    {
        //Find best pair of leafs to join with an internal node
        for (int i = 0; i < nodesRemaining-1; i++)
        {
            for (int j = i+1; j < nodesRemaining; j++)
            {
                //Recursive compute of path cost divided by length
                tmp = compute_S(i,j,nodesRemaining,tree,Dij);
           
                //If length smaller than the previous minimum, save the node pair and length
                if (tmp < min)
                {
                    min_i = i;
                    min_j = j;
                    min = tmp;
                }
            }
        }

        //Remove two nodes and join them with an internal node
        join_nodes(min_i,min_j,tree,nodesList,Dij);

        //As we remove 2 nodes from the "tree" list and add an internal one, we need to subtract a node each round
        nodesRemaining--;

        //Assign a big value for min, to prevent problems
        min = BIG_MIN_VAL;
    }

    //Pick up the remaining nodes
    TreeNode *leftNode,*rightNode;
    leftNode = (*tree)[0];
    rightNode = (*tree)[1];

    //Create an ancestor node that connects them
    TreeNode * ancestor = new TreeNode(nullptr, leftNode, rightNode, nullptr, 0.0, 0.0, 0.0, 0.0, 0.0, -2);

    //Update nodes pointers to brothers and parent
    leftNode->brother = rightNode;
    rightNode->brother = leftNode;
    leftNode->parent = rightNode->parent = ancestor;

    //Remove left and right node from node list, keeping only the ancestor node, finishing the tree structure
    (*tree).erase((*tree).begin());
    (*tree).erase((*tree).begin());

    //Emplace the ancestor on the nodes list (the tree one and real list one)
    (*tree).push_back(ancestor);
    (*nodesList).push_back(ancestor);

    //Calculate weight of left son of ancestor
    int count2 = 1;
    float len;
    len = compute_path_cost_n(leftNode,rightNode,tree,Dij);
    len -= compute_path_cost_to_leafs(leftNode,0.0,&count2) / count2;
    count2 = 1;
    len -= compute_path_cost_to_leafs(rightNode,0.0,&count2) / count2;
    ancestor->leftSon->weight = len;



}

void compute_weights_from_tree(float product, float sum, TreeNode* no, TreeNode* brother, float *** weightMatrix, TreeNode ** pN)
{
     if (no->sequenceNumber > INTERNAL_NODE)
    {
        (*weightMatrix)[(*pN)->sequenceNumber][no->sequenceNumber] = sum*product;
    }
    else if (brother==nullptr)
    {
        compute_weights_from_tree(product * no->leftSon->W, sum + no->rightSon->weight, no->rightSon, nullptr, weightMatrix, pN);
        compute_weights_from_tree(product * no->rightSon->W, sum + no->leftSon->weight, no->leftSon, nullptr, weightMatrix, pN);
    }
    else
    {
        compute_weights_from_tree(product * no->V, sum + brother->weight, brother, nullptr, weightMatrix, pN);
        if (no->sequenceNumber != TREE_ROOT)
            compute_weights_from_tree(product * brother->W, sum + no->weight, no->parent, no->brother, weightMatrix, pN);
    }
}

/*! Calculate weights for weighted sum-of-pairs usings Altschul's rationale-2 method
*/
void weightAltschulsRationale2(Sequences * seqs)
{
    //std::cout << "altschul 1" << std::endl;
    float ** weightMatrix, ** weightMatrix2, sm;
    std::vector<TreeNode*> tree;
    std::vector<TreeNode*> nodesList;
    int num_seq = Sequences::get_seq_num();
    std::vector<std::string> sequences;

    //Allocate a matrix to store weights
    weightMatrix = new float*[num_seq]();
    weightMatrix2 = new float*[num_seq]();
    HeuristicHPair * inst = HeuristicHPair::getInstance();
    inst->weightMatrix = new float*[num_seq]();

    std::string dashstr = "-";
    for (int i = 0; i < num_seq; i++)
    {
        weightMatrix[i] = new float[num_seq];
        inst->weightMatrix[i] = new float[num_seq];

        //Adding dash to all sequences to workaround Altschul's algorithm
        sequences.push_back(seqs->get_seq(i));
        sequences[i].insert(0, dashstr, 0, 1); 
        for (int j = 0; j < num_seq; j++)
        {
            weightMatrix[i][j] = 0.0;
            inst->weightMatrix[i][j] = 0.0;
        }
    }

    //Calculate path and costs
    primer(&sequences, &(inst->weightMatrix),&(weightMatrix));

    //Determine optimum sequence pairs (similar sequences) and make them neighbors
    phylogeneticThreeNeighborJoin(sequences,&tree,&nodesList, &(inst->weightMatrix));

    //Print tree
    //rpt(tree[0]);

    TreeNode ** pN, *no;
    //Compute partial weights of all leafs
    for (pN=nodesList.data(); (*pN)->sequenceNumber > INTERNAL_NODE; ++pN)
    {
        no=*pN;
        no->w = 1.0;
        no->W = no->weight;
    }

    //Compute partial weights of all internal nodes
    for (; (no= *pN)->sequenceNumber > TREE_ROOT; ++pN)
    {
        no->w = no->leftSon->w * no->rightSon->W + no->leftSon->W * no->rightSon->w;
        no->W = no->weight  * no->w     + no->leftSon->W * no->rightSon->W;
    }

    //Set root values
    no->V = 1;
    no->v = 0;

    //Continue computing of partial weights, back to the first node
    do
    {
        no= *(--pN);
        no->v = no->parent->v * no->brother->W + no->parent->V * no->brother->w;
        no->V = no->weight   * no->v      + no->parent->V * no->brother->W;
    }
    while (pN != nodesList.data());

    //Compute final Altschul weights based on the tree created
    for(; (no= *pN)->sequenceNumber > INTERNAL_NODE; ++pN)
        compute_weights_from_tree(1.0, no->weight,no->parent,no->brother, &weightMatrix, pN);

    // Scale pair weights so that smallest is about 8
    sm=1.0E+30;

    for (int j=1;j<num_seq;++j)
        for (int i=0;i<j;++i)
            if (weightMatrix[i][j]<sm)
                sm=weightMatrix[i][j];

    sm /= 7.9;

    for (int i=0;i<num_seq-1;++i)
        for (int j=i+1;j<num_seq;++j)
            inst->weightMatrix[i][j]=inst->weightMatrix[j][i]=(weightMatrix[i][j]/sm+0.5);

    //Free intermediary weight matrix
    for (int i=0; i < num_seq;i++)
    {
        delete weightMatrix[i];
    }

    delete[] weightMatrix;

}

/*
void weightGotohTreeway()
{

}

void printPhylogeneticThree()
{

}

void savePhylogeneticThree()
{

}

void printWeight()
{

}

void saveWeight()
{

}
*/
