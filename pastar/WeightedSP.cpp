/*!
 * \authors {Daniel Sundfeld, Gabriel Ferreira}
 * \copyright MIT License
 *  Mostly rewriting Altschul's WSP code and planning to implement Gotoh's alternative
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

void join_nodes(int min_i, int min_j, std::vector<TreeNode*> *tree,std::vector<TreeNode*> *nodesList)
{
    //std::cout<< "join nodes 1" << std::endl;
    TreeNode *leftNode,*rightNode,*newNode;
    leftNode = rightNode = newNode = nullptr;

    //Get pointer to left node and remove from node list
    leftNode = (*tree)[min_i];

    //Get pointer to right node and remove from node list
    rightNode = (*tree)[min_j];

    //Remove references
    (*tree).erase((*tree).begin()+min_i);

    if (min_i > min_j)
        (*tree).erase((*tree).begin()+min_j);
    else
        (*tree).erase((*tree).begin()+min_j-1);

    //std::cout<< "join nodes 2" << std::endl;
    //Create new internal node
    newNode = new TreeNode(nullptr,leftNode,rightNode,nullptr,leftNode->weight*rightNode->weight,0.0, 0.0, 0.0, 0.0, -1);

    //Update previously removed nodes
    leftNode->brother = rightNode;
    rightNode->brother = leftNode;
    leftNode->parent = rightNode->parent = newNode;

    //Create new node and emplace it on the list
    (*tree).push_back(newNode);
    (*nodesList).push_back(newNode);
    //std::cout<< "join nodes 3" << std::endl;

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
        std::cout << "I="<< I << " J=" << J;
        std::cout << " V=" << V << " H="<<H<<" M="<<M<<std::endl;
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
    register int    I,J,i,j,Gi,Gj,n,m;
    //std::cout << "primer 1" << std::endl;
    int num_seq, max_seq_length = 1000;
    num_seq = seq->size();

    if (num_seq == 0)
        exit(-1);

    std::string seqA, seqB;

    //Allocate memory for diagonal, horizontal and vertical matrices
    int **dd = nullptr, **hh = nullptr,**vv= nullptr;


    dd = new int*[max_seq_length]();
    hh = new int*[max_seq_length]();
    vv = new int*[max_seq_length]();

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
        }
    }


    //std::cout << "primer 2" << std::endl;
    //std::cout << "num_seq=" << num_seq << std::endl;

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

            //std::cout << "primer 3" << std::endl;
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

            //std::cout << "primer 4" << std::endl;

            //For each element in sequences I and J
            for (i=1; i<=n; i++)
            {
                for (Gi=i==n?EfectiveGapCost:GapCost, j=1; j<=m; j++)
                {
                    //todo: discover what that does
                    Gj= j==m ? EfectiveGapCost : GapCost;

                    //dafuck
                    dd[i][j] = minOf3(dd[i-1][j-1] , hh[i-1][j-1] , vv[i-1][j-1] ) + Cost::cost(seqA[i], seqB[j]);

                    hh[i][j] = minOf3(dd[i][j-1]+Gi, hh[i][j-1]   , vv[i][j-1]+Gi) + Cost::cost(DASH   , seqB[j]);

                    vv[i][j] = minOf3(dd[i-1][j]+Gj, hh[i-1][j]+Gj, vv[i-1][j]   ) + Cost::cost(seqA[i], DASH   );
                }
            }

            //std::cout << "primer 5" << std::endl;

            //Calculate path cost based on the path
            (*scale)[J][I] = convert_path_to_cost(I,J,n-1,m-1,dd,hh,vv,seq);

            //std::cout << "primer 6" << std::endl;

            //Altschul's rationale-2 needs distances >=1
            if ((*scale)[J][I]<=0)
                (*scale)[J][I]=1;

            //Copy calculated values to distance matrix
            (*Dij)[J][I] = (*Dij)[I][J] = (*scale)[J][I];
        }
    }
    //std::cout << "primer 7" << std::endl;

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
    //std::cout << "comp path rec 1" << std::endl;
    if (A->sequenceNumber < 0)
    {
        ++(*path_length);
        return(compute_path_cost_rec(A->leftSon,B,path_length,Dij) + compute_path_cost_rec(A->rightSon,B,path_length,Dij));
    }
    //std::cout << "comp path rec 2" << std::endl;
    else if (B->sequenceNumber < 0)
    {
        ++(*path_length);
        return(compute_path_cost_rec(A,B->leftSon,path_length,Dij) + compute_path_cost_rec(A,B->rightSon,path_length,Dij));
    }
    //std::cout << "rdist["<<A->sequenceNumber<<"]["<<B->sequenceNumber<<"]="<<((*Dij)[A->sequenceNumber][B->sequenceNumber]) << std::endl;
    return ((*Dij)[A->sequenceNumber][B->sequenceNumber]); // TODO: discover how to implement that thing
}

/*! Calculate path cost between two leafs and divide by path length linking them
*/
float compute_path_cost(int i, int j, std::vector<TreeNode*> *tree, float ***Dij)
{
    int path_length = 1;
    //std::cout << "comp path cost 1" << std::endl;
    //std::cout << "i=" << i << " j=" << j << " tree[i]=" << (*tree)[i] << " tree[j]=" << (*tree)[j] << std::endl;
    float cost = compute_path_cost_rec((*tree)[i],(*tree)[j],&path_length,Dij);

    //std::cout << "dist["<<i<<"]["<<j<<"="<<cost<< " count=" <<path_length << std::endl;
    //Calculates the path cost that links two leafs and divide by the number of nodes linking them (length)
    return (float)(cost/path_length);
}

float compute_S(int i, int j, std::vector<TreeNode*> *tree, float ***Dij)
{
	int t, tt, numNodes = (*tree).size();
	float s1=0, s2=0;

	for (t=0;t<numNodes;t++)
		if (t!=i && t!=j)
		{
			s1 += compute_path_cost(i,t,tree,Dij)+compute_path_cost(j,t,tree,Dij);
			//std::cout << "s1=" << s1 << std::endl;
        }
	s1 = s1 / (2 * (numNodes - 2));
	//std::cout << "s1 final=" << s1 << std::endl;

	for (t=0; t<numNodes; t++)
		for (tt=t+1; tt<numNodes; tt++)
	    	if (t!=i && t!=j && tt!=i && tt!=j)
	    	{
                s2 += compute_path_cost(t,tt,tree,Dij);
                //std::cout << "s2=" << s2 << std::endl;
            }
    s2 = s2 / (numNodes- 2);
    //std::cout << "s2 final=" << s2 << std::endl;
    float total = (s1 + s2 + compute_path_cost(i,j,tree,Dij) / 2);
    return total;
}

float compute_path_cost_to_leafs(TreeNode * A, float total, int * count2)
{
    if (A->sequenceNumber >= 0)
        return(total + A->weight);

    (*count2)++;
    return(compute_path_cost_to_leafs(A->leftSon,A->weight+total,count2) + compute_path_cost_to_leafs(A->rightSon,A->weight+total,count2));
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
    std::cout << "neighbors join 1" << std::endl;
    //Get number of sequences being aligned
    num_seq = nodesRemaining = seq.size();

    //Create leafs and add to tree
    for (int i = 0; i < num_seq; i++)
    {
        TreeNode * node = new TreeNode(nullptr, nullptr, nullptr, nullptr, 0.0, 0.0, 0.0, 0.0, 0.0, i);
        (*tree).push_back(node);
        (*nodesList).push_back(node);
    }



    std::cout << "neighbors join 2" << std::endl;
    //Build tree until we have only two nodes remaining (start as num_seq)
    while(nodesRemaining > 2)
    {
        //Find best pair of leafs to join with an internal node
        for (int i = 0; i < nodesRemaining-1; i++)
        {
            for (int j = i+1; j < nodesRemaining; j++)
            {
                //Recursive compute of path cost divided by length
                tmp = compute_S(i,j,tree,Dij);
                std::cout << "minimize_Sij[" << i << "][" << j << "] tmp=" <<tmp<< " min=" <<min << " Dij="<<(*Dij)[i][j]<<std::endl;
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
        join_nodes(min_i,min_j,tree,nodesList);

        //As we remove 2 nodes from the "tree" list and add an internal one, we need to subtract a node each round
        nodesRemaining--;

        //Assign a big value for min, to prevent problems
        min = BIG_MIN_VAL;
    }
    std::cout << "neighbors join 3" << std::endl;

    //Join last remaining nodes with a root node
    TreeNode * ancestor = new TreeNode(nullptr, (*tree)[0], (*tree)[1], nullptr, 0.0, 0.0, 0.0, 0.0, 0.0, -2);
    int count2 = 1;
    float len;
    len = compute_path_cost(0,1,tree,Dij);
    len -= compute_path_cost_to_leafs((*tree)[0],0.0,&count2) / count2;
    count2 = 1;
    len -= compute_path_cost_to_leafs((*tree)[1],0.0,&count2) / count2;
    ancestor->weight = len;

    TreeNode *leftNode,*rightNode;
    leftNode = rightNode = nullptr;

    //Get pointer to left node and remove from node list
    leftNode = (*tree).at(0);
    (*tree).erase((*tree).begin());

    //Get pointer to right node and remove from node list
    rightNode = (*tree).at(0);
    (*tree).erase((*tree).begin());

    //Update previously removed nodes
    leftNode->brother = rightNode;
    rightNode->brother = leftNode;
    leftNode->parent = rightNode->parent = ancestor;

    //Create new node and emplace it on the list
    (*tree).push_back(ancestor);
    (*nodesList).push_back(ancestor);
    std::cout << "neighbors join 4" << std::endl;

}

void compute_weights_from_tree(float product, float sum, TreeNode* no, TreeNode* brother, float *** weightMatrix, TreeNode ** pN)
{
    printf("trace\nno=%d\nro=%d\npN=%d\nprod=%f\nsum=%f\n\n\n", no, brother, *pN,product, sum);
    if (no->sequenceNumber > -1)
    {
        (*weightMatrix)[(*pN)->sequenceNumber][no->sequenceNumber] = sum*product;
    }
    else if (brother==NULL)
    {
        compute_weights_from_tree(product * no->leftSon->W, sum + no->rightSon->weight, no->rightSon, NULL, weightMatrix, pN);
        compute_weights_from_tree(product * no->rightSon->W, sum + no->leftSon->weight, no->leftSon, NULL, weightMatrix, pN);
    }
    else
    {
        compute_weights_from_tree(product * no->V, sum + brother->weight, brother, NULL, weightMatrix, pN);
        if (no->sequenceNumber != -1)
            compute_weights_from_tree(product * brother->W, sum + no->weight, no->parent, no->brother, weightMatrix, pN);
    }
}

/*! Calculate weights for weighted sum-of-pairs usings Altschul's rationale-2 method
*/
void weightAltschulsRationale2(Sequences * seqs)
{
    std::cout << "altschul 1" << std::endl;
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
        sequences[i].insert(0,dashstr,0,1);
        for (int j = 0; j < num_seq; j++)
        {
            weightMatrix[i][j] = 0.0;
            inst->weightMatrix[i][j] = 0.0;
        }
    }
    //std::cout << "altschul 2" << std::endl;

    //Calculate path and costs
    primer(&sequences, &(inst->weightMatrix),&(weightMatrix));
    //std::cout << "altschul 3" << std::endl;

    //Determine optimum sequence pairs (similar sequences) and make them neighbors
    phylogeneticThreeNeighborJoin(sequences,&tree,&nodesList, &(inst->weightMatrix));
    //std::cout << "altschul 4" << std::endl;

    rpt(tree[0]);

    TreeNode ** pN, *no;
    //Compute partial weights
    for (pN=nodesList.data(); (no= *pN)->sequenceNumber > -1; ++pN)
    {
        no->w = 1.0;
        no->W = no->weight;
    }
    for (; (no= *pN)->sequenceNumber > -2; ++pN)
    {
        no->w = no->leftSon->w * no->rightSon->W + no->leftSon->W * no->rightSon->w;
        no->W = no->weight  * no->w     + no->leftSon->W * no->rightSon->W;
    }
    no->V = 1;
    no->v = 0;
    do
    {
        no= *(--pN);
        no->v = no->parent->v * no->brother->W + no->parent->V * no->brother->w;
        no->V = no->weight   * no->v      + no->parent->V * no->brother->W;
    }
    while (pN != nodesList.data());

    //Compute final Altschul weights based on the tree created
    for(; (no= *pN)->sequenceNumber > -1; ++pN)
        compute_weights_from_tree(1.0, no->weight,no->parent,no->brother, &weightMatrix, pN);

    std::cout << "altschul 5" << std::endl;
    // Scale pair weights so that smallest is about 8
    sm=1.0E+30;

    for (int j=1;j<num_seq;++j)
        for (int i=0;i<j;++i)
            if (weightMatrix[i][j]<sm)
                sm=weightMatrix[i][j];
    sm /= 7.9;
    for (int i=0;i<num_seq-1;++i)
        for (int j=i+1;j<num_seq;++j)
        {
            inst->weightMatrix[i][j]=weightMatrix[i][j]/sm+0.5;
            std::cout << "finalWeightMatrix["<<i<<"]["<< j << "]="<< inst->weightMatrix[i][j];
            std::cout << " weightMatrix[" << i << "][" << j<<"]=" << weightMatrix[i][j] << std::endl;
        }

    std::cout << "altschul 6" << std::endl;
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
