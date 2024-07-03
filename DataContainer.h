#pragma once
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include "Taxon.h"
#include "NumSat.h"

using namespace std;



class PartitionByTreeNode;
string getPartitionString(bool* b, int n);

class PartitionByTreeNodeWithIndex{
public:
    PartitionByTreeNode* partitionByTreeNode;
    int index;
    
    PartitionByTreeNodeWithIndex(PartitionByTreeNode* a, int b){
        this->partitionByTreeNode = a;
        this->index = b;
    }
};


class PartitionNode{

public:
    vector<PartitionNode*>* parents;
    vector<PartitionNode*>* children;

    vector<PartitionByTreeNodeWithIndex*>* nodePartitions;
    
    bool isLeaf;
    string label;

    Data** data;

    PartitionNode(bool isLeaf){
        this->parents = new vector<PartitionNode*>();
        this->children = new vector<PartitionNode*>();
        this->isLeaf = isLeaf;
        this->nodePartitions = new vector<PartitionByTreeNodeWithIndex*>();
        // this->data = new Data*[N_THREADS];

        // for(int i = 0; i < N_THREADS; ++i){
        //     this->data[i] = new Data();
        //     this->data[i]->gainsForSubTree = new double[2];
        //     this->data[i]->branch = new Branch(0);
        
        // }
    }
    
    void addChild(PartitionNode* child){
        this->children->push_back(child);
    }

    void addParent(PartitionNode* parent){
        this->parents->push_back(parent);
    }

    void addNodePartitions(PartitionByTreeNode *p, int index){
        this->nodePartitions->push_back(new PartitionByTreeNodeWithIndex(p, index));
    }

};

class PartitionByTreeNode{
public:
    vector<PartitionNode*>* partitionNodes;

    int count;
    
    NumSatCalculatorNode** scoreCalculator;

    PartitionByTreeNode(vector<PartitionNode*>* partitionNodes){
        this->partitionNodes = partitionNodes;
        this->count = 1;
        for(int i = 0; i < partitionNodes->size(); ++i){
            PartitionNode* p = (*partitionNodes)[i];
            p->addNodePartitions(this, i);
        }
        // this->scoreCalculator = new NumSatCalculatorNode*[N_THREADS];
        // for(int i = 0; i < N_THREADS; ++i){
        //     vector<Branch*>* b = new vector<Branch*>(this->partitionNodes->size());
        //     for(int j = 0; j < this->partitionNodes->size(); ++j){
        //         b->at(j) = this->partitionNodes->at(j)->data[i]->branch;
        //     }
        //     if(this->partitionNodes->size() > 3){
        //         this->scoreCalculator[i] = new NumSatCalculatorNodeE(b, new int[0]);
        //     }
        //     else{
        //         this->scoreCalculator[i] = new NumSatCalculatorBinaryNode(b, new int[0]);
        //     }
        // }
    }

    void increaseCount(){
        this->count++;
    }

    void batchTransfer(int i, int netTransfer, int tid){
        if(netTransfer != 0){
            this->scoreCalculator[tid]->batchTransferRealTaxon(i, netTransfer);
        }
    }
    
};






class PartitionsByTreeNode {
public:
    vector<PartitionByTreeNode*>* partitions;
    unordered_map<PartitionNode*, bool*>* realTaxaInPartition;
    unordered_map<string, PartitionByTreeNode*>* stringIdToPartitionByTreeNode;

    PartitionsByTreeNode(
        unordered_map<PartitionNode*, bool*>* realTaxaInPartition
    ){
        this->partitions = new vector<PartitionByTreeNode*>();
        this->realTaxaInPartition = realTaxaInPartition;
        this->stringIdToPartitionByTreeNode = new unordered_map<string, PartitionByTreeNode*>();
    }

    void addPartitionByTreeNode(vector<PartitionNode*>* partitionNodes, int nTaxa){
        vector<string> partitionStrings;        
        for(auto p : *partitionNodes){
            partitionStrings.push_back(getPartitionString(this->realTaxaInPartition->at(p), nTaxa));
        }
        sort(partitionStrings.begin(), partitionStrings.end());
        string sb = "";
        for(auto s : partitionStrings){
            sb.append(s);
        }
        string partitionByTreeNodeString = sb;
        auto f = this->stringIdToPartitionByTreeNode->find(partitionByTreeNodeString);
        if(
            f != this->stringIdToPartitionByTreeNode->end()
        ){
            f->second->increaseCount();
        }
        else{
            PartitionByTreeNode* partitionByTreeNode = new PartitionByTreeNode(partitionNodes);
            this->partitions->push_back(partitionByTreeNode);
            this->stringIdToPartitionByTreeNode->insert({partitionByTreeNodeString, partitionByTreeNode});
        }
    }

    int getPartitionCount(){
        return this->partitions->size();
    }
    
    
};


class DataContainer {
public:

    vector<PartitionByTreeNode*>* partitionsByTreeNodes;
    vector<PartitionNode*>* topSortedPartitionNodes;
    vector<PartitionNode*>* realTaxaPartitionNodes;
    bool** realTaxaInTrees;
    vector<RealTaxon*>* taxa;
    int nTaxa;
    int nTrees;
};


class PartitionGraph {
public:

    vector<RealTaxon*>* taxa;
    vector<PartitionNode*>* taxaPartitionNodes;
    unordered_map<PartitionNode*, bool*>* realTaxaInPartition;
    unordered_map<string, PartitionNode*>* stringIdToPartition;
    vector<PartitionNode*>* partitionNodes;
    int count = 0;
    int nTaxa;
    


    PartitionGraph(vector<RealTaxon*>* taxa, int nTaxa){
        this->taxa = taxa;
        this->nTaxa = nTaxa;
        this->taxaPartitionNodes = new vector<PartitionNode*>();
        this->taxaPartitionNodes->resize(nTaxa);

        this->realTaxaInPartition = new unordered_map<PartitionNode*, bool*>();
        this->stringIdToPartition = new unordered_map<string, PartitionNode*>();
        this->partitionNodes = new vector<PartitionNode*>();


        for(int i = 0; i < nTaxa; ++i){
            this->taxaPartitionNodes->at(i) = new PartitionNode(true);
            bool* realTaxaInSubTree = new bool[nTaxa];
            for(int j = 0; j < nTaxa; ++j){
                realTaxaInSubTree[j] = false;
            }
            realTaxaInSubTree[i] = true;
            this->realTaxaInPartition->insert({this->taxaPartitionNodes->at(i), realTaxaInSubTree});
            this->stringIdToPartition->insert({getPartitionString(realTaxaInSubTree, nTaxa), this->taxaPartitionNodes->at(i)});
            this->partitionNodes->push_back(this->taxaPartitionNodes->at(i));
            this->taxaPartitionNodes->at(i)->label = taxa->at(i)->label;
        }

        count = nTaxa;

    }
    
    static string getPartitionString(bool* b, int n){
        // StringBuilder sb = new StringBuilder();
        string sb = "";
        for(int i = 0; i < n; ++i){
            if(b[i]){
                sb.append("1");
            }
            else{
                sb.append("0");
            }
        }
        return sb;
    }

    PartitionNode* getPartitionNode(RealTaxon* taxon){
        return this->taxaPartitionNodes->at(taxon->id);
    }

    PartitionNode* addPartition(vector<PartitionNode*>* childs){
        bool* b = new bool[this->nTaxa];
        for(int i = 0; i < this->nTaxa; ++i){
            b[i] = false;
        }

        for(auto child : *childs){
            bool* realTaxaInSubTree = this->realTaxaInPartition->at(child);
            for(int i = 0; i < this->nTaxa; ++i){
                b[i] = b[i] || realTaxaInSubTree[i];
            }
        }

        string partitionString = getPartitionString(b, this->nTaxa);
        auto f = this->stringIdToPartition->find(partitionString);

        // cout << "PartitionString: " << partitionString << "\n";

        if(
            f != this->stringIdToPartition->end()
        ){
            return f->second;
        }
        else{
            PartitionNode* partitionNode = new PartitionNode(false);
            for(auto child: *childs){
                partitionNode->addChild(child);
                child->addParent(partitionNode);
            }
            this->realTaxaInPartition->insert({partitionNode, b});
            this->stringIdToPartition->insert({partitionString, partitionNode});
            this->partitionNodes->push_back(partitionNode);
            count += 1;
            return partitionNode;
        }
        
    }

    vector<PartitionNode*>* getTopSortedNodes(){
        vector<PartitionNode*>* topSortedNodes = new vector<PartitionNode*>();
        queue<PartitionNode*> q;
        unordered_map<PartitionNode*, int> inDegree;
        
        for(auto partitionNode: *this->partitionNodes){
            if(partitionNode->parents->size() == 0){
                q.push(partitionNode);
            }
            else{
                inDegree.insert({partitionNode, partitionNode->parents->size()});
            }
        }

        while(!q.empty()){
            PartitionNode* partitionNode = q.front();
            q.pop();
            topSortedNodes->push_back(partitionNode);
            for(auto child: *partitionNode->children){
                inDegree[child] -= 1;
                if(inDegree[child] == 0){
                    q.push(child);
                }
            }
        }
        return topSortedNodes;
    }

};


string getPartitionString(bool* b, int n){
    // StringBuilder sb = new StringBuilder();
    string sb = "";
    for(int i = 0; i < n; ++i){
        if(b[i]){
            sb.append("1");
        }
        else{
            sb.append("0");
        }
    }
    return sb;
}
