#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <set>
#include <stack>
#include <algorithm>

#include "string"
#include "Taxon.h"
#include "DataContainer.h"

using namespace std;

class TreeNode{
public:
    int index;
    vector<TreeNode*>* childs;
    TreeNode* parent;
    RealTaxon* taxon;

    PartitionNode* partitionNode, *parentPartitionNode;
    
    
    TreeNode(){
        this->index = -1;
        this->parent = NULL;
        this->childs = NULL;
        this->taxon = NULL;
        this->partitionNode = NULL;
        this->parentPartitionNode = NULL;
    }

    bool isRoot(){
        return this->parent == NULL;
    }

    bool isLeaf(){
        return this->childs == NULL;
    }
    

};


class Tree{
private:

    


public:

    vector<TreeNode*> nodes;
    vector<TreeNode*> topSortedNodes;
    TreeNode* root;
    map<string, RealTaxon*>* taxaMap;
    vector<TreeNode*> leaves;


    int leavesCount;

    ~Tree(){
        for(auto x : nodes){
            delete x;
        }
    }

    Tree(const string& line, map<string, RealTaxon*>* realTaxons){
        this->leavesCount = 0;
        this->taxaMap = realTaxons;
        parseFromNewick(line);
    }

    TreeNode* addNode(vector<TreeNode*>* children, TreeNode* parent){
        TreeNode* nd = new TreeNode();
        nd->index = nodes.size();
        nd->childs = children;
        nd->parent = parent;
        nodes.push_back(nd);
        return nd;
    }


    TreeNode* addInternalNode(vector<TreeNode*>* children){
        TreeNode* nd = addNode(children,NULL);
        for (auto x : *children)
            x->parent = nd;
        return nd;
    }

    TreeNode* addLeaf(RealTaxon* taxon){
        TreeNode* nd = addNode(NULL, NULL);
        nd->taxon = taxon;
        return nd;
    }



    void parseFromNewick(const string& line){
        int leavesCount = 0;

        stack<TreeNode*> nodeStack;

        int n =  line.length();
    
        int i = 0, j = 0;
    
        while(i < n){
            char curr = line[i];
            if(isspace(curr)){
                continue;
            }

            if(curr == '('){
                nodeStack.push(NULL);
            }
            else if(curr == ')'){
                vector<TreeNode*>* arr = new vector<TreeNode*>();
                while( !nodeStack.empty() && nodeStack.top() != NULL){
                    arr->push_back(nodeStack.top());
                    nodeStack.pop();
                }
                if(!nodeStack.empty())
                    nodeStack.pop();
                nodeStack.push(addInternalNode(arr));
                
            }
            else if(curr == ',' || curr == ';'){
    
            }
            else{
                string taxa;
                j = i;
                TreeNode* newNode = NULL;
                while(j < n){
                    char curr_j = line[j];
                    if(curr_j == ')' || curr_j == ','){
                        RealTaxon* taxon;
                        taxon = this->taxaMap->at(taxa);
                        newNode = addLeaf(taxon);
                        leavesCount++;

                        break;
                    }
                    taxa += curr_j;
                    ++j;
                }
                if(j == n){
                    leavesCount++;
                    RealTaxon* taxon;
                    taxon = this->taxaMap->at(taxa);
                    newNode = addLeaf(taxon);
                }
                i = j - 1;
                nodeStack.push(newNode);
            }
            ++i;
        }

        this->leavesCount = leavesCount;
        this->root = nodeStack.top();
        
        if(root->childs->size() > 2){
            balanceRoot();
        }
        filterLeaves();
        topSort();

        // cout << "top sorted nodes count : " << this->topSortedNodes.size() << endl;
        // cout << "leaves count : " << this->leavesCount << endl;

        
    }

    void topSort(){
        topSortUtil(root, this->topSortedNodes);
    }

    void topSortUtil(TreeNode* node, vector<TreeNode*>& topSort){
        if(node->childs == NULL){

        }
        else{
            for(auto x : *node->childs){
                topSortUtil(x, topSort);
            }
        }
        topSort.push_back(node);
    }
    

    void filterLeaves(){
        this->leaves.resize(this->leavesCount, NULL);
        for(auto x : nodes){
            if(x->childs == NULL){
                this->leaves[x->taxon->id] = x;
            }
        }
    }

    int dfs(TreeNode* node, vector<int>& subTreeNodeCount){
        if(node->childs == NULL){
            subTreeNodeCount[node->index] = 1;
            return 1;
        }
        int res = 0;
        for(auto x : *(node->childs)){
            res += dfs(x, subTreeNodeCount);
        }
        subTreeNodeCount[node->index] = res + 1;
        return res + 1;
    }

    void balanceRoot(){

        int n = nodes.size();
        vector<int> subTreeNodeCount(n, 0);

        dfs(root, subTreeNodeCount);
        
        TreeNode* closest = root;
        int diff = n;
        int v;
        for(int i = 0; i < n; ++i){
            v = abs(n/2 - subTreeNodeCount[i]); 
            if(v < diff){
                diff = v;
                closest = nodes[i];
            }
        }
        TreeNode* closestP = closest->parent;
        closestP->childs->erase(
            find(closestP->childs->begin(), closestP->childs->end(), closest)
        );

        TreeNode* curr = closestP;
        TreeNode* currP, *temp;
        currP = curr->parent;
        while(curr != NULL && currP != NULL){
            curr->childs->push_back(currP);
            currP->childs->erase(
                find(currP->childs->begin(), currP->childs->end(), curr)
            );
            temp = currP;
            currP = currP->parent;
            temp->parent = curr;
            curr = temp;
        }

        vector<TreeNode*>* arr = new vector<TreeNode*>();
        arr->push_back(closest);
        arr->push_back(closestP);

        root = addInternalNode(arr) ;

    }

    string newickFormatUitl(TreeNode* node){
        if(node->childs == NULL){
            return node->taxon->label;
        }
        string s;
        s += "(";
        for(int i = 0; i < node->childs->size(); ++i){
            s += newickFormatUitl(node->childs->at(i));
            if(i != node->childs->size() - 1){
                s += ",";
            }
        }
        s += ")";
        return s;
    }


    string getNewickFormat(){
        return (newickFormatUitl(root) + ";");
    }
    
    

};

class GeneTrees{
public:

    vector<Tree*> trees;
    string path;
    int realTaxonCount;
    RealTaxon** taxa;    
    map<string, RealTaxon*>* realTaxons;

    GeneTrees(string path){
        this->path = path;
        this->realTaxonCount = 0;
    }

    void parseTaxa(const string& line, set<string>& taxaNames){
        int n = line.length();
        int i = 0, j = 0;
        while(i < n){
            char curr = line[i];
            if(isspace(curr)){
                continue;
            }
            if(curr == '(' || curr == ')' || curr == ',' || curr == ';'){

            }
            else{
                string taxa;
                j = i;
                while(j < n){
                    char curr_j = line[j];
                    if(curr_j == ')' || curr_j == ','){
                        taxaNames.insert(taxa);
                        break;
                    }
                    taxa += curr_j;
                    ++j;
                }
                if(j == n){
                    taxaNames.insert(taxa);
                }
                i = j - 1;
            }
            ++i;
        }
    }

    void readTaxaNames(){
        this->realTaxons = new map<string, RealTaxon*>();
        set<string> taxaNames;

        ifstream file(path);
        
        if(!file.is_open()){
            cout << "Gene Trees File not found" << endl;
            exit(1);
        }
        string line;
        while(getline(file, line)){
            if(line.length() == 0){
                continue;
            }
            parseTaxa(line, taxaNames);
        }
        
        file.close();

        this->taxa = new RealTaxon*[taxaNames.size()];

        for(auto it = taxaNames.begin(); it != taxaNames.end(); ++it){
            RealTaxon* realTaxon = new RealTaxon(*it);
            (*this->realTaxons)[*it] = realTaxon;
            this->taxa[realTaxon->id] = realTaxon;
        }

        this->realTaxonCount = taxaNames.size();

        // print contents of realTaxons
        

    }


    void readGeneTrees(){
        ifstream file(path);
        
        if(!file.is_open()){
            cout << "Gene Trees File not found" << endl;
            exit(1);
        }
        string line;
        while(getline(file, line)){
            if(line.length() == 0){
                continue;
            }
            // cout << line << endl;
            // for(auto it = this->realTaxons->begin(); it != this->realTaxons->end(); ++it){
            //     cout << it->first << " len: " << it->first.length() << endl;
            // }
            auto tree = new Tree(line, this->realTaxons);
            this->trees.push_back(tree);
            // cout << tree->getNewickFormat() << endl;
        }
        
        file.close();
    }

    DataContainer createDataContainer(){
        auto partitionGraph = new PartitionGraph(this->taxa, this->realTaxonCount);
        // cout << "realTaxonCount : " << this->realTaxonCount << endl;
        for(auto tree : this->trees){
            // cout << tree->getNewickFormat() << endl;

            for(auto node : tree->leaves){
                node->partitionNode = partitionGraph->getPartitionNode(node->taxon);
            }

            for (auto node : tree->topSortedNodes) {
                if(node->isRoot() || node->isLeaf()) continue;
                // ArrayList<PartitionNode> childs = new ArrayList<>();
                vector<PartitionNode*>* childs = new vector<PartitionNode*>();
                for(auto child : *node->childs){
                    childs->push_back(child->partitionNode);
                }
                node->partitionNode = partitionGraph->addPartition(childs);
            }
            tree->root->childs->at(0)->parentPartitionNode = tree->root->childs->at(1)->partitionNode;
            tree->root->childs->at(1)->parentPartitionNode = tree->root->childs->at(0)->partitionNode;
            
            int sz = tree->topSortedNodes.size() - 1;
            for(int i = sz - 1; i > -1; --i){
                auto node = tree->topSortedNodes[i];
                if(node->isLeaf() || node->parent == tree->root) continue;
                vector<PartitionNode*>* childs = new vector<PartitionNode*>();

                for(auto child : *node->parent->childs){
                    if(child == node) continue;
                    childs->push_back(child->partitionNode);
                }
                childs->push_back(node->parent->parentPartitionNode);
                node->parentPartitionNode = partitionGraph->addPartition(childs);
            }

        }


        auto partitions = new PartitionsByTreeNode(partitionGraph->realTaxaInPartition);
        for(auto tree : this->trees){
            for(auto node : tree->topSortedNodes){
                if(node->isLeaf() || node->isRoot()) continue;
                vector<PartitionNode*>* ps = new vector<PartitionNode*>();
                for(auto child : *node->childs){
                    ps->push_back(child->partitionNode);
                }
                ps->push_back(node->parentPartitionNode);
                partitions->addPartitionByTreeNode(ps,this->realTaxonCount);
            }
        }


        cout << "Partition graph created" << endl;
        cout << "Partition graph nodes count : " << partitionGraph->count << endl;

        cout << "Partitions created" << endl;
        cout << "Partitions count : " << partitions->getPartitionCount() << endl;
        

        DataContainer dc;

        dc.partitionsByTreeNodes = partitions->partitions;
        dc.topSortedPartitionNodes = partitionGraph->getTopSortedNodes();
        dc.realTaxaPartitionNodes = partitionGraph->taxaPartitionNodes;
        dc.taxa = this->taxa;
        dc.nTaxa = this->realTaxonCount;
        dc.nTrees = this->trees.size();

        dc.realTaxaInTrees = new bool*[this->trees.size()];
        
        for(int i = 0; i < this->trees.size(); ++i){
            dc.realTaxaInTrees[i] = new bool[this->realTaxonCount];
            for(int j = 0; j < this->realTaxonCount; ++j){
                dc.realTaxaInTrees[i][j] = this->trees[i]->leaves[j] != NULL;
            }
        }

        return dc;

    }

    ~GeneTrees(){
        for(auto tree : this->trees){
            delete tree;
        }
        delete this->realTaxons;
    }


};

