#pragma once

#include "Tree.h"


class SolutionNode {
public:

    int leftDTid;
    int rightDTid;

    SolutionNode* left, *right;

    Tree* tree;
    
    SolutionNode(){
        this->left = NULL;
        this->right = NULL;
        this->tree = NULL;
        this->leftDTid = -1;
        this->rightDTid = -1;
    }
};

class SolutionTree {
public:

    static Tree* createTree(SolutionNode* root){
        return dfs(root);
    }

    static Tree* dfs(SolutionNode* solnNode){
        if(solnNode->tree != NULL){
            return solnNode->tree;
        }
        else{
            Tree** trees = new Tree*[2];
            trees[0] = dfs(solnNode->left);
            trees[1] = dfs(solnNode->right);

            // int[] dummyIds = new int[2];
            int dummyIds[2];

            dummyIds[0] = solnNode->leftDTid;
            dummyIds[1] = solnNode->rightDTid;

            TreeNode** dtNodes = new TreeNode*[2];
            dtNodes[0] = NULL;
            dtNodes[1] = NULL;


            for(int i = 0; i < 2; ++i){
                for(auto node : trees[i]->nodes){
                    if(node->dummyTaxonId == dummyIds[i]){
                        dtNodes[i] = node;
                    }
                }
                if(dtNodes[i] == NULL){
                    std::cout << "Error: Dummy node not found" << "\n";
                    exit(-1);
                }
                if(dtNodes[i]->childs != NULL){
                    cout << "Error: Dummy node should be leaf" << "\n";
                    exit(-1);
                }
            }
            
            trees[0]->reRootTree(dtNodes[0]);
            if(dtNodes[0]->childs->size() > 1){
                // System.out.println("Error: Dummy node should have only one child after reroot");
                // System.exit(-1);
                cout << "Error: Dummy node should have only one child after reroot" << "\n";
                exit(-1);
            }
            // dtNodes[0]->childs->at(0)->setParent(dtNodes[1]->parent);
            dtNodes[0]->childs->at(0)->parent = dtNodes[1]->parent;
            dtNodes[1]->parent->childs->erase(
                find(dtNodes[1]->parent->childs->begin(), dtNodes[1]->parent->childs->end(), dtNodes[1])
            );
            dtNodes[1]->parent->childs->push_back(dtNodes[0]->childs->at(0));
            for(auto node: trees[0]->nodes){
                trees[1]->nodes.push_back(node);
            }
    
            return trees[1];
            
        }
    }

    
};
