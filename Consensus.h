#pragma once
#include <vector>
#include <string>

#include "Taxon.h"
#include "Tree.h"
#include "MakePartition.h"
#include "PerLevelDS.h"


class RandomPartition : public IMakePartition{
public:
    virtual MakePartitionReturnType* makePartition(
        vector<RealTaxon*>* rts,
        vector<DummyTaxon*>* dts, 
        int tid
    ){
        // random partition
        int pa = (rts->size() + dts->size()) / 2;

        int* rtPart = new int[rts->size()];
        for(int i = 0; i < rts->size(); ++i){
            if(pa > 0){
                rtPart[i] = 0;
                --pa;
            }
            else{
                rtPart[i] = 1;
            }
        }
        auto dtPart = new int[dts->size()];
        for(int i = 0; i < dts->size(); ++i){
            if(pa > 0){
                dtPart[i] = 0;
                --pa;
            }
            else{
                dtPart[i] = 1;
            }
        }

        return new MakePartitionReturnType(rtPart, dtPart);
    }
};

class CarryInfo{
public:
    double* weight;
    bool* isRealTaxon;
    int* inWhichDummyTaxa;
    vector<RealTaxon*>* rts;
    vector<DummyTaxon*>* dts;
    BookKeepingPerLevel* book;
    TaxaPerLevelWithPartition* taxas;
    TreeNode* minNode;
    double maxScore;
    double minDiff;
    double score;
    
    double* minNodeDummyTaxaWeights;

    int tid;

    ~CarryInfo(){
        delete[] weight;
        delete[] isRealTaxon;
        delete[] inWhichDummyTaxa;
        delete[] minNodeDummyTaxaWeights;

        if(book != NULL){
            delete taxas;
            delete book;
        }
        
    }
};


class ConsensusTreePartition : public IMakePartition{
public:
    Tree* consTree;
    DataContainer* dc;
    RandomPartition* rp;
    int taxonCount;

    ConsensusTreePartition(string file, map<string, RealTaxon*>* taxaMap, DataContainer* dc){
        this->dc = dc;
        ifstream in(file);
        string line;
        if(!in.is_open()){
            cout << "Consensus File not found" << "\n";
            return;
        }
        getline(in, line);
        this->consTree = new Tree(line, taxaMap);
        in.close();
        this->rp = new RandomPartition();
        this->taxonCount = taxaMap->size();
    }

    void assignSubTreeToPartition(TreeNode* node, int* rtsp, unordered_map<int,int>& idToIndex){
        if(node->isLeaf()){
            if(idToIndex.find(node->taxon->id) != idToIndex.end()){
                rtsp[idToIndex[node->taxon->id]] = 1;
            }
        }
        else{
            for(auto child : *node->childs){
                assignSubTreeToPartition(child, rtsp,idToIndex);
            }
        }
    }


    double scoreForPartitionByNode(TreeNode* node,CarryInfo* carryInfo, double* nodeDummyTaxaWeights){

        int* rtsP = new int[carryInfo->rts->size()];
        int* dtsp = new int[carryInfo->dts->size()];

        for(int i = 0; i < carryInfo->rts->size(); ++i){
            rtsP[i] = 0;
        }
        for(int i = 0; i < carryInfo->dts->size(); ++i){
            dtsp[i] = 0;
        }

        unordered_map<int, int> idToIndex;
        int i = 0;
        for(auto x : *carryInfo->rts){
            idToIndex[x->id] = i++;
        }
    
        assignSubTreeToPartition(node, rtsP, idToIndex);

        for(i = 0; i < carryInfo->dts->size(); ++i){
            if(nodeDummyTaxaWeights[i] >= .5){
                dtsp[i] = 1;
            }
        }

        if(carryInfo->book == NULL){
            TaxaPerLevelWithPartition* taxas = new TaxaPerLevelWithPartition(carryInfo->rts, carryInfo->dts, rtsP, dtsp, this->taxonCount);
            carryInfo->book = new BookKeepingPerLevel(this->dc, taxas, carryInfo->tid);
        }
        else{
            int rtCount = 0;
            int dtCount = 0;
            vector<int> rtIndices;
            vector<int> dtIndices;

            
            bool changed = false;
            for(i = 0; i < carryInfo->rts->size(); ++i){
                if(rtsP[i] != carryInfo->book->taxaPerLevel->inWhichPartitionRealTaxonByIndex(i)){
                    changed = true;
                    rtCount++;
                    rtIndices.push_back(i);
                }
            }

            for(i = 0; i < carryInfo->dts->size(); ++i){
                if(dtsp[i] != carryInfo->book->taxaPerLevel->inWhichPartitionDummyTaxonByIndex(i)){
                    changed = true;
                    dtCount++;
                    dtIndices.push_back(i);
                }
            }
            if(!changed){
                return carryInfo->score;
            }
            else{
                    if(rtIndices.size() > 0){
                        carryInfo->book->batchTrasferRealTaxon(rtIndices);
                    }

                    for(int j : dtIndices){
                        carryInfo->book->swapTaxon(j, true);
                    }
            }
            
        }
        carryInfo->score = carryInfo->book->calculateScore();

        return carryInfo->score;
    }




    pair<int, double*> dfs(
        TreeNode* node,
        CarryInfo* carryInfo
    ){
        if(node->isLeaf()){
            pair<int, double*> p = make_pair(0, new double[carryInfo->dts->size()]);
            if(carryInfo->isRealTaxon[node->taxon->id]){
                p.first = 1;
            }
            else{
                p.second[carryInfo->inWhichDummyTaxa[node->taxon->id]] = carryInfo->weight[node->taxon->id];
            }
            return p;
        }

        pair<int, double*> p = make_pair(0, new double[carryInfo->dts->size()]);

        for(auto child : *node->childs){

            auto childPair = dfs(child,carryInfo);

            int partASize = childPair.first;
            int partBSize = carryInfo->rts->size() - partASize;


            for(int j = 0; j < carryInfo->dts->size(); ++j){
                p.second[j] += childPair.second[j];
                if(childPair.second[j] >= .5){
                    partASize++;
                }
                else{
                    partBSize++;
                }
            }
            p.first += childPair.first;

            if(partASize >= 1 && partBSize >= 1){
                if(USE_SCORING_IN_CONSENSUS){
                    double score = scoreForPartitionByNode(child, carryInfo, childPair.second);
                    if( carryInfo->minNode == NULL || score > carryInfo->maxScore){
                        carryInfo->maxScore = score;
                        carryInfo->minNode = child;
                        carryInfo->minNodeDummyTaxaWeights = childPair.second;
                    }
                }
                else{
                    double childTotalTaxaCounts = childPair.first;
                    for(int j = 0; j < carryInfo->dts->size(); ++j){
                        childTotalTaxaCounts += childPair.second[j];
                    }
                    double diff = abs(carryInfo->rts->size() + carryInfo->dts->size() - childTotalTaxaCounts);
                    if(carryInfo->minNode == NULL || diff < carryInfo->minDiff){
                        carryInfo->minNode = child;
                        carryInfo->minDiff = diff;
                        carryInfo->minNodeDummyTaxaWeights = childPair.second;
                    }
                    else{
                        delete[] childPair.second;
                    }
                }
            }
            else{
                delete[] childPair.second;
            }
        }
        return p;
        
    }


    virtual MakePartitionReturnType* makePartition(
        vector<RealTaxon*>* rts,
        vector<DummyTaxon*>* dts, 
        int tid
    ){
         
        double* weight = new double[consTree->leavesCount];
        int* inWhichDummyTaxa = new int[consTree->leavesCount];

        bool* isRealTaxon = new bool[consTree->leavesCount];

        for(int i = 0; i < consTree->leavesCount; ++i){
            weight[i] = 0;
            isRealTaxon[i] = false;
            inWhichDummyTaxa[i] = -1;
        }

        for(auto x : *rts){
            weight[x->id] = 1;
            isRealTaxon[x->id] = true;
        }

        int i = 0;

        for(auto x : *dts){
            if(CONSENSUS_WEIGHT_TYPE == ConsensusWeightType::NESTED){
                x->calcDivCoeffs(ScoreNormalizationType::NESTED_NORMALIZATION, weight, 1);
                for(auto y : *x->flattenedRealTaxa){
                    weight[y->id] = 1 / weight[y->id];
                }
            }
            else{
                double sz = x->flattenedTaxonCount;
                for(auto y : *x->flattenedRealTaxa){
                    weight[y->id] += 1. / sz;
                }
            }
            
            for(auto y : *x->flattenedRealTaxa){
                inWhichDummyTaxa[y->id] = i;
            }
            ++i;
        }

        CarryInfo* carryInfo = new CarryInfo();
        carryInfo->weight = weight;
        carryInfo->isRealTaxon = isRealTaxon;
        carryInfo->inWhichDummyTaxa = inWhichDummyTaxa;
        carryInfo->rts = rts;
        carryInfo->dts = dts;
        carryInfo->minNode = NULL;
        carryInfo->maxScore = __DBL_MIN__;
        carryInfo->minDiff = __DBL_MAX__;
        carryInfo->score = 0;
        carryInfo->book = NULL;
        carryInfo->taxas = NULL;
        carryInfo->tid = tid;

        dfs(consTree->root, carryInfo);

        
        if(carryInfo->minNode == NULL){
            cout << "Min Node NULL" << "\n";
            return rp->makePartition(rts, dts, carryInfo->tid);
        }
        int* rtsP = new int[rts->size()];
        int* dtsp = new int[dts->size()];

        for(i = 0; i < rts->size(); ++i){
            rtsP[i] = 0;
        }
        for(i = 0; i < dts->size(); ++i){
            dtsp[i] = 0;
        }

        unordered_map<int, int> idToIndex;
        i = 0;
        for(auto x : *rts){
            idToIndex.insert({x->id, i++});
        }
    
        assignSubTreeToPartition(carryInfo->minNode, rtsP, idToIndex);

        for(i = 0; i < dts->size(); ++i){
            if(carryInfo->minNodeDummyTaxaWeights[i] >= .5){
                dtsp[i] = 1;
            }
        }


        delete carryInfo;

        return new MakePartitionReturnType(rtsP, dtsp);

        // return this->rp->makePartition(realTaxons, dummyTaxons, tid);
    }

};


