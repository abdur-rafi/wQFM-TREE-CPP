#pragma once

#include <vector>
#include <string>
#include <queue>
#include <utility>
#include <unordered_set>

#include "MakePartition.h"
#include "Taxon.h"
#include "Tree.h"
#include "NumSat.h"
#include "Util.h"

using namespace std;

class TaxaPerLevelWithPartition {
public:


    int allRealTaxaCount;
    
    vector<RealTaxon*>* realTaxa;
    vector<DummyTaxon*>* dummyTaxa;

    int* realTaxonPartition;
    int* dummyTaxonPartition;

    int realTaxonCount;
    int dummyTaxonCount;


    bool* isInRealTaxa;
    bool* isInDummyTaxa;
    double* coeffs;
    
    int* inWhichPartition;
    int* inWhichDummyTaxa;
    int* realTaxonIndex;


    int taxonCountsInPartitions[2] = {0, 0};
    int realTaxonCountsInPartitions[2] = {0, 0};
    int dummyTaxonCountsInPartitions[2] = {0, 0};
    int dummyTaxonCountsFlattenedInPartitions[2] = {0, 0};


    bool smallestUnit;

    ~TaxaPerLevelWithPartition(){
        delete[] isInRealTaxa;
        delete[] isInDummyTaxa;
        delete[] coeffs;
        delete[] inWhichPartition;
        delete[] inWhichDummyTaxa;
        delete[] realTaxonIndex;
    }

    TaxaPerLevelWithPartition( vector<RealTaxon*>* rts, vector<DummyTaxon*>* dts, int* rtp, int* dtp, int rtc){

        this->realTaxa = rts;
        this->dummyTaxa = dts;
        this->realTaxonPartition = rtp;
        this->dummyTaxonPartition = dtp;
        this->realTaxonCount = rts->size();
        this->dummyTaxonCount = dts->size();
        this->allRealTaxaCount = rtc;

        if(rts->size() + dts->size() < 4){
            this->smallestUnit = true;
            return;
        }

        this->smallestUnit = false;
        

        this->isInRealTaxa = new bool[this->allRealTaxaCount];
        this->coeffs = new double[this->allRealTaxaCount];
        this->realTaxonIndex = new int[this->allRealTaxaCount];
        this->inWhichDummyTaxa = new int[this->allRealTaxaCount];
        this->isInDummyTaxa = new bool[this->allRealTaxaCount];
        
        for(int i = 0; i < this->allRealTaxaCount; ++i){
            this->isInRealTaxa[i] = false;
            this->coeffs[i] = 0;
            this->realTaxonIndex[i] = -1;
            this->inWhichDummyTaxa[i] = -1;
            this->isInDummyTaxa[i] = false;
        }

        int i = 0;
        

        for(auto x : *realTaxa){
            isInRealTaxa[x->id] = true;
            coeffs[x->id] = 1.0;
            taxonCountsInPartitions[realTaxonPartition[i]]++;
            realTaxonIndex[x->id] = i++;
        }
        this->realTaxonCountsInPartitions[0] = taxonCountsInPartitions[0];
        this->realTaxonCountsInPartitions[1] = taxonCountsInPartitions[1];


        i = 0;
        for(auto x : *dummyTaxa){
            for(auto y : *x->flattenedRealTaxa){
                isInDummyTaxa[y->id] = true;
                inWhichDummyTaxa[y->id] = i;
            }
            x->calcDivCoeffs(NESTED_NORMALIZATION, coeffs, 1.);
            taxonCountsInPartitions[dummyTaxonPartition[i]]++;
            this->dummyTaxonCountsInPartitions[dummyTaxonPartition[i]]++;
            this->dummyTaxonCountsFlattenedInPartitions[dummyTaxonPartition[i]] += x->flattenedTaxonCount;
            
            ++i;
        }
        i = 0;
        this->inWhichPartition = new int[this->allRealTaxaCount];
        // for(int i = 0; i < this->allRealTaxaCount; ++i){
        //     this->inWhichPartition[i] = -1;
        // }
        for(auto x : *realTaxa){
            this->inWhichPartition[x->id] = realTaxonPartition[i++];
        }
        i = 0;
        for(auto x : *dummyTaxa){
            for(auto y : *x->flattenedRealTaxa){
                this->inWhichPartition[y->id] = dummyTaxonPartition[i];
            }
            ++i;
        }
    }


    double getWeight(int realTaxonId){
        return 1.L / this->coeffs[realTaxonId];
    }


    
    int getTaxonCountFlattenedInPartition(int partition){
        return (this->realTaxonCountsInPartitions[partition] + this->dummyTaxonCountsFlattenedInPartitions[partition]);
    }

    int getFlattenedCount(int index){
        return this->dummyTaxa->at(index)->flattenedTaxonCount;
    }

    void batchTransferRealTaxon(vector<int>& indices){
        int netTranser = 0;

        for(int i = 0; i < indices.size(); ++i){
            int index = indices[i];
            int partition = this->realTaxonPartition[index];
            netTranser += (partition == 0 ? 1 : -1);

            this->realTaxonPartition[index] = 1 - partition;
            this->inWhichPartition[this->realTaxa->at(index)->id] = 1 - partition;
        }

        if(netTranser == 0) return;

        int add = 1;
        int sub = 0;
        
        if(netTranser < 0){
            add = 0;
            sub = 1;
            netTranser = -netTranser;
        }

        this->taxonCountsInPartitions[sub] -= netTranser;
        this->taxonCountsInPartitions[add] += netTranser;

        this->realTaxonCountsInPartitions[sub] -= netTranser;
        this->realTaxonCountsInPartitions[add] += netTranser;


    }

    void swapPartitionRealTaxon(int index){
        int currPartition = this->realTaxonPartition[index];
        int switchedPartition = (int) (1 - currPartition);

        this->realTaxonPartition[index] = switchedPartition;
        this->inWhichPartition[this->realTaxa->at(index)->id] = switchedPartition;
        this->taxonCountsInPartitions[currPartition]--;
        this->taxonCountsInPartitions[switchedPartition]++;
        this->realTaxonCountsInPartitions[currPartition]--;
        this->realTaxonCountsInPartitions[switchedPartition]++;
    }
    

    void swapPartitionDummyTaxon(int index){
        int currPartition = this->dummyTaxonPartition[index];
        int switchedPartition = (1 - currPartition);

        this->dummyTaxonPartition[index] = switchedPartition;
        
        for(auto x : *this->dummyTaxa->at(index)->flattenedRealTaxa){
            this->inWhichPartition[x->id] = switchedPartition;
        }
        
        this->taxonCountsInPartitions[currPartition]--;
        this->taxonCountsInPartitions[switchedPartition]++;
        this->dummyTaxonCountsInPartitions[currPartition]--;
        this->dummyTaxonCountsInPartitions[switchedPartition]++;

        this->dummyTaxonCountsFlattenedInPartitions[currPartition] -= this->dummyTaxa->at(index)->flattenedTaxonCount;
        this->dummyTaxonCountsFlattenedInPartitions[switchedPartition] += this->dummyTaxa->at(index)->flattenedTaxonCount;

    }
    int inWhichPartitionDummyTaxonByIndex(int index){
        return this->dummyTaxonPartition[index];
    }

    int inWhichPartitionRealTaxonByIndex(int index){
        return this->realTaxonPartition[index];
    }

    int getRealTaxonCountInPartition(int partition){
        return this->realTaxonCountsInPartitions[partition];
    }
    int getDummyTaxonCountInPartition(int partition){
        return this->dummyTaxonCountsInPartitions[partition];
    }

    int getTaxonCountInPartition(int partition){
        return this->taxonCountsInPartitions[partition];
    }


    Tree* createStar(){
        if(!smallestUnit){
            cout << "Create Star should be called only on smallest unit\n";
            exit(-1);
        }

        // cout << "Creating Star\n";

        Tree* t = new Tree();
        vector<TreeNode*>* childs = new vector<TreeNode*>();
        for(auto x : *this->realTaxa){
            auto y = t->addLeaf(x);
            y->dummyTaxonId = -1;
            childs->push_back(y);
        }

        for(auto x : *this->dummyTaxa){
            auto y = t->addLeaf(NULL);
            y->dummyTaxonId = x->id;
            childs->push_back(y);
        }

        t->root = t->addInternalNode(childs);
        t->root->dummyTaxonId = -1;

        // cout << "Star Created\n";

        return t;
    }

    void printTaxa(){
        cout << "Real Taxa\n";
        for(auto x : *this->realTaxa){
            cout << x->label << " ";
        }
        cout << "\n";
        cout << this->dummyTaxa->size() << " Dummy Taxa\n";
        for(auto x : *this->dummyTaxa){
            cout << "dt id: " << x->id << "\n";
            for(auto y : *x->flattenedRealTaxa){
                cout << y->label << " ";
            }
        }
        cout << "\n";
    
    }
};


class BookKeepingPerTree {
public:
    bool* realTaxaInTree;
    TaxaPerLevelWithPartition* taxaPerLevel;
    double pairsFromPart[2] = {0, 0};
    double realTaxaCountsInPartitions[2] = {0, 0};
    double* dummyTaxonWeightsIndividual;
    double dummyTaxonCountsInPartitions[2] = {0, 0};

    ~BookKeepingPerTree(){
        delete[] dummyTaxonWeightsIndividual;
    }

    BookKeepingPerTree(bool* realTaxaInTree, TaxaPerLevelWithPartition* taxaPerLevel){
        this->realTaxaInTree = realTaxaInTree;
        this->taxaPerLevel = taxaPerLevel;
        double totalTaxon[2] = {0, 0};
        this->dummyTaxonWeightsIndividual = new double[taxaPerLevel->dummyTaxonCount];

        for(int i = 0; i < taxaPerLevel->dummyTaxonCount; ++i){
            this->dummyTaxonWeightsIndividual[i] = 0;
        }
        
        for(int i = 0; i < this->taxaPerLevel->allRealTaxaCount; ++i){
            if(this->realTaxaInTree[i]){
                int partition = this->taxaPerLevel->inWhichPartition[i];
                totalTaxon[partition] += this->taxaPerLevel->getWeight(i);
                if(this->taxaPerLevel->isInDummyTaxa[i]){
                    this->dummyTaxonWeightsIndividual[this->taxaPerLevel->inWhichDummyTaxa[i]] += this->taxaPerLevel->getWeight(i);
                }
                else{
                    this->realTaxaCountsInPartitions[partition]++;
                }
            }
        }

        this->pairsFromPart[0] = totalTaxon[0] * totalTaxon[0];
        this->pairsFromPart[1] = totalTaxon[1] * totalTaxon[1];

        for(int i = 0; i < this->taxaPerLevel->dummyTaxonCount; ++i){
            int partition = this->taxaPerLevel->inWhichPartitionDummyTaxonByIndex(i);
            this->pairsFromPart[partition] -= this->dummyTaxonWeightsIndividual[i] * this->dummyTaxonWeightsIndividual[i];
            this->dummyTaxonCountsInPartitions[partition] += this->dummyTaxonWeightsIndividual[i];
            
        }
        this->pairsFromPart[0] -= this->realTaxaCountsInPartitions[0];
        this->pairsFromPart[1] -= this->realTaxaCountsInPartitions[1];

        this->pairsFromPart[0] /= 2;
        this->pairsFromPart[1] /= 2;        

    }


    double totalQuartets(){
        return this->pairsFromPart[0] * this->pairsFromPart[1];
    }

    double totalQuartetsAfterDummySwap(int dIndex, int toPartition){
        double a = this->pairsFromPart[1 - toPartition] - this->dummyTaxonWeightsIndividual[dIndex] * (
            this->getTotalTaxon(1-toPartition) - 
            this->getDummyTaxonIndiWeight(dIndex)
        );
        double b = this->pairsFromPart[toPartition] + this->dummyTaxonWeightsIndividual[dIndex] * (this->getTotalTaxon(toPartition));
        return a * b;        
    
    }

    double getTotalTaxon(int p){
        return this->realTaxaCountsInPartitions[p] + this->dummyTaxonCountsInPartitions[p];
    }

    double getDummyTaxonIndiWeight(int index){
        return this->dummyTaxonWeightsIndividual[index];
    }

    double totalQuartetsAfterSwap(int rtId, int toPartition){
        if(!this->realTaxaInTree[rtId]){
            return this->totalQuartets();
        }
        return (this->pairsFromPart[1 - toPartition] - (this->getTotalTaxon(1 - toPartition) - 1)) * 
        (this->pairsFromPart[toPartition] + (this->getTotalTaxon(toPartition)));
    }

    void swapRealTaxon(int rtId, int partition){
        if(!this->realTaxaInTree[rtId]) return;

        this->pairsFromPart[ 1 - partition] += this->getTotalTaxon(1 - partition);
        this->pairsFromPart[partition] -= this->getTotalTaxon(partition) - 1;
        
        this->realTaxaCountsInPartitions[partition]--;
        this->realTaxaCountsInPartitions[1 - partition]++;

    }

    void swapDummyTaxon(int index, int partition){

        this->pairsFromPart[ 1 - partition] += this->dummyTaxonWeightsIndividual[index] * this->getTotalTaxon(1 - partition);
        this->pairsFromPart[partition] -= this->dummyTaxonWeightsIndividual[index] * (this->getTotalTaxon(partition) - this->dummyTaxonWeightsIndividual[index]);
        
        this->dummyTaxonCountsInPartitions[partition] -= this->dummyTaxonWeightsIndividual[index];
        this->dummyTaxonCountsInPartitions[1 - partition] += this->dummyTaxonWeightsIndividual[index];

    }

    void batchTranserRealTaxon(vector<int>& rtIds, vector<int>& currPartition){
        int netTranser = 0;
        for(int i = 0; i < rtIds.size(); ++i){
            int rtId = rtIds[i];
            int partition = currPartition.at(i);
            if(this->realTaxaInTree[rtId]){
                netTranser += (partition == 0 ? 1 : -1);
            }
        }

        if(netTranser == 0) return;

        int add = 1;
        int sub = 0;
        
        if(netTranser < 0){
            add = 0;
            sub = 1;
            netTranser = -netTranser;
        }

        this->pairsFromPart[sub] -= (this->getTotalTaxon(sub) - netTranser) * netTranser + netTranser * (netTranser - 1) / 2;
        this->pairsFromPart[add] += (this->getTotalTaxon(add)) * netTranser + netTranser * (netTranser - 1) / 2;
        
        this->realTaxaCountsInPartitions[sub] -= netTranser;
        this->realTaxaCountsInPartitions[add] += netTranser;
        

    }

};

class BookKeepingPerLevel{
public:

    DataContainer* dc;
    TaxaPerLevelWithPartition* taxaPerLevel;
    vector<BookKeepingPerTree*>* bookKeepingPerTreeDCs;
    int tid;

    BookKeepingPerLevel(DataContainer* dc, TaxaPerLevelWithPartition* taxaPerLevelWithPartition, int tid){
        if(DEBUG){
            cerr << "BookKeepingPerLevel constructor\n";
        }
        this->dc = dc;
        this->tid = tid;
        this->taxaPerLevel = taxaPerLevelWithPartition;
        if(this->taxaPerLevel->smallestUnit)
            return;

        this->initialBookKeeping();
        this->bookKeepingPerTreeDCs = new vector<BookKeepingPerTree*>(dc->nTrees);
        for(int i = 0; i < dc->nTrees; ++i){
            this->bookKeepingPerTreeDCs->at(i) = new BookKeepingPerTree(dc->realTaxaInTrees[i], this->taxaPerLevel);
        }
        if(DEBUG){
            cerr << "BookKeepingPerLevel constructor end\n";
        }

    }

    ~BookKeepingPerLevel(){
        if(DEBUG){
            cerr << "BookKeepingPerLevel destructor\n";
        }

        if(this->taxaPerLevel->smallestUnit)
            return;
        for(auto x : *this->bookKeepingPerTreeDCs){
            delete x;
        }
        delete this->bookKeepingPerTreeDCs;

        if(DEBUG){
            cerr << "BookKeepingPerLevel destructor end\n";
        }
    }

    void initialBookKeeping(){

        if(DEBUG){
            cerr << "inside initialBookKeeping\n";
        }
        
        // cout << "inside initialBookKeeping\n";
        // cout << "dummy taxon count " << this->taxaPerLevel->dummyTaxonCount << "\n";
        // for(int i = 0; i < this->taxaPerLevel->realTaxonCount; ++i){
        //     RealTaxon* rt = this->taxaPerLevel->realTaxa->at(i);
        //     cout << "rt id: " << rt->id << " label " << rt->label << "\n";
        // }
        // for(int i = 0; i < this->taxaPerLevel->dummyTaxonCount; ++i){
        //     DummyTaxon* dt = this->taxaPerLevel->dummyTaxa->at(i);
        //     cout << "dt id: " << dt->id << "\n";
        //     for(auto x : *dt->flattenedRealTaxa){
        //         cout << "rt id: " << x->id << " label " << x->label << "\n";
        //     }
        // }

        for(int i = 0; i < this->dc->realTaxaPartitionNodes->size(); ++i){
            PartitionNode* p = this->dc->realTaxaPartitionNodes->at(i);
            // cout << "before reset\n";
            if(DEBUG){
                cerr << "before reset\n";
            }
            p->data[tid]->branch->reset(this->taxaPerLevel->dummyTaxonCount);
            // cout << "after reset\n";
            if(DEBUG){
                cerr << "after reset\n";
            }

            if(this->taxaPerLevel->isInDummyTaxa[i]){
                // cout << "any dummy\n";
                int dtid = this->taxaPerLevel->inWhichDummyTaxa[i];
                // cout << "dtid: " << dtid << "\n";
                int partition = this->taxaPerLevel->inWhichPartitionDummyTaxonByIndex(dtid);
                // cout << "partition: " << partition << "\n";
                p->data[tid]->branch->dummyTaxaWeightsIndividual[dtid] = this->taxaPerLevel->getWeight(i);
                p->data[tid]->branch->totalTaxaCounts[partition] += this->taxaPerLevel->getWeight(i);

            }
            else{
                int partition = this->taxaPerLevel->inWhichPartition[i];
                // cout << "partition: " << partition << "\n";

                


                p->data[tid]->branch->realTaxaCounts[partition] = 1;
                p->data[tid]->branch->totalTaxaCounts[partition] = 1;
            }

        }

        // cout << "before topsorted initialBookKeeping\n";


        int sz = this->dc->topSortedPartitionNodes->size();
        for(int i = sz - 1; i >  -1; --i){
            PartitionNode* p = this->dc->topSortedPartitionNodes->at(i);
            if(p->isLeaf){
                continue;
            }
            else{
                p->data[tid]->branch->reset(this->taxaPerLevel->dummyTaxonCount);
                for(PartitionNode* child : *p->children){
                    p->data[tid]->branch->addToSelf(child->data[tid]->branch);
                }
            }
        }

        // cout << "after topsorted initialBookKeeping\n";


        for(PartitionByTreeNode* p : *this->dc->partitionsByTreeNodes){
            p->scoreCalculator[tid]->reset(this->taxaPerLevel->dummyTaxonPartition);
        }

        // cout << "exit initialBookKeeping\n";

        // SubProblemsQueue.instance.initScoreCalculators(this->tid, this->taxaPerLevel->dummyTaxonPartition);

        if(DEBUG){
            cerr << "exit initialBookKeeping\n";
        }

    }


    double calculateScore(){
        if(DEBUG){
            cerr << "inside calculateScore\n";
        }
        double score = 0;
        double totalQuartets = 0;

        for(PartitionByTreeNode* p : *this->dc->partitionsByTreeNodes){
            score += p->scoreCalculator[tid]->score() * p->count;
        }
        // score = SubProblemsQueue.instance.calcScores(tid);

        for(BookKeepingPerTree* bt : *this->bookKeepingPerTreeDCs){
            totalQuartets += bt->totalQuartets();
        }

        if(DEBUG){
            cerr << "exit calculateScore\n";
        }
        

        return scoreEqn(totalQuartets, score);
    }

    double calculateScoreAndGains(double** realTaxaGains, double* dummyTaxaGains){

        if(DEBUG){
            cerr << "inside calculateScoreAndGains\n";
        }

        double totalSat = 0;
        
        for(PartitionNode* p : *this->dc->topSortedPartitionNodes){
            p->data[tid]->gainsForSubTree[0] = 0;
            p->data[tid]->gainsForSubTree[1] = 0;

        }

        if(DEBUG){
            // cerr << "before calculateScoreAndGains\n";
            cerr << this->dc->partitionsByTreeNodes->size() << "\n";
            cerr << "after gain for sub tree init\n";
        }

        for(PartitionByTreeNode* p : *this->dc->partitionsByTreeNodes){
            double score = p->scoreCalculator[tid]->score();
            double** branchGainsForRealTaxa = p->scoreCalculator[tid]->gainRealTaxa(score, p->count);
            
            p->scoreCalculator[tid]->gainDummyTaxa(score, p->count, dummyTaxaGains);
            score *= p->count;

            totalSat += score;

            for(int i = 0; i < p->partitionNodes->size(); ++i){
                if(DEBUG){
                    cerr << "before addArrayToFirst\n";
                }
                addArrayToFirst(p->partitionNodes->at(i)->data[tid]->gainsForSubTree, branchGainsForRealTaxa[i], 2);
                if(DEBUG){
                    cerr << "after addArrayToFirst\n";
                }
            }
            
            // for(int i = 0; i < p->partitionNodes->size(); ++i){
            //     delete[] branchGainsForRealTaxa[i];
            // }
            // delete[] branchGainsForRealTaxa;
        }

        if(DEBUG){
            cerr << "after score and gain loop\n";
        }

        // cout << "totalSat: " << totalSat << "\n";

        // var x = SubProblemsQueue.instance.calcGains(tid, this->taxaPerLevel->dummyTaxonCount);
        // totalSat = x.first;
        // Utility.addArrayToFirst(dummyTaxaGains, x.second);

        for(PartitionNode* p : *this->dc->topSortedPartitionNodes){
            // if(p->isLeaf){
            //     continue;
            // }
            for(PartitionNode* childs : *p->children){
                addArrayToFirst(childs->data[tid]->gainsForSubTree, p->data[tid]->gainsForSubTree, 2);
            }
        }

        if(DEBUG){
            cerr << "after adding gains for sub tree\n";
        }

        double currTotalQuartets = 0;
        double* dtTotals = new double[this->taxaPerLevel->dummyTaxonCount];
        
        for(BookKeepingPerTree* bkpt : *this->bookKeepingPerTreeDCs){
            currTotalQuartets += bkpt->totalQuartets();
            for(int i = 0;i < this->taxaPerLevel->dummyTaxonCount; ++i){
                dtTotals[i] += bkpt->totalQuartetsAfterDummySwap(i, 1 - this->taxaPerLevel->inWhichPartitionDummyTaxonByIndex(i));
            }
        }

        double totalScore = scoreEqn(currTotalQuartets, totalSat);

        for(int i = 0; i < this->taxaPerLevel->realTaxonCount; ++i){
            RealTaxon* rt = this->taxaPerLevel->realTaxa->at(i);
            int partition = this->taxaPerLevel->inWhichPartitionRealTaxonByIndex(i);
            addArrayToFirst(realTaxaGains[i], this->dc->realTaxaPartitionNodes->at(rt->id)->data[tid]->gainsForSubTree, 2);
            double totalQuartetsAfterTransferringi = 0;
            for(BookKeepingPerTree* bkpt : *this->bookKeepingPerTreeDCs){
                totalQuartetsAfterTransferringi += bkpt->totalQuartetsAfterSwap(i, 1 - partition);
            }
            realTaxaGains[i][partition] += totalSat;
            realTaxaGains[i][partition] = scoreEqn(totalQuartetsAfterTransferringi, realTaxaGains[i][partition]);
            realTaxaGains[i][partition] -= totalScore;   
        }

        for(int i = 0; i < this->taxaPerLevel->dummyTaxonCount; ++i){
            
            dummyTaxaGains[i] = scoreEqn(
                dtTotals[i],
                dummyTaxaGains[i] + totalSat
            ) - totalScore;


        }

        delete[] dtTotals;

        if(DEBUG){
            cerr << "exit calculateScoreAndGains\n";
        }

        return totalScore;
    }


    void batchTrasferRealTaxon(vector<int>& realTaxonIndices){

        if(DEBUG){
            cerr << "inside batchTrasferRealTaxon\n";
        }

        vector<int> currPartitions;
        vector<int> realTaxonIds;

        for(int i = 0; i < realTaxonIndices.size(); ++i){
            int index = realTaxonIndices.at(i);
            int partition = this->taxaPerLevel->inWhichPartitionRealTaxonByIndex(index);
            currPartitions.push_back(partition);
            realTaxonIds.push_back(this->taxaPerLevel->realTaxa->at(index)->id);
        }
        

        for(BookKeepingPerTree* bkpt : *this->bookKeepingPerTreeDCs){
            bkpt->batchTranserRealTaxon(realTaxonIds, currPartitions);
        }

        {   

            queue<pair<PartitionNode*, int>> q;


            for(int rtId : realTaxonIds){
                q.push(make_pair(this->dc->realTaxaPartitionNodes->at(rtId), this->taxaPerLevel->inWhichPartition[rtId]));
            }

            unordered_set<PartitionNode*> pst;

            while(!q.empty()){
                auto f = q.front();
                q.pop();
                for(auto x : *f.first->parents){
                    q.push(make_pair(x, f.second));
                }
                f.first->data[tid]->branch->cumulateTransfer(f.second);
                pst.insert(f.first);
            }


            for(auto x : pst){
                for(auto y : *x->nodePartitions){
                    y->partitionByTreeNode->batchTransfer(y->index, x->data[tid]->branch->netTranser, tid);
                }
                x->data[tid]->branch->batchTransferRealTaxon();
            }

        }

        this->taxaPerLevel->batchTransferRealTaxon(realTaxonIndices);

        if(DEBUG){
            cerr << "exit batchTrasferRealTaxon\n";
        }
    }


    void swapRealTaxon(int index){

        if(DEBUG){
            cerr << "inside swapRealTaxon\n";
        }
        
        int partition = this->taxaPerLevel->inWhichPartitionRealTaxonByIndex(index);
        this->taxaPerLevel->swapPartitionRealTaxon(index);
        int rtId = this->taxaPerLevel->realTaxa->at(index)->id;

        for(BookKeepingPerTree* bkpt : *this->bookKeepingPerTreeDCs){
            bkpt->swapRealTaxon(rtId, partition);
        }
        
        queue<PartitionNode*> q;
        q.push(this->dc->realTaxaPartitionNodes->at(rtId));

        while(!q.empty()){
            PartitionNode* f = q.front();
            q.pop();

            for(PartitionByTreeNodeWithIndex* p : *f->nodePartitions){
                p->partitionByTreeNode->scoreCalculator[tid]->swapRealTaxon(
                    p->index,
                    partition
                );

            }
            for(PartitionNode* x : *f->parents){
                q.push(x);
            }
            f->data[tid]->branch->swapRealTaxa(partition);
        }

        if(DEBUG){
            cerr << "exit swapRealTaxon\n";
        }
    }

    void swapDummyTaxon(int index){

        if(DEBUG){
            cerr << "inside swapDummyTaxon\n";
        }

        int partition = this->taxaPerLevel->inWhichPartitionDummyTaxonByIndex(index);
        this->taxaPerLevel->swapPartitionDummyTaxon(index);
        
        for(BookKeepingPerTree* bkpt : *this->bookKeepingPerTreeDCs){
            bkpt->swapDummyTaxon(index, partition);
        }

        for(PartitionByTreeNode* p : *this->dc->partitionsByTreeNodes){
            p->scoreCalculator[tid]->swapDummyTaxon(index, partition);
        }

        // SubProblemsQueue.instance.swapDT(this->tid, index, partition);

        unordered_set<PartitionNode*> st;
        queue<PartitionNode*> q;

        DummyTaxon* dt = this->taxaPerLevel->dummyTaxa->at(index);

        for(RealTaxon* rt : *dt->flattenedRealTaxa){
            this->dc->realTaxaPartitionNodes->at(rt->id)->data[tid]->branch->swapDummyTaxon(index, partition);
            for(PartitionNode* p : *this->dc->realTaxaPartitionNodes->at(rt->id)->parents){
                if(st.find(p) == st.end()){
                    st.insert(p);
                    q.push(p);
                }
            }
        }

        while(!q.empty()){
            PartitionNode* f = q.front();
            q.pop();
            f->data[tid]->branch->swapDummyTaxon(index, partition);
            for(PartitionNode* p : *f->parents){
                if(st.find(p) == st.end()){
                    st.insert(p);
                    q.push(p);
                }
            }
        }

        if(DEBUG){
            cerr << "exit swapDummyTaxon\n";
        }

    }

    void swapTaxon(int index, bool isDummy){
        if(isDummy) this->swapDummyTaxon(index);
        else this->swapRealTaxon(index);
    }

    TaxaPerLevelWithPartition** divide(IMakePartition* makePartition){
        if(DEBUG){
            cerr << "inside divide\n";
        }

        vector<RealTaxon*>** rts = new vector<RealTaxon*>*[2];
        vector<DummyTaxon*>** dts = new vector<DummyTaxon*>*[2];
        for(int i = 0; i < 2; ++i){
            // cout << this->taxaPerLevel->getRealTaxonCountInPartition(i) << "\n";
            // cout << this->taxaPerLevel->getDummyTaxonCountInPartition(i) << "\n";
            rts[i] = new vector<RealTaxon*>(this->taxaPerLevel->getRealTaxonCountInPartition(i), NULL);
            dts[i] = new vector<DummyTaxon*>(this->taxaPerLevel->getDummyTaxonCountInPartition(i), NULL);
        }

        int index[2];
        index[0] = 0;
        index[1] = 0;


        for(auto x : *this->taxaPerLevel->realTaxa){
            int part = this->taxaPerLevel->inWhichPartition[x->id];
            rts[part]->at(index[part]) = x;
            index[part]++;
        }
        index[0] = 0;
        index[1] = 0;
        int i = 0;
        for(auto x : *this->taxaPerLevel->dummyTaxa){
            int part = this->taxaPerLevel->inWhichPartitionDummyTaxonByIndex(i++);
            // cout << "part: " << part << "\n";
            dts[part]->at(index[part]) = x;
            // cout << "dts size " << dts[part]->size() << "\n";
            index[part]++;
        }

        // cout << "here\n";


        // ith dummy taxon for ith partition
        DummyTaxon** newDt = new DummyTaxon*[2];
        
        

        // BookKeepingPerLevel[] bookKeepingPerLevels = new BookKeepingPerLevel[2];
        TaxaPerLevelWithPartition** taxaPerLevelWithPartitions = new TaxaPerLevelWithPartition*[2];

        for( i = 0; i < 2; ++i){
            newDt[i] = new DummyTaxon(rts[1 - i], dts[1 - i]);
            
            // cout << "dts size " << dts[i]->size() << "\n";
            vector<DummyTaxon*>* dtsWithNewDt = new vector<DummyTaxon*>(dts[i]->size() + 1, NULL);

            for(int j = 0; j < dts[i]->size(); ++j){
                dtsWithNewDt->at(j) = dts[i]->at(j);
            }
            dtsWithNewDt->at(dtsWithNewDt->size() - 1) = newDt[i];

            if(rts[i]->size() + dtsWithNewDt->size() > 3){

                auto y = makePartition->makePartition(rts[i], dtsWithNewDt, this->tid);
                taxaPerLevelWithPartitions[i] = new TaxaPerLevelWithPartition(
                    rts[i], dtsWithNewDt, 
                    y->realTaxonPartition, 
                    y->dummyTaxonPartition, 
                    this->dc->taxa->size()
                );
            }
            else{
                taxaPerLevelWithPartitions[i] = new TaxaPerLevelWithPartition(
                    rts[i], dtsWithNewDt, 
                    NULL, NULL,
                    this->dc->taxa->size()
                );
            }
            
            // bookKeepingPerLevels[i] = new BookKeepingPerLevel(this->geneTrees,x);
        }



        delete[] rts;
        delete[] dts;
        delete[] newDt;

        // cout << "taxa divide \n";

        // taxaPerLevelWithPartitions[0]->printTaxa();
        // taxaPerLevelWithPartitions[1]->printTaxa();

        if(DEBUG){
            cerr << "exit divide\n";
        }


        return taxaPerLevelWithPartitions;  
    }
    
    

};

