#pragma once

#include <vector>
#include <DataContainer.h>


class NumSatCalculatorNode {
public:
    virtual double score() = 0;

    virtual double** gainRealTaxa(double originalScore, double multiplier) = 0;

    virtual void swapRealTaxon(int branchIndex, int currPartition) = 0;

    virtual void swapDummyTaxon(int dummyIndex, int currPartition) = 0;

    virtual void gainDummyTaxa(double originalScore, double multiplier, double* dummyTaxaGains) = 0;

    virtual void batchTransferRealTaxon(int branchIndex, int netTranser) = 0;    
    
};


class NumSatCalculatorBinaryNode: public NumSatCalculatorNode {
public:
    vector<Branch*>* branches;

    double subs[3][2];
    int nDummyTaxa;
    int* dummyTaxaPartition;
    double** gainsOfBranches;



    NumSatCalculatorBinaryNode(vector<Branch*>* b, int* dummyTaxaToPartitionMap) {
        this->dummyTaxaPartition = dummyTaxaToPartitionMap;
        this->branches = b;
        this->nDummyTaxa = b->at(0)->dummyTaxaCount;

        for(int i = 0; i < 3; ++i){
            subs[i][0] = 0;
            subs[i][1] = 0;
            for(int j = 0; j < this->nDummyTaxa; ++j){
                int pIndex = this->dummyTaxaPartition[j];
                if(pIndex == 0)
                    subs[i][0] += b->at(i)->dummyTaxaWeightsIndividual[j] * b->at((i+1) % 3)->dummyTaxaWeightsIndividual[j];
                else if(pIndex == 1){
                    subs[i][1] += (b->at(i)->dummyTaxaWeightsIndividual[j] * (b->at(i)->dummyTaxaWeightsIndividual[j]) ); 
                }
                else{
                    // System.out.println("error");
                }
            }
            subs[i][1] += b->at(i)->realTaxaCounts[1];
        }

        gainsOfBranches = new double*[3];
        for(int i = 0; i < 3; ++i){
            gainsOfBranches[i] = new double[2];
        }

    }

    double scoreOf2Branch(int i) {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;

        double csubs[2];
        csubs[0] = subs[i][0];
        csubs[1] = subs[k][1];

        double score = satisfiedEqn(
            branches->at(i)->totalTaxaCounts[0], 
            branches->at(j)->totalTaxaCounts[0],
            branches->at(k)->totalTaxaCounts[1], 
            csubs
        );


        return score;

    }

    double satisfiedEqn(double a1, double a2, double b3, double* subs) {
        return (a1 * a2 - subs[0]) * ((( b3 * b3 ) - subs[1] ) / 2);
    }

    double scoreAfterRTSwap(int branchIndex, int currPartition){
        double res = 0;

        double adjust[3][2];

        for(int i = 0; i < 3; ++i){
            adjust[i][0] = 0;
            adjust[i][1] = 0;
        }

        adjust[branchIndex][currPartition] = -1;
        adjust[branchIndex][1 - currPartition] = 1;

        for(int i = 0; i < 3; ++i){
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;

            double csubs[2];            
            csubs[0] = subs[i][0];
            csubs[1] = subs[k][1];

            res += satisfiedEqn(
                branches->at(i)->totalTaxaCounts[0] + adjust[i][0], 
                branches->at(j)->totalTaxaCounts[0] + adjust[j][0],
                branches->at(k)->totalTaxaCounts[1] + adjust[k][1], 
                csubs
            );

        }

        return res;
    }

    


    double scoreAfterDTSwap(int dummyIndex, int currPartition){
        double res = 0;

        // double[][] adjust = new double[3][2];
        double adjust[3][2];

        for(int i = 0; i < 3; ++i){
            adjust[i][currPartition] = -branches->at(i)->dummyTaxaWeightsIndividual[dummyIndex];
            adjust[i][1 - currPartition] = branches->at(i)->dummyTaxaWeightsIndividual[dummyIndex];
        }

        for(int i = 0; i < 3; ++i){
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;

            double csubs[2];
            csubs[0] = subs[i][0];
            csubs[1] = subs[k][1];

            res += satisfiedEqn(
                branches->at(i)->totalTaxaCounts[0] + adjust[i][0], 
                branches->at(j)->totalTaxaCounts[0] + adjust[j][0],
                branches->at(k)->totalTaxaCounts[1] + adjust[k][1], 
                csubs
            );

        }

        return res;
    }


    


    virtual double score() {
        double res = 0;

        for (int i = 0; i < 3; ++i) {
            res += scoreOf2Branch(i);
        }

        return res;
    }


    
    virtual double** gainRealTaxa(double originalScore, double multiplier){
        for(int i = 0; i < 3; ++i){
            gainOf1BranchRealTaxa(i, originalScore, multiplier);
        }
        return this->gainsOfBranches;
    }

    

    virtual void swapRealTaxon(int branchIndex, int currPartition){
        
        if(currPartition == 0){
            this->subs[branchIndex][1] += 1;
        }
        else{
            this->subs[branchIndex][1] -= 1;
        }

    }
    

    virtual void batchTransferRealTaxon(int branchIndex, int netTranser){
        // negative if transfering from 1 to 0
        // positive if transfering from 0 to 1
        if(netTranser < 0){
            this->subs[branchIndex][1] += netTranser;
        }
        else{
            this->subs[branchIndex][1] += netTranser;
        }

        // this->branches->at(branchIndex)->batchTransferRealTaxon(netTranser);
    }

    virtual void swapDummyTaxon(int dummyIndex, int currPartition){
            
        int switchedPartition = 1 - currPartition;
        double currDummyCountCurrBranch, currDummyCountNextBranch;

        for(int i = 0; i < 3; ++i){

            currDummyCountCurrBranch = branches->at(i)->dummyTaxaWeightsIndividual[dummyIndex];
            currDummyCountNextBranch = branches->at((i + 1) % 3)->dummyTaxaWeightsIndividual[dummyIndex];

            if(switchedPartition == 1){                 
                this->subs[i][0] -= currDummyCountCurrBranch * currDummyCountNextBranch;
                // 
                this->subs[i][1] += currDummyCountCurrBranch * currDummyCountCurrBranch;
            }
            else{
                this->subs[i][0] += currDummyCountCurrBranch * currDummyCountNextBranch;
                this->subs[i][1] -= currDummyCountCurrBranch * currDummyCountCurrBranch;
            }

        }

            
    }


    void gainOf1BranchRealTaxa(int i, double originalScore, double multiplier){

        Branch* curr = branches->at(i);
        for(int p = 0; p < 2; ++p){
            if(curr->realTaxaCounts[p] > 0){
                this->swapRealTaxon(i, p);
                // curr.swapRealTaxa(p);
                gainsOfBranches[i][p] = multiplier * (scoreAfterRTSwap(i, p) - originalScore);
                this->swapRealTaxon(i, 1 - p);
                // curr.swapRealTaxa(1 - p);
            }
        }
    }
    
    virtual void gainDummyTaxa(double originalScore, double multiplier, double* dummyTaxaGains){

        for(int i = 0; i < this->nDummyTaxa; ++i){
            int currPartition = this->dummyTaxaPartition[i];
            int switchedPartition = 1 - currPartition;

            this->swapDummyTaxon(i, currPartition);
            // for(Branch b : this->branches){
            //     b.swapDummyTaxon(i, currPartition);
            // }

            double newScore = scoreAfterDTSwap(i, currPartition);

            dummyTaxaGains[i] +=  multiplier * (newScore - originalScore);
            
            this->swapDummyTaxon(i, switchedPartition);
            // for(Branch b : this->branches){
            //     b.swapDummyTaxon(i, switchedPartition);
            // }
        }
    }


    ~NumSatCalculatorBinaryNode(){
        for(int i = 0; i < 3; ++i){
            delete[] gainsOfBranches[i];
        }
        delete[] gainsOfBranches;
    }


};


class NumSatCalculatorNodeE : public NumSatCalculatorNode{

public:
    vector<Branch*>* branches;
    
    double* pairsBFromSingleBranch;
    double** sumPairsBranch;
    double* dummyTaxaWeightsIndividual;
    // double totalTaxaA;
    // double totalTaxaB;
    double totalTaxa[2];
    double*** pairs;
    // double[][][] subs;
    double sumPairs[2];

    double nonQuartets;
    // double sumPairsB;
    double sumPairsBSingleBranch;
    

    int nDummyTaxa;



    double calcNonQuartets(){
        double q = 0;
        for(int i = 0; i < branches->size(); ++i){
            for(int j = i + 1; j < branches->size(); ++j){
                q += pairs[i][j][0] * (sumPairs[1] - sumPairsBranch[i][1] - sumPairsBranch[j][1] + pairs[i][j][1]);
            }
        }
        return q;
    }

    double changeAmount(int branchIndex){
        double q = 0;
        for(int i = 0; i < branches->size(); ++i){
            if(i == branchIndex) continue;
            int mni = branchIndex > i ? i : branchIndex;
            int mxi = branchIndex > i ? branchIndex : i;
            q += pairs[mni][mxi][0] * (sumPairs[1] - sumPairsBranch[mni][1] - sumPairsBranch[mxi][1] + pairs[mni][mxi][1]);
            q += pairs[mni][mxi][1] * (sumPairs[0] - sumPairsBranch[mni][0] - sumPairsBranch[mxi][0] + pairs[mni][mxi][0]);
            
        }
        return q;
    }

    int* dummyTaxaPartition;

    NumSatCalculatorNodeE(vector<Branch*>* b, int* dummyTaxaToPartitionMap) {

        // System.out.println("kjasdfj");

        this->dummyTaxaPartition = dummyTaxaToPartitionMap;
        this->branches = b;
        this->dummyTaxaWeightsIndividual = new double[b->at(0)->dummyTaxaCount];
        this->totalTaxa[0] = 0;
        this->totalTaxa[1] = 0;

        this->nDummyTaxa = b->at(0)->dummyTaxaCount;

        pairsBFromSingleBranch = new double[b->size()];
        // sumPairsBranch = new double[b->size()][2];
        this->sumPairsBranch = new double*[b->size()];
        for(int i = 0; i < b->size(); ++i){
            this->sumPairsBranch[i] = new double[2];
            this->sumPairsBranch[i][0] = 0;
            this->sumPairsBranch[i][1] = 0;
        }
        this->sumPairs[0] = 0;
        this->sumPairs[1] = 0;

        this->sumPairsBSingleBranch = 0;
        // this->pairs = new double[b->size()][b->size()][2];
        this->pairs = new double**[b->size()];
        for(int i = 0; i < b->size(); ++i){
            this->pairs[i] = new double*[b->size()];
            for(int j = 0; j < b->size(); ++j){
                this->pairs[i][j] = new double[2];
                this->pairs[i][j][0] = 0;
                this->pairs[i][j][1] = 0;
            }
        }
        
        this->nonQuartets = 0;

        for(int i = 0; i < b->size(); ++i){
            for(int j = i + 1; j < b->size(); ++j){
                this->pairs[i][j][0] = b->at(i)->totalTaxaCounts[0] * b->at(j)->totalTaxaCounts[0];
                this->pairs[i][j][1] = b->at(i)->totalTaxaCounts[1] * b->at(j)->totalTaxaCounts[1];

                for(int k = 0; k < this->nDummyTaxa; ++k){
                    int partition = this->dummyTaxaPartition[k];
                    // subs[i][j][partition] += b->at(i)->dummyTaxaWeightsIndividual[k] * b->at(j).dummyTaxaWeightsIndividual[k];
                    this->pairs[i][j][partition] -= b->at(i)->dummyTaxaWeightsIndividual[k] * b->at(j)->dummyTaxaWeightsIndividual[k];
                }
                // this->pairs[i][j][0] -= subs[i][j][0];
                // this->pairs[i][j][1] -= subs[i][j][1];

                sumPairsBranch[i][0] += this->pairs[i][j][0];
                sumPairsBranch[j][0] += this->pairs[i][j][0];
                sumPairsBranch[i][1] += this->pairs[i][j][1];
                sumPairsBranch[j][1] += this->pairs[i][j][1];

                this->sumPairs[0] += this->pairs[i][j][0];
                this->sumPairs[1] += this->pairs[i][j][1];
                
            }

            pairsBFromSingleBranch[i] = b->at(i)->totalTaxaCounts[1] * b->at(i)->totalTaxaCounts[1];
            for(int k = 0; k < this->nDummyTaxa; ++k){
                int partition = this->dummyTaxaPartition[k];
                this->dummyTaxaWeightsIndividual[k] += b->at(i)->dummyTaxaWeightsIndividual[k];
                if(partition == 1){
                    pairsBFromSingleBranch[i] -= b->at(i)->dummyTaxaWeightsIndividual[k] * b->at(i)->dummyTaxaWeightsIndividual[k];
                }
            }
            pairsBFromSingleBranch[i] -= b->at(i)->realTaxaCounts[1];
            pairsBFromSingleBranch[i] /= 2;
            this->sumPairsBSingleBranch += pairsBFromSingleBranch[i];

            this->totalTaxa[0] += b->at(i)->totalTaxaCounts[0];
            this->totalTaxa[1] += b->at(i)->totalTaxaCounts[1];
            
        }

        this->nonQuartets = this->calcNonQuartets();

    }

    virtual double score(){
        double res = 0;
        for(int i = 0; i < this->branches->size(); ++i){
            res -=  pairsBFromSingleBranch[i] * sumPairsBranch[i][0];
        }
        res += this->sumPairs[0] * this->sumPairsBSingleBranch;
        res += (this->nonQuartets / 2);
        return res;
    }

    virtual void swapRealTaxon(int branchIndex, int currPartition){
        
        this->nonQuartets -= changeAmount(branchIndex);
        for(int i = 0; i < this->branches->size(); ++i){
            if(branchIndex == i){
                this->sumPairsBranch[i][1 - currPartition] += (this->totalTaxa[1-currPartition] - this->branches->at(i)->totalTaxaCounts[1-currPartition]);
                this->sumPairsBranch[i][currPartition] -= (this->totalTaxa[currPartition] - this->branches->at(i)->totalTaxaCounts[currPartition]);

                this->sumPairs[1 - currPartition] += (this->totalTaxa[1-currPartition] - this->branches->at(i)->totalTaxaCounts[1-currPartition]);
                this->sumPairs[currPartition] -= (this->totalTaxa[currPartition] - this->branches->at(i)->totalTaxaCounts[currPartition]);

            }
            else{
                int mni = branchIndex > i ? i : branchIndex;
                int mxi = branchIndex > i ? branchIndex : i;


                this->pairs[mni][mxi][currPartition] -= this->branches->at(i)->totalTaxaCounts[currPartition];
                this->pairs[mni][mxi][1 - currPartition] += this->branches->at(i)->totalTaxaCounts[1 - currPartition];


                this->sumPairsBranch[i][currPartition] -= this->branches->at(i)->totalTaxaCounts[currPartition];
                this->sumPairsBranch[i][1 - currPartition] += this->branches->at(i)->totalTaxaCounts[1 - currPartition];
            }
            

        }

        if(currPartition == 1){
            pairsBFromSingleBranch[branchIndex] -= this->branches->at(branchIndex)->totalTaxaCounts[1] - 1;
            this->sumPairsBSingleBranch -= this->branches->at(branchIndex)->totalTaxaCounts[1] - 1;                
        }
        else{
            pairsBFromSingleBranch[branchIndex] += this->branches->at(branchIndex)->totalTaxaCounts[1];
            this->sumPairsBSingleBranch += this->branches->at(branchIndex)->totalTaxaCounts[1];                
        }
        this->totalTaxa[currPartition] -= 1;
        this->totalTaxa[1 - currPartition] += 1;
        this->nonQuartets += changeAmount(branchIndex);


    }


    void undoSwapRealTaxon(int branchIndex, int currPartition){

        this->nonQuartets -= changeAmount(branchIndex);
        for(int i = 0; i < this->branches->size(); ++i){
            if(branchIndex == i){
                this->sumPairsBranch[i][1 - currPartition] += (this->totalTaxa[1-currPartition] - (this->branches->at(i)->totalTaxaCounts[1-currPartition] - 1 ) );
                this->sumPairsBranch[i][currPartition] -= (this->totalTaxa[currPartition] - (this->branches->at(i)->totalTaxaCounts[currPartition] + 1));

                this->sumPairs[1 - currPartition] += (this->totalTaxa[1-currPartition] - (this->branches->at(i)->totalTaxaCounts[1-currPartition] - 1));
                this->sumPairs[currPartition] -= (this->totalTaxa[currPartition] - (this->branches->at(i)->totalTaxaCounts[currPartition] + 1));

            }
            else{
                int mni = branchIndex > i ? i : branchIndex;
                int mxi = branchIndex > i ? branchIndex : i;


                this->pairs[mni][mxi][currPartition] -= this->branches->at(i)->totalTaxaCounts[currPartition];
                this->pairs[mni][mxi][1 - currPartition] += this->branches->at(i)->totalTaxaCounts[1 - currPartition];


                this->sumPairsBranch[i][currPartition] -= this->branches->at(i)->totalTaxaCounts[currPartition];
                this->sumPairsBranch[i][1 - currPartition] += this->branches->at(i)->totalTaxaCounts[1 - currPartition];
            }
            

        }

        if(currPartition == 1){
            pairsBFromSingleBranch[branchIndex] -= (this->branches->at(branchIndex)->totalTaxaCounts[1] + 1) - 1;
            this->sumPairsBSingleBranch -= (this->branches->at(branchIndex)->totalTaxaCounts[1] + 1) - 1;                
        }
        else{
            pairsBFromSingleBranch[branchIndex] += (this->branches->at(branchIndex)->totalTaxaCounts[1] - 1);
            this->sumPairsBSingleBranch += (this->branches->at(branchIndex)->totalTaxaCounts[1] - 1);                
        }
        this->totalTaxa[currPartition] -= 1;
        this->totalTaxa[1 - currPartition] += 1;
        this->nonQuartets += changeAmount(branchIndex);


    }


    virtual void swapDummyTaxon(int dummyIndex, int currPartition){
        

        for(int i = 0; i < this->branches->size(); ++i){
            double wi = this->branches->at(i)->dummyTaxaWeightsIndividual[dummyIndex];

            for(int j = i + 1; j < this->branches->size(); ++j){
                
                double wj = this->branches->at(j)->dummyTaxaWeightsIndividual[dummyIndex];
                
                double inc = wi * this->branches->at(j)->totalTaxaCounts[1 - currPartition] + wj * this->branches->at(i)->totalTaxaCounts[1-currPartition];
                this->pairs[i][j][1 - currPartition] += inc;
                
                double dec = wi * (this->branches->at(j)->totalTaxaCounts[currPartition] - wj) + wj * (this->branches->at(i)->totalTaxaCounts[currPartition] - wi);
                this->pairs[i][j][currPartition] -= dec;
                
                this->sumPairsBranch[i][1 - currPartition] += inc;
                this->sumPairsBranch[j][1 - currPartition] += inc;
                this->sumPairsBranch[i][currPartition] -= dec;
                this->sumPairsBranch[j][currPartition] -= dec;

                this->sumPairs[1 - currPartition] += inc;
                this->sumPairs[currPartition] -= dec;

                // this->subs[i][j][currPartition] -= wi * wj;
                // this->subs[i][j][1 - currPartition] += wi * wj;
            }

            if(currPartition == 1){
                this->pairsBFromSingleBranch[i] -= (this->branches->at(i)->totalTaxaCounts[1] - wi) * wi;
                this->sumPairsBSingleBranch -= (this->branches->at(i)->totalTaxaCounts[1] - wi) * wi; 

            }
            else{

                this->pairsBFromSingleBranch[i] += (this->branches->at(i)->totalTaxaCounts[1]) * wi;
                this->sumPairsBSingleBranch += (this->branches->at(i)->totalTaxaCounts[1]) * wi;                     
            }
        }

        this->totalTaxa[1 - currPartition] += this->dummyTaxaWeightsIndividual[dummyIndex];
        this->totalTaxa[currPartition] -= this->dummyTaxaWeightsIndividual[dummyIndex];

        this->nonQuartets = calcNonQuartets();
        
    }

    void undoSwapDummyTaxon(int dummyIndex, int currPartition){
        

        for(int i = 0; i < this->branches->size(); ++i){
            double wi = this->branches->at(i)->dummyTaxaWeightsIndividual[dummyIndex];

            for(int j = i + 1; j < this->branches->size(); ++j){
                
                double wj = this->branches->at(j)->dummyTaxaWeightsIndividual[dummyIndex];
                
                double inc = wi * (this->branches->at(j)->totalTaxaCounts[1 - currPartition] - wj) + wj * (this->branches->at(i)->totalTaxaCounts[1-currPartition] - wi);
                this->pairs[i][j][1 - currPartition] += inc;
                
                double dec = wi * (this->branches->at(j)->totalTaxaCounts[currPartition]) + wj * (this->branches->at(i)->totalTaxaCounts[currPartition]);
                this->pairs[i][j][currPartition] -= dec;
                
                this->sumPairsBranch[i][1 - currPartition] += inc;
                this->sumPairsBranch[j][1 - currPartition] += inc;
                this->sumPairsBranch[i][currPartition] -= dec;
                this->sumPairsBranch[j][currPartition] -= dec;

                this->sumPairs[1 - currPartition] += inc;
                this->sumPairs[currPartition] -= dec;

                // this->subs[i][j][currPartition] -= wi * wj;
                // this->subs[i][j][1 - currPartition] += wi * wj;
            }

            if(currPartition == 1){
                this->pairsBFromSingleBranch[i] -= (this->branches->at(i)->totalTaxaCounts[1]) * wi;
                this->sumPairsBSingleBranch -= (this->branches->at(i)->totalTaxaCounts[1]) * wi; 

            }
            else{

                this->pairsBFromSingleBranch[i] += (this->branches->at(i)->totalTaxaCounts[1] - wi) * wi;
                this->sumPairsBSingleBranch += (this->branches->at(i)->totalTaxaCounts[1] - wi) * wi;                     
            }
        }

        this->totalTaxa[1 - currPartition] += this->dummyTaxaWeightsIndividual[dummyIndex];
        this->totalTaxa[currPartition] -= this->dummyTaxaWeightsIndividual[dummyIndex];

        this->nonQuartets = calcNonQuartets();
        
    }


    virtual double** gainRealTaxa(double originalScore, double multiplier) {
        // double** gainsOfBranches = new double[this->branches->size()][2];
        double** gainsOfBranches = new double*[this->branches->size()];
        for(int i = 0; i < this->branches->size(); ++i){
            gainsOfBranches[i] = new double[2];
            gainsOfBranches[i][0] = 0;
            gainsOfBranches[i][1] = 0;
        }
        for(int i = 0; i < branches->size(); ++i){
            for(int p = 0; p < 2; ++p){
                if(this->branches->at(i)->realTaxaCounts[p] > 0){
                    this->swapRealTaxon(i, p);
                    // this->branches->at(i)->swapRealTaxa(p);
                    gainsOfBranches[i][p] = multiplier * (this->score() - originalScore);
                    this->undoSwapRealTaxon(i, 1 - p);
                    // this->branches->at(i)->swapRealTaxa(1 - p);
                }
            }
        }
        return gainsOfBranches;
    }


    virtual void gainDummyTaxa(double originalScore, double multiplier, double* dummyTaxaGains) {
        for(int i = 0; i < this->nDummyTaxa; ++i){
            int currPartition = this->dummyTaxaPartition[i];
            this->swapDummyTaxon(i, currPartition);
            // for(int j = 0; j < this->branches->size(); ++j){
            //     this->branches->at(j)->swapDummyTaxon(i, currPartition);
            // }
            dummyTaxaGains[i] += multiplier * (this->score() - originalScore);
            this->undoSwapDummyTaxon(i, 1 - currPartition);
            // for(int j = 0; j < this->branches->size(); ++j){
            //     this->branches->at(j)->swapDummyTaxon(i, 1 - currPartition);
            // }
        }
    }

    virtual void batchTransferRealTaxon(int branchIndex, int netTranser){
        // negative if transfering from 1 to 0
        // positive if transfering from 0 to 1

        int currPartition = netTranser > 0 ? 0 : 1;
        
        netTranser = abs(netTranser);

        this->nonQuartets -= changeAmount(branchIndex);

        for(int i = 0; i < this->branches->size(); ++i){
            if(branchIndex == i){
                this->sumPairsBranch[i][1 - currPartition] += netTranser * (this->totalTaxa[1-currPartition] - this->branches->at(i)->totalTaxaCounts[1-currPartition]);
                this->sumPairsBranch[i][currPartition] -= netTranser * (this->totalTaxa[currPartition] - this->branches->at(i)->totalTaxaCounts[currPartition]);

                this->sumPairs[1 - currPartition] += netTranser * (this->totalTaxa[1-currPartition] - this->branches->at(i)->totalTaxaCounts[1-currPartition]);
                this->sumPairs[currPartition] -= netTranser * (this->totalTaxa[currPartition] - this->branches->at(i)->totalTaxaCounts[currPartition]);

            }
            else{
                int mni = branchIndex > i ? i : branchIndex;
                int mxi = branchIndex > i ? branchIndex : i;


                this->pairs[mni][mxi][currPartition] -= netTranser * this->branches->at(i)->totalTaxaCounts[currPartition];
                this->pairs[mni][mxi][1 - currPartition] += netTranser * this->branches->at(i)->totalTaxaCounts[1 - currPartition];


                this->sumPairsBranch[i][currPartition] -= netTranser * this->branches->at(i)->totalTaxaCounts[currPartition];
                this->sumPairsBranch[i][1 - currPartition] += netTranser * this->branches->at(i)->totalTaxaCounts[1 - currPartition];
            }
            

        }

        if(currPartition == 1){
            pairsBFromSingleBranch[branchIndex] -= netTranser * (this->branches->at(branchIndex)->totalTaxaCounts[1] - netTranser) + (netTranser * (netTranser - 1)) / 2;
            this->sumPairsBSingleBranch -= netTranser * (this->branches->at(branchIndex)->totalTaxaCounts[1] - netTranser) + (netTranser * (netTranser - 1)) / 2;                
        }
        else{
            pairsBFromSingleBranch[branchIndex] += netTranser * (this->branches->at(branchIndex)->totalTaxaCounts[1]) + (netTranser * (netTranser - 1)) / 2;
            this->sumPairsBSingleBranch += netTranser * (this->branches->at(branchIndex)->totalTaxaCounts[1]) + (netTranser * (netTranser - 1)) / 2;
        }
        this->totalTaxa[currPartition] -= netTranser;
        this->totalTaxa[1 - currPartition] += netTranser;
        this->nonQuartets += changeAmount(branchIndex);

        // this->branches->at(branchIndex)->batchTransferRealTaxon(netTranser);
        
    }

    
    ~NumSatCalculatorNodeE(){
        delete[] pairsBFromSingleBranch;
        for(int i = 0; i < branches->size(); ++i){
            delete[] sumPairsBranch[i];
        }
        delete[] sumPairsBranch;
        for(int i = 0; i < branches->size(); ++i){
            for(int j = 0; j < branches->size(); ++j){
                delete[] pairs[i][j];
            }
            delete[] pairs[i];
        }
        delete[] pairs;
    }


    
};

