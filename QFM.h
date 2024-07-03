#pragma once

#include "DataContainer.h"
#include "PerLevelDS.h"
#include "Config.h"
#include "MakePartition.h"
#include "Taxon.h"
#include "SolutionTree.h"


class Swap{
public:
    int index;
    bool isDummy;
    double gain;

    Swap(int i, bool id, double g){
        this->index = i;
        this->isDummy = id;
        this->gain = g;
    }
};

class QFM {
public:
    static double EPS;
    static void recurse(
        DataContainer* dc,
        TaxaPerLevelWithPartition* taxa,
        int tid,
        SolutionNode* solnNode, 
        int level,
        IMakePartition* initPartition
    ){

        int itrCount = 0;


        BookKeepingPerLevel* book = new BookKeepingPerLevel(dc, taxa, tid);
        
        // cout << "book address " << book << "\n";

        cout << "Level: " << level << "\n";

        while(oneInteration(book) ){
            itrCount++;
            if(itrCount > MAX_ITERATION){
                // System.out.println("Max iteration reached");
                cout << "Max iteration reached" << "\n";
                break;
            }
            
        }

        cout << "#iterations : " << itrCount << "\n";
        // cout << "here\n";
        // cout << "book address " << book << "\n";

        auto x = book->divide(initPartition); 

        SolutionNode** children = new SolutionNode*[2];
        


        // cout << "Divided\n";
        // taxa->printTaxa();


        delete book;
        delete taxa;

        int i = 0;
        int dummyIds[2];
        
        for(i = 0; i < 2; ++i){
            
            auto taxaWPart = x[i];
            dummyIds[i] = taxaWPart->dummyTaxa->at(taxaWPart->dummyTaxonCount - 1)->id;
            children[i] = new SolutionNode();
            if(taxaWPart->smallestUnit){
                children[i]->tree = taxaWPart->createStar();
            }
            else{
                // trees[i] = recurse(childBooks, children[i]);
                // recurse(childBooks, children[i], level + 1, initPartition);
                // SubProblemsQueue.instance.addItem(new Item(taxaWPart, children[i], level + 1));
                recurse(
                    dc,
                    taxaWPart,
                    0,
                    children[i],
                    level + 1,
                    initPartition
                );
            }
        }

        solnNode->leftDTid = dummyIds[0];
        solnNode->rightDTid = dummyIds[1];
        solnNode->left = children[0];
        solnNode->right = children[1];
        
        // SubProblemsQueue.instance.free(tid);

    }

    
    

    static Swap* swapMax(BookKeepingPerLevel* book, double** rtGains, double* dtGains, bool* rtLocked, bool* dtLocked){

        int maxGainIndex = -1;
        double maxGain = 0;

        for(int i = 0; i < book->taxaPerLevel->realTaxonCount; ++i){
            if(rtLocked[i]) continue;
            int partition = book->taxaPerLevel->inWhichPartitionRealTaxonByIndex(i);
            if((book->taxaPerLevel->getTaxonCountInPartition(partition) > 2) ){
                if(maxGainIndex == -1){
                    maxGain = rtGains[i][partition];
                    maxGainIndex = i;
                }
                else if(maxGain < rtGains[i][partition]){
                    maxGain = rtGains[i][partition];
                    maxGainIndex = i;
                }
            }
        }

        bool dummyChosen = false;

        for(int i = 0; i < book->taxaPerLevel->dummyTaxonCount; ++i){
            if(dtLocked[i]) continue;
            int partition = book->taxaPerLevel->inWhichPartitionDummyTaxonByIndex(i);
            if(book->taxaPerLevel->getTaxonCountInPartition(partition) > 2 ){
                if(maxGainIndex == -1){
                    maxGain = dtGains[i];
                    maxGainIndex = i;
                    dummyChosen = true;
                }
                else if(maxGain < dtGains[i]){
                    maxGain = dtGains[i];
                    maxGainIndex = i;
                    dummyChosen = true;
                }
            }
        }

        if(maxGainIndex == -1) return NULL;

        book->swapTaxon(maxGainIndex, dummyChosen);
        if(dummyChosen){
            dtLocked[maxGainIndex] = true;
        }
        else{
            rtLocked[maxGainIndex] = true;
        }

        
        return new Swap(maxGainIndex, dummyChosen, maxGain);


    }

    static bool oneInteration(BookKeepingPerLevel* book){
        
        double cg = 0;
        int maxCgIndex = -1;
        double maxCg = 0;


        bool singletonPartition = book->taxaPerLevel->getTaxonCountInPartition(0) == 1  || book->taxaPerLevel->getTaxonCountInPartition(1) == 1;

        bool* rtLocked = new bool[book->taxaPerLevel->realTaxonCount];
        bool* dtLocked = new bool[book->taxaPerLevel->dummyTaxonCount];

        for(int i = 0; i < book->taxaPerLevel->realTaxonCount; ++i){
            rtLocked[i] = false;
        }
        for(int i = 0; i < book->taxaPerLevel->dummyTaxonCount; ++i){
            dtLocked[i] = false;
        }

        double** rtGains;
        double* dtGains;

        // ArrayList<Swap> swaps = new ArrayList<Swap>();
        vector<Swap*> swaps;

        // ArrayList<double> cgs = new ArrayList<double>();

        // System.out.println("Iteration started");

        rtGains = new double*[book->taxaPerLevel->realTaxonCount];
        for(int i = 0; i < book->taxaPerLevel->realTaxonCount; ++i){
            rtGains[i] = new double[2];
        }
        dtGains = new double[book->taxaPerLevel->dummyTaxonCount];

        while(true){
            for(int i = 0; i < book->taxaPerLevel->realTaxonCount; ++i){
                rtGains[i][0] = 0;
                rtGains[i][1] = 0;
            }
            for(int i = 0; i < book->taxaPerLevel->dummyTaxonCount; ++i){
                dtGains[i] = 0;
            }
            
            // cout << "initialized gains\n";

            book->calculateScoreAndGains(rtGains, dtGains);

            // cout << "calculated gains\n";

            // var x = swapMax(book, rtGains, dtGains, rtLocked,dtLocked);
            Swap* x = swapMax(book, rtGains, dtGains, rtLocked,dtLocked);

            // cout << "swapped max\n";
            
            if(x != NULL){
                swaps.push_back(x);
                
                double gain = x->gain;

                // if(gain < 0) break;

                // System.out.println("Gain: " + gain + " CG: " + cg);
                // cout << "Gain: " << gain << " CG: " << cg << "\n";
                cg += gain;

                // cgs.add(cg);

                if(singletonPartition){
                    if(maxCgIndex == -1 ){ // && book->taxas.getTaxonCountInPartition(0) > 1 && book->taxas.getTaxonCountInPartition(1) > 1 ){
                        maxCg = cg;
                        maxCgIndex = swaps.size() - 1;
                    }
                }
                
                if(cg > maxCg && abs(maxCg - cg) > EPS ){ // && book->taxas.getTaxonCountInPartition(0) > 1 && book->taxas.getTaxonCountInPartition(1) > 1 ){
                    maxCg = cg;
                    maxCgIndex = swaps.size() - 1;
                }

                // if(!singletonPartition){
                //     if(cg < 0){
                //         break;
                //     }
                // }

                // if(cg < 0){
                //     break;
                // }
            }
            else{
                break;
            }

            // cout << "one itr of loop\n"; 
            
        }

        delete[] dtGains;
        for(int i = 0; i < book->taxaPerLevel->realTaxonCount; ++i){
            delete[] rtGains[i];
        }
        delete[] rtGains;

        // System.out.println("Cg : " + cg);

        // System.out.println("Iteration ended, maxCgIndex: " + maxCgIndex + " maxCg: " + maxCg);

        // cout << "Cg : " << cg << "\n";
        

        if(maxCgIndex == -1){
            if(swaps.size() != (book->taxaPerLevel->realTaxonCount + book->taxaPerLevel->dummyTaxonCount)){
                for(int i = swaps.size() - 1; i >= 0; --i){
                    // var x = swaps.get(i);
                    Swap* x = swaps[i];
                    // book->swapTaxon(x.index, x.isDummy);
                    if(x->isDummy){
                        book->taxaPerLevel->swapPartitionDummyTaxon(x->index);
                    }
                    else{
                        book->taxaPerLevel->swapPartitionRealTaxon(x->index);
                    }
                }
            }
            for(Swap* s : swaps){
                delete s;
            }
            return false;
        }
        vector<int> rtIds;
        vector<int> dtIds;

        for(int i = swaps.size() - 1; i > maxCgIndex; --i){
            Swap* s = swaps[i];
            if(s->isDummy){
                dtIds.push_back(s->index);
            }
            else{
                rtIds.push_back(s->index);
            }
        }

        if(rtIds.size() > 0){
            if(rtIds.size() > 1){
                book->batchTrasferRealTaxon(rtIds);
            }
            else{
                book->swapRealTaxon(rtIds[0]);
            }

        }
        for(int i : dtIds){
            book->swapDummyTaxon(i);
        }

        for(Swap* s : swaps){
            delete s;
        }

        return true;

    }

};
