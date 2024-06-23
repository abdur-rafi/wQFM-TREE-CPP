#pragma once
#include <vector>
#include <string>

#include "Taxon.h"
#include "Tree.h"

class MakePartitionReturnType{
public:
    int* realTaxonPartition;
    int* dummyTaxonPartition;
    MakePartitionReturnType(int* realTaxonPartition, int* dummyTaxonPartition){
        this->realTaxonPartition = realTaxonPartition;
        this->dummyTaxonPartition = dummyTaxonPartition;
    }
};

class IMakePartition{
public:
    virtual MakePartitionReturnType* makePartition(
        vector<RealTaxon*>* realTaxons,
        vector<DummyTaxon*>* dummyTaxons, 
        int tid
    ) = 0;
};

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




class ConsensusTreePartition : public IMakePartition{
public:
    Tree* tree;
    DataContainer* dc;
    RandomPartition* rp;

    ConsensusTreePartition(string file, map<string, RealTaxon*>* taxaMap, DataContainer* dc){
        this->dc = dc;
        ifstream in(file);
        string line;
        if(!in.is_open()){
            cout << "Consensus File not found" << endl;
            return;
        }
        getline(in, line);
        this->tree = new Tree(line, taxaMap);
        in.close();
        this->rp = new RandomPartition();
    }

    virtual MakePartitionReturnType* makePartition(
        vector<RealTaxon*>* realTaxons,
        vector<DummyTaxon*>* dummyTaxons, 
        int tid
    ){
        return this->rp->makePartition(realTaxons, dummyTaxons, tid);
    }

};


