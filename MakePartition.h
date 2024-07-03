#pragma once
#include <vector>
#include <string>

#include "Taxon.h"


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
