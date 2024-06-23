#pragma once

#include <string>

#include "Config.h"

using namespace std;


class RealTaxon{
public:
    static int count;

    int id;
    string label;

    RealTaxon(int i, string label){
        this->id = i;
        this->label = label;
    }

    RealTaxon(string label){
        this->id = count++;
        this->label = label;
    }
};

class DummyTaxon {
public:
    static int idCounter;
    vector<RealTaxon*>* realTaxa;
    vector<DummyTaxon*>* dummyTaxa;
    vector<RealTaxon*>* flattenedRealTaxa;

    int taxonCount;

    int realTaxonCount;

    int flattenedTaxonCount;

    int id;

    int nestedLevel;


    DummyTaxon(vector<RealTaxon*>* rts, vector<DummyTaxon*>* dts){
        this->realTaxonCount = rts->size();
        this->taxonCount = this->realTaxonCount + dts->size();
        this->realTaxa = rts;
        this->dummyTaxa = dts;
        this->nestedLevel = 0;

        for(auto x : *dts){
            this->flattenedTaxonCount += x->flattenedTaxonCount;
            this->nestedLevel = max(this->nestedLevel, x->nestedLevel);
        }
        
        int i = 0;
        this->flattenedTaxonCount += this->realTaxonCount;
        this->flattenedRealTaxa = new vector<RealTaxon*>(this->flattenedTaxonCount);

        for(auto x : *rts){
            this->flattenedRealTaxa->at(i) = x;
            i++;
        }

        for(auto x : *dts){
            for(auto y : *x->flattenedRealTaxa){
                this->flattenedRealTaxa->at(i) = y;
                i++;
            }
        }

        this->id = idCounter++;
        this->nestedLevel += 1;
    }

    void calcDivCoeffs(ScoreNormalizationType normalizationType, double* coeffs, double multiplier){
        if(normalizationType == NO_NORMALIZATION){
            for(auto x : *this->realTaxa){
                coeffs[x->id] = 1;
            }
        }
        else if(normalizationType == FLAT_NORMALIZATION){
            for(auto x : *this->realTaxa){
                coeffs[x->id] = this->flattenedTaxonCount * multiplier;
            }
        }
        else if(normalizationType == NESTED_NORMALIZATION){
            double sz = this->realTaxonCount + this->dummyTaxa->size();
            for(auto x : *this->realTaxa){
                coeffs[x->id] = sz * multiplier;
            }
            for(auto x : *this->dummyTaxa){
                x->calcDivCoeffs(normalizationType, coeffs, sz * multiplier);
            }
        }
    }
};
