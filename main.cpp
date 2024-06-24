#include "Tree.h"
#include "Taxon.h"
#include "DataContainer.h"
#include "Consensus.h"
#include "PerLevelDS.h"


int RealTaxon::count;
int DummyTaxon::idCounter;


int main(){
    string path = "input.txt";
    string consPath = "input.txt";
    RealTaxon::count = 0;
    DummyTaxon::idCounter = 0;

    GeneTrees* geneTrees = new GeneTrees(path);

    geneTrees->readTaxaNames();

    geneTrees->readGeneTrees();

    auto dc = geneTrees->createDataContainer();

    ConsensusTreePartition* consensusTreePartition = new ConsensusTreePartition(
        consPath,
        geneTrees->taxaMap,
        dc
    );

    auto y = consensusTreePartition->makePartition(
        geneTrees->taxa,
        new vector<DummyTaxon*>(),
        0
    );

    delete geneTrees;

    // for(int i = 0; i < dc->nTaxa; ++i ){
    //     cout << dc->taxa->at(i)->label << " " << y->realTaxonPartition[i] << endl;
    // }

    auto x = new TaxaPerLevelWithPartition(
        dc->taxa,
        new vector<DummyTaxon*>(),
        y->realTaxonPartition,
        y->dummyTaxonPartition,
        dc->nTaxa
    );


    auto book = new BookKeepingPerLevel(
        dc,
        x,
        0
    );


    double** realTaxaGains = new double*[dc->nTaxa];
    for(int i = 0; i < dc->nTaxa; ++i){
        realTaxaGains[i] = new double[dc->nTaxa];
        for(int j = 0; j < dc->nTaxa; ++j){
            realTaxaGains[i][j] = 0;
        }
    }

    double* dummyTaxaGains = new double[0];

    book->calculateScore();
    book->calculateScoreAndGains(realTaxaGains, dummyTaxaGains);


    

    // for(auto it = geneTrees->realTaxons->begin(); it != geneTrees->realTaxons->end(); ++it){
    //     cout << it->first << endl;
    // }
}