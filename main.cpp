#include "Tree.h"
#include "Taxon.h"


int RealTaxon::count;


int main(){
    string path = "input.txt";

    RealTaxon::count = 0;

    GeneTrees* geneTrees = new GeneTrees(path);

    geneTrees->readTaxaNames();

    geneTrees->readGeneTrees();

    auto dc = geneTrees->createDataContainer();

    delete geneTrees;
    // for(auto it = geneTrees->realTaxons->begin(); it != geneTrees->realTaxons->end(); ++it){
    //     cout << it->first << endl;
    // }
}