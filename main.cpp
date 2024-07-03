#include <chrono>


#include "Tree.h"
#include "Taxon.h"
#include "DataContainer.h"
#include "Consensus.h"
#include "PerLevelDS.h"
#include "QFM.h"


int RealTaxon::count;
int DummyTaxon::idCounter;
double QFM::EPS = 1e-3 ;


void initMemory(DataContainer* dc){
    for(int i = 0; i < dc->realTaxaPartitionNodes->size(); ++i){
        PartitionNode* p = dc->realTaxaPartitionNodes->at(i);
        p->data = new Data*[N_THREADS];
        for(int j = 0; j < N_THREADS; ++j){
            p->data[j] = new Data();
            p->data[j]->gainsForSubTree = new double[2];
            p->data[j]->branch = new Branch(0);
        }
    }
    int sz = dc->topSortedPartitionNodes->size();
    for(int i = sz - 1; i > -1; --i){
        PartitionNode* p = dc->topSortedPartitionNodes->at(i);
        if(p->isLeaf){
            continue;
        }
        p->data = new Data*[N_THREADS];
        for(int j = 0; j < N_THREADS; ++j){
            p->data[j] = new Data();
            p->data[j]->gainsForSubTree = new double[2];
            p->data[j]->branch = new Branch(0);
        }
    }

    for(PartitionByTreeNode* p : *dc->partitionsByTreeNodes){
        p->scoreCalculator = new NumSatCalculatorNode*[N_THREADS];
        for(int i = 0; i < N_THREADS; ++i){
            vector<Branch*>* b = new vector<Branch*>(p->partitionNodes->size());
            for(int j = 0; j < p->partitionNodes->size(); ++j){
                b->at(j) = p->partitionNodes->at(j)->data[i]->branch;
            }
            if(p->partitionNodes->size() > 3){
                p->scoreCalculator[i] = new NumSatCalculatorNodeE(b, new int[0]);
            }
            else{
                p->scoreCalculator[i] = new NumSatCalculatorBinaryNode(b, new int[0]);            
            }
        }
    }
}


void tc3(){
        // GeneTrees trees = new GeneTrees("./input/gtree_11tax_est_5genes_R1.tre");
    GeneTrees* trees = new GeneTrees("./input.txt");
    trees->readTaxaNames();
    trees->readGeneTrees();
    // DataContainer dc = trees->readGeneTrees(null);
    DataContainer* dc = trees->createDataContainer();

    // SubProblemsQueue.setInstance(dc, 1, new RandPartition());


    // RealTaxon** rt = new RealTaxon*[2];
    vector<RealTaxon*>* rt = new vector<RealTaxon*>(2);
    // rt[0] = trees->taxaMap.get("3");
    // rt[1] = trees->taxaMap.get("4");
    rt->at(0) = trees->taxaMap->find("3")->second;
    rt->at(1) = trees->taxaMap->find("4")->second;
    
    // 11, {5, 6, 7, 8, 9, 10}
    // RealTaxon[] dt0r = new RealTaxon[6];
    // dt0r[0] = trees->taxaMap.get("5");
    // dt0r[1] = trees->taxaMap.get("6");
    // dt0r[2] = trees->taxaMap.get("7");
    // dt0r[3] = trees->taxaMap.get("8");
    // dt0r[4] = trees->taxaMap.get("9");
    // dt0r[5] = trees->taxaMap.get("10");

    // RealTaxon** dt0r = new RealTaxon*[6];
    vector<RealTaxon*>* dt0r = new vector<RealTaxon*>(6);
    // dt0r[0] = trees->taxaMap->find("5")->second;
    // dt0r[1] = trees->taxaMap->find("6")->second;
    // dt0r[2] = trees->taxaMap->find("7")->second;
    // dt0r[3] = trees->taxaMap->find("8")->second;
    // dt0r[4] = trees->taxaMap->find("9")->second;
    // dt0r[5] = trees->taxaMap->find("10")->second;
    dt0r->at(0) = trees->taxaMap->find("5")->second;
    dt0r->at(1) = trees->taxaMap->find("6")->second;
    dt0r->at(2) = trees->taxaMap->find("7")->second;
    dt0r->at(3) = trees->taxaMap->find("8")->second;
    dt0r->at(4) = trees->taxaMap->find("9")->second;
    dt0r->at(5) = trees->taxaMap->find("10")->second;



    // DummyTaxon internalDT = new DummyTaxon(dt0r, new DummyTaxon[0]);
    DummyTaxon* internalDT = new DummyTaxon(dt0r, new vector<DummyTaxon*>());

    // dt0r = new RealTaxon[1];
    dt0r = new vector<RealTaxon*>(1);
    // dt0r[0] = trees->taxaMap.get("11");
    dt0r->at(0) = trees->taxaMap->find("11")->second;

    // DummyTaxon[] a = new DummyTaxon[1];
    vector<DummyTaxon*>* a = new vector<DummyTaxon*>(1);
    // a[0] = internalDT;
    a->at(0) = internalDT;

    DummyTaxon* dt0 = new DummyTaxon(dt0r,a);

    // RealTaxon[] dt1r = new RealTaxon[2];
    // dt1r[0] = trees->taxaMap.get("1");
    // dt1r[1] = trees->taxaMap.get("2");
    vector<RealTaxon*>* dt1r = new vector<RealTaxon*>(2);
    dt1r->at(0) = trees->taxaMap->find("1")->second;
    dt1r->at(1) = trees->taxaMap->find("2")->second;

    // DummyTaxon dt1 = new DummyTaxon(dt1r, new DummyTaxon[0]);
    DummyTaxon* dt1 = new DummyTaxon(dt1r, new vector<DummyTaxon*>());

    // a = new DummyTaxon[2];
    // a[0] = dt0;
    // a[1] = dt1;
    a = new vector<DummyTaxon*>(2);
    a->at(0) = dt0;
    a->at(1) = dt1;


    // int[] realTaxaPartition = new int[2];
    // realTaxaPartition[0] = 0;
    // realTaxaPartition[1] = 0;
    int* realTaxaPartition = new int[2];
    realTaxaPartition[0] = 0;
    realTaxaPartition[1] = 0;


    // int[] dummyTaxaPartition = new int[2];
    // dummyTaxaPartition[0] = 1;
    // dummyTaxaPartition[1] = 1;
    int* dummyTaxaPartition = new int[2];
    dummyTaxaPartition[0] = 1;
    dummyTaxaPartition[1] = 1;



    TaxaPerLevelWithPartition* taxa = new TaxaPerLevelWithPartition(rt, a, realTaxaPartition, dummyTaxaPartition, 11);

    BookKeepingPerLevel* bookKeepingPerLevel = new BookKeepingPerLevel(dc, taxa, 0);


    // double[][] rtGains = new double[2][2];
    // double[] dtGains = new double[2];
    double** rtGains = new double*[2];
    for(int i = 0; i < 2; ++i){
        rtGains[i] = new double[2];
        rtGains[i][0] = 0;
        rtGains[i][1] = 0;
    }
    double* dtGains = new double[2];
    dtGains[0] = 0;
    dtGains[1] = 0;


    // System.out.println(bookKeepingPerLevel.calculateScoreAndGains(rtGains, dtGains));
    cout << bookKeepingPerLevel->calculateScoreAndGains(rtGains, dtGains) << "\n";

    // for(int i = 0; i < rtGains.length; ++i){
    //     System.out.println("Taxon " + (i+1) + ": " + rtGains[i][0] + ", " + rtGains[i][1]);
    // }
    // for(int i = 0; i < dtGains.length; ++i){
    //     System.out.println("Dummy " + (i+1) + ": " + dtGains[i]);
    // }

    for(int i = 0; i < 2; ++i){
        cout << "Taxon " << (i+1) << ": " << rtGains[i][0] << ", " << rtGains[i][1] << "\n";
    }
    for(int i = 0; i < 2; ++i){
        cout << "Dummy " << (i+1) << ": " << dtGains[i] << "\n";
    }

    // rtGains = new double[2][2];
    // dtGains = new double[2];
    rtGains = new double*[2];
    for(int i = 0; i < 2; ++i){
        rtGains[i] = new double[2];
        rtGains[i][0] = 0;
        rtGains[i][1] = 0;
    }
    dtGains = new double[2];
    dtGains[0] = 0;
    dtGains[1] = 0;

    // BookKeepingPerLevelDC bookDc = new BookKeepingPerLevelDC(dc, taxa, 0);
    BookKeepingPerLevel* book = new BookKeepingPerLevel(dc, taxa, 0);

    // System.out.println( "score : " + bookDc.calculateScoreAndGains(rtGains, dtGains));
    cout << book->calculateScoreAndGains(rtGains, dtGains) << "\n";

    for(int i = 0; i < 2; ++i){
        cout << "Taxon " << (i+1) << ": " << rtGains[i][0] << ", " << rtGains[i][1] << "\n";
    }
    for(int i = 0; i < 2; ++i){
        cout << "Dummy " << (i+1) << ": " << dtGains[i] << "\n";
    }
    // for(int i = 0; i < rtGains.length; ++i){
    //     System.out.println("Taxon " + (i+1) + ": " + rtGains[i][0] + ", " + rtGains[i][1]);
    // }
    // for(int i = 0; i < dtGains.length; ++i){
    //     System.out.println("Dummy " + (i+1) + ": " + dtGains[i]);
    // }
}



int main1(int argc, char** argv){
    string path, consPath, output;
    if(argc < 4){
        cout << "Usage: ./main <gene_tree_file> <consensus_tree_file> <output_file>\n";
        return 1;
    }
    path = argv[1];
    consPath = argv[2];
    output = argv[3];

    // string path = "input.txt";
    // string consPath = "input.txt";
    // string output = "output.txt";

    auto start = chrono::high_resolution_clock::now();

    RealTaxon::count = 0;
    DummyTaxon::idCounter = 0;

    GeneTrees* geneTrees = new GeneTrees(path);

    geneTrees->readTaxaNames();

    geneTrees->readGeneTrees();

    auto dc = geneTrees->createDataContainer();

    initMemory(dc);

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

    

    auto x = new TaxaPerLevelWithPartition(
        dc->taxa,
        new vector<DummyTaxon*>(),
        y->realTaxonPartition,
        y->dummyTaxonPartition,
        dc->nTaxa
    );

    SolutionNode* root = new SolutionNode();

    QFM::recurse(
        dc,
        x,
        0,
        root,
        0,
        consensusTreePartition
    );

    cout << "FM Done\n";

    // cout << SolutionTree::createTree(root) << "\n";
    auto tree = SolutionTree::createTree(root);

    cout << "Tree Created\n";
    cout << tree->getNewickFormat() << "\n";

    ofstream out(output);
    out << tree->getNewickFormat();

    out.close();    

    auto stop = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);


    cout << "Time taken: " << duration.count() << " seconds\n";

    // auto book = new BookKeepingPerLevel(
    //     dc,
    //     x,
    //     0
    // );


    // double** realTaxaGains = new double*[dc->nTaxa];
    // for(int i = 0; i < dc->nTaxa; ++i){
    //     realTaxaGains[i] = new double[dc->nTaxa];
    //     for(int j = 0; j < dc->nTaxa; ++j){
    //         realTaxaGains[i][j] = 0;
    //     }
    // }

    // double* dummyTaxaGains = new double[0];

    // cout << book->calculateScore() << "\n";
    // book->calculateScoreAndGains(realTaxaGains, dummyTaxaGains);

    // for(int i = 0; i < dc->nTaxa; ++i ){
    //     cout << dc->taxa->at(i)->label << " " << y->realTaxonPartition[i] << "\n";
    //     cout << realTaxaGains[i][0] << " " << realTaxaGains[i][1] << "\n";
    // }
    

    // for(auto it = geneTrees->realTaxons->begin(); it != geneTrees->realTaxons->end(); ++it){
    //     cout << it->first << "\n";
    // }

    return 0;
}


int main(int argc, char** argv){
    // tc3();
    main1(argc, argv);
    return 0;

}