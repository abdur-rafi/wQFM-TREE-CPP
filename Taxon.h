#pragma once

#include <string>
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