#pragma once

void addArrayToFirst(double* a, double* b, int n){
    for (int i = 0; i < n; ++i) {
        a[i] += b[i];
    }
}