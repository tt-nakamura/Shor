#include<iostream>
#include "QComp.h"

main() {
    long n(15),d,i,k;
    srand(100);
    for(i=0; i<10; i++) {
        k = shor(d,n);
        if(k==0)     std::cout << "factor found by shor: ";
        else if(k>0) std::cout << "factor found not by shor: ";
        else         std::cout << "factor not found: ";
        std::cout << d << std::endl;
    }
}