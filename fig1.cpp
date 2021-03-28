#include<fstream>
#include "QComp.h"
#include "tools.h"

main() {
    std::ofstream f("fig1.txt");
    
    long n(21),a(2),l(bitlen(n)),i;
    long m(n==15 ? l : l<<1);// Mermin, p82
    QBits x(m),y(l);
    mat_double p;

    srand(100);
    clear(x);
    set(y,1);// y=1
    Hadamar(x);// x=0,1,..,M-1 superposed
    PowMulMod(y,a,x,n);// y*=b^x mod n
    QFT(x);
    measure(y);// select one of p(x)
    probability(p,x);
    for(i=0; i<x.dim(); i++) {
        std::cout << i << ' ' << p[0][i] << std::endl;
        f << i << ' ' << p[0][i] << std::endl;
    }
}
