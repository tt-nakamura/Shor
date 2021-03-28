#include "tools.h"

long gcd(long a, long b) {
    if(a<0) a=-a;
    if(b<0) b=-b;
    if(b==0) return a;
    for(;;) {
        if((a%=b)==0) return b;
        if((b%=a)==0) return a;
    }
}

long gcd(long& s, long& t, long a, long b) {// sa+tb=gcd
    bool c(a<0),d(b<0);
    long u(0),v(1),w;
    ldiv_t qr;
    if(c) a=-a;
    if(d) b=-b;
    s=1; t=0;
    while(b) {
        qr = ldiv(a,b);
        a = b;
        b = qr.rem;
        w = u;
        u = s - u*qr.quot;
        s = w;
        w = v;
        v = t - v*qr.quot;
        t = w;
    }
    if(c) s=-s;
    if(d) t=-t;
    return a;
}

long InvMod(long a, long n)// a^{-1} mod n
// assume n>0, 0<=a<n
{
    long s,t;
    if(gcd(s,t,a,n)!=1) error("gcd!=1 in InvMod");
    if(s<0) s+=n;
    return s;
}

long PowerMod(long a, long e, long n)// a^e mod n
// assume n>0, 0<=a<n
{
    if(a==0) return 0;
    if(e==0 || a==1) return 1;
    if(e<0) return PowerMod(InvMod(a,n),-e,n);
    long m(1L<<bitlen(e)-1), x(a);
    for(m>>=1; m; m>>=1) {
        x*=x; x%=n;
        if(e&m) { x*=a; x%=n; }
    }
    return x;
}

void ContFrac(long& p, long& q, long a, long b, long n)
// a/b is expanded as continued fraction to obtain 
// p/q = rational approximation of a/b with largest q<n
//       |p/q - a/b| < 1/q^2, p and q are coprime, q>0
// if n==0, n is set to b
{
    bool c(a<0),d(b<0);
    long s(0),t(1),u;
    ldiv_t qr;
    if(c) a=-a;
    if(d) b=-b;
    if(n==0) n=b;
    p=1; q=0;
    while(b && q<n) {
        qr = ldiv(a,b);
        a = b;
        b = qr.rem;
        u = p;
        p = s + p*qr.quot;
        s = u;
        u = q;
        q = t + q*qr.quot;
        t = u;        
    }
    if(q>=n) { p=s; q=t; }
    if(c^d) p=-p;
}