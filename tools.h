#ifndef __tools_h__
#define __tools_h__

#include <cstdlib>
#include <string>
#include <iostream>

template<class T>
inline const T MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

template<class T>
inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline void SWAP(T& a, T& b) { T c(a); a=b; b=c; }

inline void error(const std::string error_text)
{
    std::cerr << error_text << std::endl;
	exit(1);
}

inline long bitlen(unsigned long n) {
    for(long i=0;; i++, n>>=1)
        if(n==0) return i;
}

inline long weight(unsigned long n) {// number of ones
    for(long i=0;; n&1 && i++, n>>=1)
        if(n==0) return i;
}

long gcd(long, long);
long gcd(long&, long&, long, long);
long InvMod(long, long);
long PowerMod(long, long, long);
void ContFrac(long&, long&, long, long, long=0);

#define frand() (rand()/((double)RAND_MAX+1))
#define PI 3.141592653589793
#define RSQRT2 0.707106781186548// 1/sqrt(2)

#endif // __tools_h__