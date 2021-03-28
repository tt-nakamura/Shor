#ifndef __QComp_h__
#define __QComp_h__

#include "QBits.h"

void operator+=(const QBits&, const QBits&);
void operator-=(const QBits&, const QBits&);

void negate(const QBits&);
void add(const QBits&, long, long=0);
void sub(const QBits&, long, long=0);
void AddMod(const QBits&, long, long, long=0);
void SubMod(const QBits&, long, long, long=0);
void MulAddMod(const QBits&, const QBits&, long, long, long=0);
void MulSubMod(const QBits&, const QBits&, long, long, long=0);
void PowMulMod(const QBits&, long, const QBits&, long);
void swap(const QBits&, const QBits&, long=0);

inline void operator+=(const QBits& a, long b) { add(a,b); }
inline void operator-=(const QBits& a, long b) { sub(a,b); }
inline void operator++(const QBits& q) { q+=1; }
inline void operator--(const QBits& q) { q-=1; }
inline void operator++(const QBits& q, int) { q+=1; }
inline void operator--(const QBits& q, int) { q-=1; }

long shor(long&, long);

#endif // __QComp_h__