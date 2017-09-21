#include "hbond.h"

Hbond::Hbond()
{
}

Hbond::~Hbond()
{
}

Hbond::Hbond(const Hbond& h)
{
  O1 = h.O1;
  O2 = h.O2;
  H1 = h.H1;
  H2 = h.H2;
  W1 = h.W1;
  W2 = h.W2;
  w = h.w;
}

Hbond& Hbond::operator= (const Hbond& h)
{
  if (this == &h) return *this;
  O1 = h.O1;
  O2 = h.O2;
  H1 = h.H1;
  H2 = h.H2;
  W1 = h.W1;
  W2 = h.W2;
  w = h.w;
  return *this;
}

Hbond::Hbond(int o1, int o2, int h1, int h2, int w1, int w2)
{
  O1 = o1;
  O2 = o2;
  H1 = h1;
  H2 = h2;
  W1 = w1;
  W2 = w2;
  w.insert(o1);
  w.insert(o2);
}