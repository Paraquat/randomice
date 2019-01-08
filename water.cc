#include "water.h"

Water::Water()
{
}

Water::~Water()
{
}

Water::Water(const Water& w)
{
  O = w.O;
  H1 = w.H1;
  H2 = w.H2;
  H3 = w.H3;
  H4 = w.H4;
}

Water& Water::operator= (const Water& w)
{
  if (this == &w) return *this;
  O = w.O;
  H1 = w.H1;
  H2 = w.H2;
  H3 = w.H3;
  H4 = w.H4;
  return *this;
}

bool Water::operator== (const Water& w)
{
  if (O == w.O) return true;
  else return false;
}

Water::Water(int o, int h1, int h2, int h3, int h4)
{
  O = o;
  H1 = h1;
  H2 = h2;
  H3 = h3;
  H4 = h4;
}

void Water::add_hbond(int hb)
{
  hbonds.push_back(hb);
}