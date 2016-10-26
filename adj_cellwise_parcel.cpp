#include<stdio.h>
#include <array>
#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>

extern "C"
{
  void adj_cellwise_parcel_(float *, float *, float *, float *);
}

void adj_cellwise_parcel_(float *rho, float *th, float *rv, float *rc)
{
  libcloudphxx::blk_1m::opts_t<float> opts;
  opts.conv=false;
  opts.accr=false;
  opts.sedi=false;
  opts.r_eps = 1e-6;
  float rr = 0;
  float dt = 1;
  
  std::array<float, 1> Arho({*rho});
  std::array<float, 1> Ath({*th});
  std::array<float, 1> Arv({*rv});
  std::array<float, 1> Arc({*rc});
  std::array<float, 1> Arr({rr});
  libcloudphxx::blk_1m::adj_cellwise<float>(opts, Arho, Ath, Arv, Arc, Arr, dt);
  *th = Ath[0];
  *rv = Arv[0];
  *rc = Arc[0];
}
