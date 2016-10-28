#include<stdio.h>
#include <array>
#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adj_cellwise.hpp>
#include <libcloudph++/common/theta_dry.hpp>

extern "C"
{
  void adj_cellwise_parcel_(float *, float *, float *, float *, float *);
}

void adj_cellwise_parcel_(float *rho, float *th, float *rv, float *rc, float *pre)
// total rho, th_std, rv, rc, total pre
{
  libcloudphxx::blk_1m::opts_t<float> opts;
  opts.conv=false;
  opts.accr=false;
  opts.sedi=false;
  opts.r_eps = 1e-6;
  float rr = 0; // no rain!
  float dt = 2; // could be anything?
  
  // calc th_dry
  quantity<si::dimensionless, float> rv_si = *rv;
  auto th_dry = libcloudphxx::common::theta_dry::std2dry((*th) * si::kelvins, rv_si);

  // calc rho_dry
  auto rho_dry = libcloudphxx::common::theta_std::rhod((*pre) * si::pascals, (*th) * si::kelvins, rv_si);

  auto T = libcloudphxx::common::theta_dry::T(th_dry, rho_dry);
  auto pre_calc = libcloudphxx::common::theta_dry::p(rho_dry, rv_si, T);

//  std::cout << "pre: " << (*pre) << "pre_calc: " <<  pre_calc << std::endl;

  std::array<float, 1> Arho({rho_dry / si::kilograms * si::cubic_metres});
  std::array<float, 1> Ath({th_dry / si::kelvins});
  std::array<float, 1> Arv({*rv});
  std::array<float, 1> Arc({*rc});
  std::array<float, 1> Arr({rr});
  libcloudphxx::blk_1m::adj_cellwise<float>(opts, Arho, Ath, Arv, Arc, Arr, dt);
  rv_si = Arv[0];
  *th = libcloudphxx::common::theta_dry::dry2std(Ath[0] * si::kelvins, rv_si) / si::kelvins;
  *rv = Arv[0];
  *rc = Arc[0];
}
