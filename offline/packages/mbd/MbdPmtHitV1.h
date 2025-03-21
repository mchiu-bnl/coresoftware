#ifndef __MBD_MBDPMTHITV1_H__
#define __MBD_MBDPMTHITV1_H__

#include "MbdPmtHit.h"

#include <cmath>
#include <iostream>
#include <limits>

class MbdPmtHitV1 : public MbdPmtHit
{
 public:
  MbdPmtHitV1() = default;
  ~MbdPmtHitV1() override = default;

  //! Just does a clear
  void Reset() override;

  //! Clear is used by TClonesArray to reset the tower to initial state without calling destructor/constructor
  void Clear(Option_t* = "") override;

  //! PMT number
  Short_t get_pmt() const override { return bpmt; }

  //! Effective Nch in PMT
  Float_t get_q() const override { return bq; }

  //! Time from time channel
  Float_t get_tt() const override { return btt; }
  Float_t get_time() const override { return btt; }

  //! Time from charge channel
  Float_t get_tq() const override { return btq; }

  void set_pmt(const Short_t pmt, const Float_t q, const Float_t tt, const Float_t tq) override
  {
    bpmt = pmt;
    bq = q;
    btt = tt;
    btq = tq;
  }

  //! Prints out exact identity of object
  void identify(std::ostream& out = std::cout) const override;

  //! isValid returns non zero if object contains valid data
  virtual int isValid() const override
  {
    if (std::isnan(get_time())) return 0;
    return 1;
  }

 private:
  //Short_t bpmt{-1};
  //Float_t bq{std::numeric_limits<float>::quiet_NaN()};
  //Float_t btt{std::numeric_limits<float>::quiet_NaN()};
  //Float_t btq{std::numeric_limits<float>::quiet_NaN()};
  Short_t bpmt;
  Float_t bq;
  Float_t btt;
  Float_t btq;

  ClassDefOverride(MbdPmtHitV1, 1)
};

#endif
