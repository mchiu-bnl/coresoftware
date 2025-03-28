#ifndef TRIGGER_MINIMUMBIASINFOV1_H
#define TRIGGER_MINIMUMBIASINFOV1_H

#include "MinimumBiasInfo.h"

#include <iostream>

class MinimumBiasInfov1 : public MinimumBiasInfo
{
 public:
  MinimumBiasInfov1() = default;
  virtual ~MinimumBiasInfov1() override = default;

  void identify(std::ostream &os = std::cout) const override;

  void Reset() override { _isMinBias = false; }

  int isValid() const override { return 1; }

  PHObject *CloneMe() const override { return new MinimumBiasInfov1(*this); }

  void CopyTo(MinimumBiasInfo *mbinfo) override;

  void setIsAuAuMinimumBias(bool is_min_bias) override { _isMinBias = is_min_bias; }
  bool isAuAuMinimumBias() const override { return _isMinBias; }

 private:
  bool _isMinBias{false};

  ClassDefOverride(MinimumBiasInfov1, 1)
};

#endif
