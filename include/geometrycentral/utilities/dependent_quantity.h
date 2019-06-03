#pragma once

#include <functional>
#include <vector>

namespace geometrycentral {

// Helper function which clears out a buffer of data. Can be overridden with clearing methods as appropriate.
template <typename D>
void clearBuffer(D* buff){};

class DependentQuantity {

public:
  DependentQuantity(std::vector<DependentQuantity*> dependencies_, std::function<void()> evaluateFunc_)
      : dependencies(dependencies_), evaluateFunc(evaluateFunc_) {}

  virtual ~DependentQuantity(){};

  std::vector<DependentQuantity*> dependencies;
  std::function<void()> evaluateFunc;
  bool computed = false;
  int requireCount = 0;

  // Compute the quantity, if we don't have it already
  void ensureHave();

  // Compute the quantity if we need it and don't have it already
  void ensureHaveIfRequired();

  // Note that something will reqiure this quantity (increments a count of such requirements),
  // and ensure that we have this quantity
  void require();

  // Decrement the count of requirements of this quantity
  void unrequire();

  // Clear out the underlying quantity to reduce memory usage
  virtual void clearIfNotRequired() = 0;
};

// Wrapper class which manages a dependency graph of quantities. Templated on the underlying type of the data.
template <typename D>
class DependentQuantityD : public DependentQuantity {

public:
  DependentQuantityD(){};
  virtual ~DependentQuantityD(){};

  DependentQuantityD(D* dataBuffer_, std::function<void()> evaluateFunc_, std::vector<DependentQuantity*> dependencies_)
      : DependentQuantity(dependencies, evaluateFunc_), dataBuffer(dataBuffer_) {}

  D* dataBuffer = nullptr;

  // Clear out the underlying quantity to reduce memory usage
  virtual void clearIfNotRequired() override;
};

} // namespace geometrycentral

#include "geometrycentral/utilities/dependent_quantity.ipp"
