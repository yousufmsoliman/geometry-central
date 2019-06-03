namespace geometrycentral {

inline void DependentQuantity::ensureHaveIfRequired() {
  if (requireCount > 0) {
    ensureHave();
  }
}

inline void DependentQuantity::ensureHave() {

  // If the quantity is already populated, early out
  if (computed) {
    return;
  }

  // Resolve all of the dependencies
  for (auto x : dependencies) {
    x->ensureHave();
  }

  // Compute this quantity
  evaluateFunc();

  computed = true;
};

inline void DependentQuantity::require() {
  requireCount++;
  ensureHave();
}

inline void DependentQuantity::unrequire() {
  requireCount--;

  if (requireCount < 0) {
    throw std::logic_error("Quantity was unrequire()'d more than than it was require()'d");
    requireCount = 0;
  }
}

template <typename D>
void DependentQuantityD<D>::clearIfNotRequired() {
  if (requireCount <= 0 && dataBuffer != nullptr) {
    clearBuffer(dataBuffer);
  }
}

} // namespace geometrycentral
