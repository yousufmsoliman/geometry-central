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

// Helper template
// Note: if/when we start using more types in these quantities, we might need to generalize this mechanism. But for the
// current set of uses (scalars and MeshData<>), it works just fine.
namespace {

template <typename T>
void clearBuffer(T* buffer) {
  buffer->clear();
}

template <>
void clearBuffer(double* buffer) {}
template <>
void clearBuffer(size_t* buffer) {}
template <>
void clearBuffer(int* buffer) {}

} // namespace

template <typename D>
void DependentQuantityD<D>::clearIfNotRequired() {
  if (requireCount <= 0 && dataBuffer != nullptr && computed) {
    clearBuffer(dataBuffer);
    computed = false;
  }
}

} // namespace geometrycentral
