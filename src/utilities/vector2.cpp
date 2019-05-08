#include "geometrycentral/utilities/vector2.h"
#include <geometrycentral/utilities/utilities.h>

#include <cmath>
#include <iostream>

namespace geometrycentral {

std::ostream& operator<<(std::ostream& output, const Vector2& v) {
  output << "<" << v.x << ", " << v.y << ">";
  return output;
}

} // namespace geometrycentral
