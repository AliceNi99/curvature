#pragma once
#include "mesh.hpp"

namespace io {
    bool loadSTL(const std::string& filename, Mesh& mesh);
}
