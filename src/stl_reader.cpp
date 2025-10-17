#include "stl_reader.hpp"
#include <fstream>
#include <sstream>
#include <cctype>
#include <iostream>

using Eigen::Vector3d;

namespace io {

    static bool isBinarySTL(std::ifstream& in) {
        char header[80];
        in.read(header, 80);
        if (!in) return false;

        // "solid" at start -> probably ASCII
        std::string head(header, header + 5);
        if (head == "solid") {
            in.seekg(0, std::ios::beg);
            return false;
        }
        return true;
    }

    bool loadSTL(const std::string& filename, Mesh& mesh) {
        std::ifstream in(filename, std::ios::binary);
        if (!in) {
            std::cerr << "Cannot open " << filename << "\n";
            return false;
        }

        mesh.vertices.clear();
        mesh.faces.clear();

        if (isBinarySTL(in)) {
            uint32_t triCount;
            in.read(reinterpret_cast<char*>(&triCount), 4);

            mesh.vertices.reserve(triCount * 3);
            mesh.faces.reserve(triCount);

            for (uint32_t i = 0; i < triCount; ++i) {
                float n[3];
                float v[9];
                uint16_t attr;
                in.read(reinterpret_cast<char*>(n), 12);   // normal
                in.read(reinterpret_cast<char*>(v), 36);   // 3 vertices
                in.read(reinterpret_cast<char*>(&attr), 2);

                std::array<int, 3> ids;
                for (int j = 0; j < 3; ++j) {
                    Vertex vert;
                    vert.pos = Eigen::Vector3d(v[j * 3], v[j * 3 + 1], v[j * 3 + 2]);
                    vert.removed = false;
                    ids[j] = static_cast<int>(mesh.vertices.size());
                    mesh.vertices.push_back(vert);
                }
                mesh.faces.push_back({ ids, false });
            }
        }
        else {
            // ASCII STL
            in.seekg(0);
            std::string line;
            std::array<int, 3> ids;
            int count = 0;

            while (std::getline(in, line)) {
                std::istringstream ss(line);
                std::string word;
                ss >> word;
                if (word == "vertex") {
                    double x, y, z;
                    ss >> x >> y >> z;
                    Vertex vert;
                    vert.pos = Eigen::Vector3d(x, y, z);
                    vert.removed = false;
                    ids[count++] = static_cast<int>(mesh.vertices.size());
                    mesh.vertices.push_back(vert);

                    if (count == 3) {
                        mesh.faces.push_back({ ids, false });
                        count = 0;
                    }
                }
            }
        }

        std::cout << "Loaded STL: " << mesh.vertices.size()
            << " vertices, " << mesh.faces.size() << " faces\n";
        return true;
    }

} // namespace io

