#include <omp.h>
#include "mesh.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>

   bool Mesh::loadOFF(const std::string &filename) {
        std::ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Failed to open " << filename << "\n";
            return false;
        }

        std::string line;
        // Read header
        std::getline(infile, line);
        if (line != "OFF") {
            std::cerr << "Not a valid OFF file\n";
            return false;
        }

        // Read counts
        size_t nVertices = 0, nFaces = 0, nEdges = 0;
        while (std::getline(infile, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            ss >> nVertices >> nFaces >> nEdges;
            break;
        }

        if (nVertices == 0 || nFaces == 0) {
            std::cerr << "OFF file has no vertices or faces\n";
            return false;
        }

        // Read vertices
        vertices.resize(nVertices);
        size_t vcount = 0;
        while (vcount < nVertices && std::getline(infile, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            double x, y, z;
            ss >> x >> y >> z;
            vertices[vcount].pos = Eigen::Vector3d(x, y, z);

            // Optional color
            double r, g, b;
            if (ss >> r >> g >> b) {
                vertices[vcount].color = Eigen::Vector3d(r, g, b);
                // if colors in [0,255], normalize to [0,1]
                if (r > 1.0 || g > 1.0 || b > 1.0)
                    vertices[vcount].color /= 255.0;
            }

            ++vcount;
        }

        // Read faces
        faces.resize(nFaces);
        size_t fcount = 0;
        while (fcount < nFaces && std::getline(infile, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            int n, i0, i1, i2;
            ss >> n;
            if (n != 3) {
                std::cerr << "Only triangular faces supported\n";
                return false;
            }
            ss >> i0 >> i1 >> i2;
            faces[fcount].v = {i0, i1, i2};
            ++fcount;
        }

        std::cout << "Loaded OFF mesh: " << nVertices << " vertices, "
                  << nFaces << " faces\n";
        return true;
    }


bool Mesh::loadOBJ(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Cannot open OBJ file: " << path << "\n";
        return false;
    }
    vertices.clear();
    faces.clear();

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "v") {
            double x, y, z;
            ss >> x >> y >> z;
            Vertex v;
            v.pos = Eigen::Vector3d(x, y, z);
            vertices.push_back(v);
        }
        else if (tag == "f") {
            Face f;
            for (int i = 0; i < 3; i++) {
                std::string vert;
                ss >> vert;
                size_t slash = vert.find('/');
                int idx = (slash == std::string::npos)
                    ? std::stoi(vert)
                    : std::stoi(vert.substr(0, slash));
                f.v[i] = idx - 1;
            }
            faces.push_back(f);
        }
    }
    vertexFaces.resize(vertices.size());
    for (int fi = 0; fi < faces.size(); ++fi) {
        for (int vi : faces[fi].v)
            vertexFaces[vi].insert(fi);
    }

    std::cout << "Loaded OBJ: " << vertices.size() << " vertices, "
        << faces.size() << " faces\n";
    return true;
}

/*
void Mesh::VVAdj() {
    edges.clear();

    std::set<std::pair<int, int>> edgeSet; // to avoid duplicates

    for (const auto& f : faces) {
        if (f.removed) continue;

        int v0 = f.v[0];
        int v1 = f.v[1];
        int v2 = f.v[2];

        // Add adjacency
        vertices[v0].neighbors.insert(v1);
        vertices[v0].neighbors.insert(v2);

        vertices[v1].neighbors.insert(v0);
        vertices[v1].neighbors.insert(v2);

        vertices[v2].neighbors.insert(v0);
        vertices[v2].neighbors.insert(v1);

        // Add edges (sorted to avoid duplicates)
        std::array<std::pair<int, int>, 3> fEdges = { {
            {std::min(v0,v1), std::max(v0,v1)},
            {std::min(v1,v2), std::max(v1,v2)},
            {std::min(v2,v0), std::max(v2,v0)}
        } };

        for (auto& e : fEdges) {
            if (edgeSet.insert(e).second) { // insert returns true if new
                EdgeCollapse ec;
                ec.v1 = e.first;
                ec.v2 = e.second;
                ec.cost = 0.0; // will compute later
                ec.optimal.setZero();
                edges.push_back(ec);
            }
        }
    }

    std::cout << "Built adjacency and " << edges.size() << " unique edges.\n";
}
    */


void Mesh::VVAdj2() {
    VV.clear();
    VF.clear();
    VV.resize(vertices.size());
    VF.resize(vertices.size());
    for (const auto& f : faces) {
        if (f.removed) continue;

        int v0 = f.v[0];
        int v1 = f.v[1];
        int v2 = f.v[2];

        // Add adjacency
        VV[v0].push_back(v1);
        VV[v0].push_back(v2);
        VF[v0].push_back(&f - &faces[0]);

        VV[v1].push_back(v0);
        VV[v1].push_back(v2);
        VF[v1].push_back(&f - &faces[0]);

        VV[v2].push_back(v0);
        VV[v2].push_back(v1);
        VF[v2].push_back(&f - &faces[0]);
    }
    // Remove neighboring duplicates
    for (auto& neighbors : VV) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }
    //Remove face duplicates
    for (auto& faceList : VF) {
    std::sort(faceList.begin(), faceList.end());                     // sort IDs
    faceList.erase(std::unique(faceList.begin(), faceList.end()),faceList.end());
}

    std::cout << "Built VV and VF adjacency list.\n";
}

void Mesh::EFAdj() {
    EF.clear();
    EF.reserve(faces.size() * 3);

    for (size_t f = 0; f < faces.size(); ++f) {
        const Face &face = faces[f];
        if (face.removed) continue;

        for (int e = 0; e < 3; ++e) {
            int v0 = face.v[e];
            int v1 = face.v[(e + 1) % 3];
            std::pair<int,int> edge(std::min(v0, v1), std::max(v0, v1));
            EF[edge].push_back(static_cast<int>(f));
        }
    }
    
    std::cout << "Built EF adjacency map.\n";
}


void Mesh::curvature_computing_uniform() {
    meanCurvature.resize(vertices.size());
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < vertices.size(); ++i) {
        const auto &vi = vertices[i].pos;
        Eigen::Vector3d Hn = Eigen::Vector3d::Zero();
        double area = VV[i].size(); // uniform weights
        
        for (int j : VV[i]) {
            
            const auto &vj = vertices[j].pos;
            Hn += (vj - vi);

        }
        // calculate the normal to vi by averaging adjacent face normals
        Eigen::Vector3d N = Eigen::Vector3d::Zero();
        for (int f : VF[i]) {
            Eigen::Vector3d v0 = vertices[faces[f].v[0]].pos;
            Eigen::Vector3d v1 = vertices[faces[f].v[1]].pos;
            Eigen::Vector3d v2 = vertices[faces[f].v[2]].pos;
            Eigen::Vector3d faceNormal = (v1 - v0).cross(v2 - v0);
            N += faceNormal;
        }
        auto sign = (N.dot(-Hn) > 0) ? 1.0 : -1.0;
        meanCurvature[i] = sign * Hn.norm() / area;
    }

    std::cout << "mean curvature computed for vertices\n";
}

void Mesh::curvature_computing(const std::string& weight) {
    meanCurvature.resize(vertices.size());
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < vertices.size(); ++i) {
    const auto &vi = vertices[i].pos;
    Eigen::Vector3d Hn = Eigen::Vector3d::Zero();
    Eigen::Vector3d N = Eigen::Vector3d::Zero();
    double area = 0.0;

        for (int j : VV[i]) {
            const auto &vj = vertices[j].pos;
            auto faces_list = EF[{std::min(i,j), std::max(i,j)}]; // faces adjacent to edge (i,j)
            
            // compute cotAlpha, cotBeta using the opposite vertex in each face
            const Face &face0 = faces[faces_list[0]];
            int vk1_index = (face0.v[0] != i && face0.v[0] != j) ? face0.v[0] :
                (face0.v[1] != i && face0.v[1] != j) ? face0.v[1] :
                face0.v[2];
            Eigen::Vector3d vk1 = vertices[vk1_index].pos;
            double cotAlpha = ((vi - vk1).dot(vj - vk1)) / ((vi - vk1).cross(vj - vk1)).norm();

            double cotBeta = 0.0;
            if (faces_list.size() == 2) {
                const Face &face1 = faces[faces_list[1]];
                int vk2_index = (face1.v[0] != i && face1.v[0] != j) ? face1.v[0] :
                                (face1.v[1] != i && face1.v[1] != j) ? face1.v[1] :
                                face1.v[2];
                Eigen::Vector3d vk2 = vertices[vk2_index].pos;
                cotBeta = ((vi - vk2).dot(vj - vk2)) / ((vi - vk2).cross(vj - vk2)).norm();
            }

            Hn += (cotAlpha + cotBeta) * (vj - vi);

            // area contribution (add triangle areas; duplicate counting is OK)
            for (auto f : faces_list) {
                Eigen::Vector3d v0 = vertices[faces[f].v[0]].pos;
                Eigen::Vector3d v1 = vertices[faces[f].v[1]].pos;
                Eigen::Vector3d v2 = vertices[faces[f].v[2]].pos;
                Eigen::Vector3d faceNormal = (v1 - v0).cross(v2 - v0);
                N += (weight == "angle") ? std::acos((v1 - v0).normalized().dot((v2 - v0).normalized())) * faceNormal.normalized() : faceNormal;
                area += 0.5 * faceNormal.norm();
            }
        }

        area /= 6.0; // account for double-counting
        double d = N.dot(-Hn);
        auto sign = (d > 0) ? 1.0 : -1.0;
        if (std::abs(d) < 1e-12) sign = 0;
        meanCurvature[i] = 2 * sign * Hn.norm() / area;
    }

    std::cout << "mean curvature computed for vertices\n";
}

void Mesh::color_computing() {
    std::vector<float> tmp = meanCurvature;
    std::sort(tmp.begin(), tmp.end());
    float p2 = tmp[size_t(0.05f * tmp.size())];
    float p98 = tmp[size_t(0.95f * tmp.size())];
    float range = p98 - p2;
    if (range < 1e-8f) range = 1.0f;

    #pragma omp parallel for
    for (size_t i = 0; i < meanCurvature.size(); ++i) {
        float t = (meanCurvature[i] - p2) / range;
        t = std::clamp(t, 0.0f, 1.0f);
        t = 1.0f - t;
        Eigen::Vector3d color;
        // use the same Jet-like mapping as before
        if (t <= 0.25f) {
            float f = t / 0.25f;
            color << 0.0, f, 1.0;              // blue → cyan
        } else if (t <= 0.5f) {
            float f = (t - 0.25f) / 0.25f;
            color << 0.0, 1.0, 1.0 - f;        // cyan → green
        } else if (t <= 0.75f) {
            float f = (t - 0.5f) / 0.25f;
            color << f, 1.0, 0.0;              // green → yellow
        } else {
            float f = (t - 0.75f) / 0.25f;
            color << 1.0, 1.0 - f, 0.0;        // yellow → red
        }

        vertices[i].color = color;
    }
    std::cout << "colors computed for vertices\n";
}

void Mesh::color_computing_signed() {
    if (meanCurvature.empty()) return;

    // Robust min/max using 5th and 95th percentiles
    std::vector<float> tmp = meanCurvature;
    std::sort(tmp.begin(), tmp.end());
    float Hmin = tmp[size_t(0.10f * tmp.size())];
    float Hmax = tmp[size_t(0.90f * tmp.size())];
    float Habs = std::max(std::abs(Hmin), std::abs(Hmax));
    if (Habs < 1e-8f) Habs = 1.0f;

    #pragma omp parallel for
    for (size_t i = 0; i < meanCurvature.size(); ++i) {
        // Normalize in [-1, 1]
        float t = std::clamp(-meanCurvature[i] / Habs, -1.0f, 1.0f);

        Eigen::Vector3d color;

        if (t <= -0.5f) {
            // blue -> cyan
            float f = (t + 1.0f) / 0.5f; // 0 -> 1
            color << 0.0, f, 1.0;        // blue->cyan
        } else if (t < 0.0f) {
            // cyan -> green
            float f = (t + 0.5f) / 0.5f; // 0 -> 1
            color << 0.0, 1.0, 1.0 - f;  // cyan->green
        } else if (t < 0.5f) {
            // green -> yellow
            float f = t / 0.5f;          // 0 -> 1
            color << f, 1.0, 0.0;        // green->yellow
        } else {
            // yellow -> red
            float f = (t - 0.5f) / 0.5f; // 0 -> 1
            color << 1.0, 1.0 - f, 0.0;  // yellow->red
        }

        // Convert to [0,255] for PLY
        vertices[i].color = color;
    }

    std::cout << "Colors computed for signed curvature.\n";
}



void Mesh::exportOBJ(const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Cannot open file " << filename << " for writing.\n";
        return;
    }
	std::cout << "Exporting mesh to " << filename << "...\n";
    // Write vertices
    for (const auto& v : vertices) {
        if (!v.removed)
            out << "v " << v.pos.x() << " " << v.pos.y() << " " << v.pos.z() << " " << v.color.x() << " " << v.color.y() << " " << v.color.z() << "\n";
    }
    std::cout << "Exporting mesh to " << filename << "...\n";
    // Mapping from old indices to new indices
    std::vector<int> indexMap(vertices.size(), -1);
    int idx = 1;
    for (size_t i = 0; i < vertices.size(); ++i) {
        if (!vertices[i].removed)
            indexMap[i] = idx++;
    }
    std::cout << "Exporting mesh to " << filename << "...\n";
    // Write faces
    for (const auto& f : faces) {
        if (f.removed) continue;
        if (vertices[f.v[0]].removed ||
            vertices[f.v[1]].removed ||
            vertices[f.v[2]].removed)
            continue;
        out << "f "
            << indexMap[f.v[0]] << " "
            << indexMap[f.v[1]] << " "
            << indexMap[f.v[2]] << "\n";
    }
    std::cout << "Exporting mesh to " << filename << "...\n";
    out.close();
    std::cout << "Mesh exported to " << filename << "\n";
}

#include <fstream>
#include <iostream>
#include <iomanip>

void Mesh::exportPLY(const std::string& filename) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        return;
    }

    // Header
    ofs << "ply\n";
    ofs << "format ascii 1.0\n";
    ofs << "element vertex " << vertices.size() << "\n";
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "property uchar red\n";
    ofs << "property uchar green\n";
    ofs << "property uchar blue\n";
    ofs << "property float quality\n"; // custom scalar
    ofs << "element face " << faces.size() << "\n";
    ofs << "property list uchar int vertex_indices\n";
    ofs << "end_header\n";

    // Vertices
    for (size_t i = 0; i < vertices.size(); ++i) {
        const auto& v = vertices[i];
        unsigned char r = static_cast<unsigned char>(std::clamp(v.color.x(), 0.0, 1.0) * 255.0);
        unsigned char g = static_cast<unsigned char>(std::clamp(v.color.y(), 0.0, 1.0) * 255.0);
        unsigned char b = static_cast<unsigned char>(std::clamp(v.color.z(), 0.0, 1.0) * 255.0);

        ofs << std::fixed << std::setprecision(6)
            << v.pos.x() << " " << v.pos.y() << " " << v.pos.z() << " "
            << (int)r << " " << (int)g << " " << (int)b << " "
            << meanCurvature[i] << "\n"; // optional quality
    }

    // Faces
    for (const auto& f : faces) {
        ofs << "3 " << f.v[0] << " " << f.v[1] << " " << f.v[2] << "\n";
    }

    ofs.close();
    std::cout << "Exported PLY: " << filename << std::endl;
}
