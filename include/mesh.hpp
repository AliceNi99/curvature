#pragma once

#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <set>
#include <vector>
#include <array>
#include <Eigen/Dense>


struct Vertex {
    Eigen::Vector3d pos;
    Eigen::Vector3d color;
    bool removed = false;   
};

struct Edge {
    int v0, v1;                 
    std::vector<int> faces;     
};

struct Face {
    std::array<int, 3> v;
    bool removed = false;
};

struct UndirectedEdgeHash {
    std::size_t operator()(const std::pair<int,int>& e) const noexcept {
        int a = std::min(e.first, e.second);
        int b = std::max(e.first, e.second);
        // 64-bit mix using large primes
        return ((std::size_t)a * 73856093u) ^ ((std::size_t)b * 19349663u);
    }
};

struct Mesh {
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::vector<std::unordered_set<int>> vertexFaces;
    std::vector<Edge> edges;     // all unique edges
    std::vector<std::vector<int>> VV;
    std::vector<std::vector<int>> VF;
    std::unordered_map<std::pair<int,int>, std::vector<int>, UndirectedEdgeHash> EF;
    std::vector<float> meanCurvature;

    bool loadOFF(const std::string & path);
    bool loadOBJ(const std::string& path);
    void VVAdj();
    void VVAdj2();
    void VFAdj();
    void EFAdj();
    void curvature_computing_uniform();
    void curvature_computing(const std::string& weight = "area");
    void color_computing();
    void color_computing_signed();
    void exportOBJ(const std::string& filename);
    void exportPLY(const std::string& filename);
};

