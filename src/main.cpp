#include "mesh.hpp"
#include <iostream>
#include "stl_reader.hpp"
#include <chrono>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " meshes/input.obj meshes/output.obj \n";
        return 1;
    }
    std::string method = "";
    if (argc == 4) {
        method = argv[3];
    }
    Mesh mesh;
	auto filename = std::string(argv[1]);
    auto pos = filename.find_last_of('.');
 
    auto ext = filename.substr(pos + 1);

    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == "obj")
        mesh.loadOBJ(filename);
    else if (ext == "stl")
        io::loadSTL(filename, mesh);
    else if (ext == "off")
        mesh.loadOFF(filename);
    else
        std::cerr << "Unsupported file type: " << ext << "\n";

    mesh.VVAdj2();
    mesh.EFAdj();
    auto start = std::chrono::high_resolution_clock::now();
    if(method == "uniform") 
        mesh.curvature_computing_uniform();
    else
        mesh.curvature_computing(method);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_ms = end - start;
    std::cout << "Elapsed time for curvature computation: " << duration_ms.count() << " ms\n";
    mesh.color_computing_signed();
	mesh.exportPLY(argv[2]);
    
    return 0;
}

