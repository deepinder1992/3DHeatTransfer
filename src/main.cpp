#include "grid.hpp"
#include "solverCPU.hpp"
#include "solverCUDA.hpp"
#include "boundaryConditions.hpp"
#include "outputWriter.hpp"
#include "simGlobals.hpp"
#include "voxelReader.hpp"
#include <chrono>
#include <iostream>
#include <memory>
#include <string>

// Forward declaration for CLI parser
void parseCLI(int argc, char** argv, SimulationGlobals& g);

// Templated simulation runner
template<typename Solver>
void runSimulation(
    Solver& solver,
    Grid3D& current,
    Grid3D& next,
    SimulationGlobals& globs,
    BoundaryConditions& bc)
{
    BinaryWriter binWriter("../BinaryOutput", "temperature");
    VTKWriter vtkWriter("../VTKOutput", "temperature");

    double maxDiff = 0.0;                      
    size_type gridSize = current.size();

    auto start = std::chrono::high_resolution_clock::now();

    for (globs.t = 0; globs.t < globs.steps; ++globs.t) {
        solver.step(current, next, globs, bc);

        if (globs.t % globs.writeInterval == 0) {
            binWriter.write(next, globs.t);
            vtkWriter.write(next, globs.t);
        }

        std::swap(current, next);

        maxDiff = 0.0;
        for (size_type i = 0; i < gridSize; ++i) {
            double diff = std::abs(current.data()[i] - next.data()[i]);
            if (diff > maxDiff)
                maxDiff = diff;
        }

        std::cout << "Step " << globs.t + 1
                  << ": center " << current(globs.nx/2, globs.ny/2, globs.nz/2)
                  << ", Max Diff: " << maxDiff << std::endl;

        if (maxDiff < globs.globalTol) break;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();

    std::cout << solver.name() << " Total simulation time: " << elapsed << " seconds\n";

    if (globs.verbosity & SimulationGlobals::VERB_HIGH) {
        std::cout << solver.name() << " Total Internal Iters: " << globs.totalIters << std::endl;
    }
}

int main(int argc, char** argv) {
    SimulationGlobals globs;

    // Parse CLI arguments
    parseCLI(argc, argv, globs);

    size_type nx = globs.nx;
    size_type ny = globs.ny;
    size_type nz = globs.nz;
    double dx = globs.dx;

    Grid3D current(nx, ny, nz, dx);
    Grid3D next(nx, ny, nz, dx);
        //load stl file
    VoxelReader(globs.stlFileloc, current);

    current.fill(75.0);

    BoundaryConditions bc(globs.types, globs.values);
    bc.applyBCsToStencil(current, globs.dx, globs.k);

    LinearAlgebra linAlgebra(globs.maxIters);

    // Select solver based on CLI input
    switch (globs.solver) {
        case SolverType::CPU_STENCIL: {
            HeatSolverCPUStencil solver(globs.alpha, globs.dx, globs.dt, linAlgebra);
            runSimulation(solver, current, next, globs, bc);
            break;
        }
        case SolverType::CPU_MATRIX: {
            HeatSolverCPUMatrix solver(nx, ny, nz, globs.alpha, globs.dx, globs.dt, globs.k, bc, linAlgebra);
            runSimulation(solver, current, next, globs, bc);
            break;
        }
        case SolverType::CUDA_STENCIL: {
            HeatSolverCUDAStencil solver(globs.alpha, globs.dx, globs.dt, linAlgebra);
            runSimulation(solver, current, next, globs, bc);
            break;
        }
        case SolverType::CUDA_MATRIX: {
            HeatSolverCUDAMatrix solver(nx, ny, nz, globs.alpha, globs.dx, globs.dt, globs.k, bc, linAlgebra);
            runSimulation(solver, current, next, globs, bc);
            break;
        }
        default:
            std::cerr << "Unknown solver selected!" << std::endl;
            return 1;
    }

    std::cout << "Simulation completed in " << globs.t << " iterations!" << std::endl;
    return 0;
}

// ------------------ CLI parser ------------------
void parseCLI(int argc, char** argv, SimulationGlobals& g) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--solver") {
            int s = std::stoi(argv[++i]);
            if (s < 1 || s > 4) {
                std::cerr << "Solver must be between 1 and 4\n";
                exit(1);
            }
            g.solver = static_cast<SolverType>(s);
        }
        else if (arg == "--nx"){
                g.nx = std::stoi(argv[++i]);
                g.ny = g.nx;  // enforce equal dimensions. the solver is capable to handle different
                g.nz = g.nx;    // but we need dx dy and dz separately. This will be enchanced later
            }
        else if (arg == "--ny") g.ny = std::stoi(argv[++i]);
        else if (arg == "--nz") g.nz = std::stoi(argv[++i]);
        else if (arg == "--steps") g.steps = std::stoi(argv[++i]);
        else if (arg == "--dt") g.dt = std::stod(argv[++i]);
        else if (arg == "--jacobiTol") g.tol = std::stod(argv[++i]);
        else if (arg == "--globalTol") g.globalTol = std::stod(argv[++i]);
        else if (arg == "--lx") g.lx = std::stod(argv[++i]);
        else if (arg == "--maxIters") g.maxIters = std::stoi(argv[++i]);
        else if (arg == "--verbosity") g.verbosity = std::stoi(argv[++i]);
        else if (arg == "--writeInterval") g.writeInterval = std::stoi(argv[++i]);
        else if (arg == "--blockx") g.blockDimX = std::stoi(argv[++i]);
        else if (arg == "--blocky") g.blockDimY = std::stoi(argv[++i]);
        else if (arg == "--blockz") g.blockDimZ = std::stoi(argv[++i]);
        else if (arg == "--bcTypes"){g.types[0]  = static_cast<BCType>(std::stoi(argv[++i]));
                                    g.types[1]  = static_cast<BCType>(std::stoi(argv[++i]));
                                    g.types[2]  = static_cast<BCType>(std::stoi(argv[++i]));}
        else if (arg == "--bcVals") {g.values[0]  = std::stoi(argv[++i]);
                                    g.values[1]  = std::stoi(argv[++i]);}
        else if (arg == "--stlFilePath") g.stlFileloc = std::stoi(argv[++i]);
        else if (arg == "--help") {
            std::cout << "Usage: ./heatSolver [options]\n\n";
            std::cout << "Options:\n";
            std::cout << "  --solver [1-4]       Select solver type:\n";
            std::cout << "                       1 = CPU_STENCIL\n";
            std::cout << "                       2 = CPU_MATRIX\n";
            std::cout << "                       3 = CUDA_STENCIL\n";
            std::cout << "                       4 = CUDA_MATRIX\n";
            std::cout << "  --nx N               Number of grid points in x (ny and nz are set equal to nx by default)\n";
            std::cout << "  --ny N               Number of grid points in y (optional, will override default nx)\n";
            std::cout << "  --nz N               Number of grid points in z (optional, will override default nx)\n";
            std::cout << "  --steps N            Maximum number of time steps\n";
            std::cout << "  --dt VALUE           Time step size\n";
            std::cout << "  --jacobiTol VALUE    Tolerance for Jacobi iterations (if using iterative solver)\n";
            std::cout << "  --globalTol VALUE    Global convergence tolerance for the simulation\n";
            std::cout << "  --lx VALUE           Physical domain size in x (length units)\n";
            std::cout << "  --maxIters N         Maximum number of iterations per time step\n";
            std::cout << "  --verbosity LEVEL    Output verbosity level (bitmask):\n";
            std::cout << "                       1 = low, 2 = medium, 4 = high (can combine e.g., 3 = low+medium)\n";
            std::cout << "  --writeInterval N    Write output every N steps\n";
            std::cout << "  --blockx N           CUDA block size in x (for GPU solvers)\n";
            std::cout << "  --blocky N           CUDA block size in y\n";
            std::cout << "  --blockz N           CUDA block size in z\n";
            std::cout << "  --bcTypes            respective types of BC at inlet, wall and outlet patches\n"
                      << "                       1 = Dirichlet \n"
                      << "                       2 = Neumann \n"
                      << "                       eg. --bcTypes 2 1 2 for inlet, wall and outlet respectively \n";
            std::cout << "  --bcVals             respective values of BC at inlet, wall and outlet patches\n"
                      << "                       eg. --bcVals 1000 100 -1000 for inlet, wall and outlet respectively \n";
            std::cout << "  --binaryStlFilePath  Location of the stl file, just provide the stl geometry file name \n"
                      << "                       make sure the inlet,outlet wall files are present in samelocation with \n"
                      << "                       _inlet, _outlet,_wall appedned to the name of geometry file. \n"
                      << "                       eg. cylinder.stl, cylinder_inlet.stl, cylinder_outlet.stl, cylinder_wall.stl \n";
            std::cout << "  --help               Show this help message\n\n";
            std::cout << "Cmd Examples:\n";
            std::cout << "  ./heatSolver --solver 2 --nx 128 --steps 10000 --dt 0.01 --verbosity 3\n";
            std::cout << "  ./heatSolver --solver 4 --nx 256 --writeInterval 100 --blockx 16 --blocky 16 --blockz 16\n";
            exit(0);
        }
    }

    // recompute dependent values
    g.dx = g.lx / g.nx;
}