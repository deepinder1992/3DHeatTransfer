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
#include <CLI/CLI.hpp>

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
    std::cout << "Running: "  << solver.name() << " Solver!"<< std::endl;
    
    BinaryWriter binWriter("../BinaryOutput", "temperature");
    VTKWriter vtkWriter("../VTKOutput", "temperature");

    double maxDiff = 0.0;                      

    auto start = std::chrono::high_resolution_clock::now();

    for (globs.t = 0; globs.t < globs.steps; ++globs.t) {
        solver.step(current, next, globs, bc);

        if (globs.t % globs.writeInterval == 0) {
            binWriter.write(next, globs.t);
            vtkWriter.write(next, globs.t);
        }

        std::swap(current, next);
        std::size_t imax = 0, jmax = 0, kmax = 0;
        maxDiff = 0.0;
        for (auto& cell:current.activeIndices()) {
            auto [i,j,k] = cell;
            double diff = std::abs(current(i,j,k) - next(i,j,k));
            if (diff > maxDiff){
                imax = i;
                jmax = j;
                kmax = k;
                maxDiff = diff;}
        }

        std::cout << "Step " << globs.t + 1
                  << ": center " << current(globs.nx/2, globs.ny/2, globs.nz/2)
                  << ", Max Diff: " << maxDiff << " at cell: " <<imax<<" "<<jmax<<" "<<kmax<< std::endl;

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
        //load stl file
    VoxelReader(globs.stlFileloc, current);
    current.diagnostics();
    current.assignNoneCells();
    current.constructNeigbourMap(globs.solver);
    
    current.fill(75.0);
    // auto deep copy
    Grid3D next = current;

    

    BoundaryConditions bc(globs.types, globs.values);
    //bc.applyBCsToStencil(current, current.dx(), globs.k);

    LinearAlgebra linAlgebra(globs.maxIters);

    // Select solver based on CLI input
    switch (globs.solver) {
        case SolverType::CPU_STENCIL: {
            HeatSolverCPUStencil solver(globs.alpha, current.dx(), globs.dt, linAlgebra);
            runSimulation(solver, current, next, globs, bc);
            break;
        }
        case SolverType::CPU_MATRIX: {
            HeatSolverCPUMatrix solver(current, nx, ny, nz, globs.alpha, current.dx(), globs.dt, globs.k, bc, linAlgebra);
            runSimulation(solver, current, next, globs, bc);
            break;
        }
        case SolverType::CUDA_STENCIL: {
            HeatSolverCUDAStencil solver(globs.alpha, current.dx(), globs.dt, linAlgebra);
            runSimulation(solver, current, next, globs, bc);
            break;
        }
        case SolverType::CUDA_MATRIX: {
            HeatSolverCUDAMatrix solver(current, nx, ny, nz, globs.alpha, current.dx(), globs.dt, globs.k, bc, linAlgebra);
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
void parseCLI(int argc, char** argv, SimulationGlobals& globs) {
    CLI::App app{"Heat3D Solver"};
    
    app.add_option("--solver", globs.solver, "Select solver type:\n"
                                        "  1 = CPU_STENCIL\n"
                                        "  2 = CPU_MATRIX\n"
                                        "  3 = CUDA_STENCIL\n"
                                        "  4 = CUDA_MATRIX")->check(CLI::Range(1,4));
                                        
    app.add_option("--nx", globs.nx, "Grid size (uniform in x, y, z).");

    app.add_option("--steps", globs.steps, "Maximum number of time steps.");

    app.add_option("--dt", globs.dt, "Time step size (s).");

    app.add_option("--jacobiTol", globs.tol, "Tolerance for Jacobi iterations.");

    app.add_option("--globalTol", globs.globalTol, "Global convergence tolerance for the simulation.");

    app.add_option("--verbosity", globs.verbosity, "Output verbosity level (bitmask):\n" 
                                                    "1 = low, 2 = medium, 4 = high (can combine e.g., 3 = low+medium)\n");

    app.add_option("--writeInterval", globs.writeInterval, "Write output every N steps.");

    app.add_option("--blockDim", globs.blockDim, "CUDA block size (for GPU solvers).");
    

    app.add_option("--bcTypeInlet", globs.types[0], "BC type at Inlet face: \n"
                                                    "  1 = Dirichlet \n"
                                                    "  2 = Neumann")->check(CLI::Range(1,2));                

    app.add_option("--bcTypeOutlet", globs.types[1], "BC type at Outlet face: \n"
                                                    "  1 = Dirichlet \n"
                                                    "  2 = Neumann")->check(CLI::Range(1,2));
        
    app.add_option("--bcTypeWall", globs.types[2], "BC type at Wall face: \n"
                                                    "  1 = Dirichlet \n"
                                                    "  2 = Neumann")->check(CLI::Range(1,2));
    
    app.add_option("--bcValInlet", globs.values[0], "BC value at Inlet face.");              

    app.add_option("--bcValOutlet", globs.values[1], "BC value at Outelet face");
    
    app.add_option("--bcValWall", globs.values[2], "BC value at Wall face");  

    app.add_option("--stlPath", globs.stlFileloc, "Location of the stl file, just provide the stl geometry file name \n"
                                                  "make sure the inlet,outlet wall files are present in samelocation with \n"
                                                  "_inlet, _outlet,_wall appedned to the name of geometry file. \n"
                                                  "eg. cylinder.stl, cylinder_inlet.stl, cylinder_outlet.stl, cylinder_wall.stl \n");

    try {
        app.parse(argc, argv);}
    catch (const CLI::ParseError &e) {
        std::exit(app.exit(e)); }

    globs.solver = static_cast<SolverType>(globs.solver);

    // enforce uniform grid
    globs.ny = globs.nx;
    globs.nz = globs.nx;
    for (auto& t : globs.types) {
    t = static_cast<BCType>(static_cast<int>(t)-1);}

}