#include "boundaryConditions.hpp"
        
BoundaryConditions::BoundaryConditions(std::array<BCType,3>types,std::array<double,3>values)
           :types_(types),values_(values){};

void BoundaryConditions::applyBCsToStencil(Grid3D& grid,double dx, double cond) const{
    const std::vector<std::array<std::size_t,3>>& boundaryIdxs = grid.boundaryIndices();

    for (const std::array<std::size_t,3>& idx:boundaryIdxs){
        auto [i,j,k] = idx;
        int faceNum = static_cast<int>(grid.faceType(i,j,k))-1;
        const std::vector<NeighbourType> solidNeighbors = grid.getSolidNeighbours(i, j, k);
        int numSolidNeigbours = solidNeighbors.size();

        if(numSolidNeigbours==0)continue;
        double weightBc = 1.0/numSolidNeigbours;
        
        grid(i,j,k) = 0.0; //reset so we dont accumulate previous values
        
        for (NeighbourType neighbour :solidNeighbors)
        {
            if (types_[faceNum]==BCType::Dirichlet)
            {
                grid(i,j,k)+=weightBc*values_[faceNum];
            }
                
            else if (types_[faceNum]==BCType::Neumann)
            {
                int  s  = static_cast<int>(neighbour);
                auto ic = i + interiorOffsets[s][0];
                auto jc = j + interiorOffsets[s][1];
                auto kc = k + interiorOffsets[s][2];
                grid(i,j,k)+= weightBc*((2*dx*values_[faceNum]/cond)+grid(ic,jc,kc));
            }
        }
    }
 }

void BoundaryConditions::applyBCsToRhsMatrix(const Grid3D& grid, std::size_t nx, std::size_t ny, double dx,
                                                 double coeff, double cond, std::vector<double>& b) const{

    auto applyBC = [&](std::size_t faceInx, std::size_t row){
                    if(types_[faceInx]==BCType::Dirichlet)b[row] += 2.0*coeff*values_[faceInx];
                    else if(types_[faceInx]==BCType::Neumann)b[row] += (2.0*coeff*dx/cond)*values_[faceInx]; };  
    
     
    const std::vector<std::size_t>& compactLookup = grid.compactLookup();

    for (const std::array<std::size_t,3>& cell:grid.boundaryIndices()){
        auto [i,j,k] = cell;
        std::size_t idx = i + nx*(j+ny*k);
        std::size_t row =  compactLookup[idx]; 
        const std::vector<NeighbourType> solidNeighbors = grid.findSolidNeighbours(i, j, k);       
        int faceNum = static_cast<int>(grid.faceType(i,j,k))-1;

        for (NeighbourType neighbour :solidNeighbors) applyBC(faceNum, row);
    }
}
