#include "boundaryConditions.hpp"
        
BoundaryConditions::BoundaryConditions(std::array<BCType,3>types,std::array<double,3>values)
           :types_(types),values_(values){};

int BoundaryConditions::sign(Vector cellNormal, const std::array<size_type,3>&  cellIdxs,
                                 std::size_t i, std::size_t j, std::size_t k) const {
    Vector radial {static_cast<float>(cellIdxs[0]) - static_cast<float>(i),
                    static_cast<float>(cellIdxs[1]) - static_cast<float>(j),
                    static_cast<float>(cellIdxs[2]) - static_cast<float>(k)};
    
    if (radial.dot(cellNormal)>0) return 1;
    return -1; 
};

void BoundaryConditions::applyBCsToStencil(Grid3D& grid,double dx, double cond) const{
    const std::vector<std::array<size_type,3>>& boundaryIdxs = grid.boundaryIndices();

    for (const std::array<size_type,3>& idx:boundaryIdxs){
        auto [i,j,k] = idx;
        int faceNum = static_cast<int>(grid.faceType(i,j,k))-1;
        const std::vector<NeighbourType> solidNeighbors = grid.getSolidNeighbours(i, j, k);
        int numSolidNeigbours = solidNeighbors.size();

        if(numSolidNeigbours==0)continue;
        double weightBc = 1.0/numSolidNeigbours;
        
        Vector norm = grid.cellFaceNormalized(i,j,k);
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
                int sign_ = sign(norm,idx,ic,jc,kc);
                grid(i,j,k)+= weightBc*((sign_*2*dx*values_[faceNum]/cond)+grid(ic,jc,kc));
            }
        }
    }
 }

void BoundaryConditions::applyBCsToRhsMatrix(const Grid3D& grid, size_type nx, size_type ny, size_type nz, double dx,
                                                 double coeff, double cond, std::vector<double>& b) const{

    auto applyBC = [&](size_type faceInx, size_type row, int sign){
                    if(types_[faceInx]==BCType::Dirichlet)b[row] += 2.0*coeff*values_[faceInx];
                    else if(types_[faceInx]==BCType::Neumann)b[row] += (sign*2.0*coeff*dx/cond)*values_[faceInx]; };  
    
     
    const std::vector<std::size_t>& compactLookup = grid.compactLookup();

    for (const std::array<size_type,3>& cell:grid.boundaryIndices()){
        auto [i,j,k] = cell;
        size_type idx = i + nx*(j+ny*k);
        size_type row =  compactLookup[idx];

        const std::vector<NeighbourType> solidNeighbors = grid.findSolidNeighbours(i, j, k);
        
        int faceNum = static_cast<int>(grid.faceType(i,j,k))-1;
        Vector norm = grid.cellFaceNormalized(i,j,k);

        for (NeighbourType neighbour :solidNeighbors)
        {
            int  s  = static_cast<int>(neighbour);
            auto ic = i + interiorOffsets[s][0];
            auto jc = j + interiorOffsets[s][1];
            auto kc = k + interiorOffsets[s][2];
            int sign_ = sign(norm,cell,ic,jc,kc);
            applyBC(faceNum, row, sign_);
        }
    }
}

void BoundaryConditions::applyBCsToStencilCUDA(double* grid, double* oldGrid, double dx, size_type nx, 
    size_type ny, size_type nz, size_type(*bcIndices)[3], FaceType* faceTypes, std::size_t nBcCells,
    NeighbourType* devNbrTypes, std::size_t* devNbrOffset, float (*devCellNormals)[3], double cond, dim3 gridCuda, dim3 blockCuda) const 
    {    
        int numFaces = 3;
        BCType types[numFaces];
        for (int i = 0; i < numFaces; ++i)
            types[i] = types_[i];

        double values[numFaces];
        for (int i = 0; i < numFaces; ++i)
            values[i] = values_[i];

        applyBCsToStencilKern<<<gridCuda, blockCuda>>>(grid, oldGrid, nx, ny, nz, dx, bcIndices, faceTypes,
             nBcCells, devNbrTypes, devNbrOffset, devCellNormals, cond, types, values);
        cudaDeviceSynchronize();
    }