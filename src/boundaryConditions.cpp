#include "boundaryConditions.hpp"

static const std::array<int,3> interiorOffsets[6] = {{+2,0,0}, {-3,0,0}, {0,+2,0},
                                                        {0,-3,0}, {0,0,+2}, {0,0,-3}};
        

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
    const size_type nx = grid.nx();
    const size_type ny = grid.ny();
    const size_type nz = grid.nz();
    const std::vector<std::array<size_type,3>>& boundaryIdxs = grid.boundaryIndices();

    for (const std::array<size_type,3>& idx:boundaryIdxs){
        auto [i,j,k] = idx;
        int faceNum = static_cast<int>(grid.faceType(i,j,k))-1;
        const std::vector<NeighbourType> solidNeighbors = grid.findSolidNeigbour(i, j, k);
        int numSolidNeigbours = solidNeighbors.size();

        if(numSolidNeigbours==0){
            std::cout<<"Cell centered at "<<i<<", "<<j<<", "<<k
            <<" is boundary but no solid neigbours found!"<<std::endl;
            continue;}

        float weightBc = 1.0/numSolidNeigbours;

        // std::array<size_t,3> interiorOffsets[6] = {{i+2,j,k}, {i-3,j,k}, {i,j+2,k},
        //                                     {i,j-3,k}, {i,j,k+2}, {i,j,k-3}};
        
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
                grid(i,j,k)+= (weightBc*sign_*2*dx*values_[faceNum]/cond)+grid(ic,jc,kc);
            }
        }
    }
 }

void BoundaryConditions::applyBCsToRhsMatrix(const Grid3D& grid, size_type nx, size_type ny, size_type nz, double dx,
                                                 double coeff, double cond, std::vector<double>& b) const{

    auto applyBC = [&](size_type faceInx, size_type row, int sign){
                    if(types_[faceInx]==BCType::Dirichlet)b[row] += 2.0*coeff*values_[faceInx];
                    else if(types_[faceInx]==BCType::Neumann)b[row] += (sign*2.0*coeff*dx/cond)*values_[faceInx]; };  
    
    const std::vector<std::array<size_type,3>>& boundaryIdxs = grid.boundaryIndices();
    const std::vector<std::size_t>& compactLookup = grid.compactLookup();

    for (const std::array<size_type,3>& cell:boundaryIdxs){
        auto [i,j,k] = cell;
        size_type idx = i + nx*(j+ny*k);
        size_type row =  compactLookup[idx];

        const std::vector<NeighbourType> solidNeighbors = grid.findSolidNeigbour(i, j, k);
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

void BoundaryConditions::applyBCsToStencilCUDA(double* grid, double dx, size_type nx, 
    size_type ny, size_type nz, double cond, dim3 gridCuda, dim3 blockCuda) const
    {    
        //hard coded num 6 works for cubes and cuboids, however more elegant handling 
        //is needed once we generalize to arbitrary shapes      
        int numFaces = 6;
        BCType types[numFaces];
        for (int i = 0; i < numFaces; ++i)
            types[i] = types_[i];

        double values[numFaces];
        for (int i = 0; i < numFaces; ++i)
            values[i] = values_[i];

        applyBCsToStencilKern<<<gridCuda, blockCuda>>>(grid, nx, ny, nz, dx, cond, types, values);
        cudaDeviceSynchronize();
    }