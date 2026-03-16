#include "boundaryConditions.hpp"

BoundaryConditions::BoundaryConditions(std::array<BCType,3>types,std::array<double,3>values)
           :types_(types),values_(values){};


void BoundaryConditions::applyBCsToStencil(Grid3D& grid,double dx, double cond) const{
    const size_type nx = grid.nx();
    const size_type ny = grid.ny();
    const size_type nz = grid.nz();
    const std::vector<std::array<size_type,3>>& boundaryIdxs = grid.boundaryIndices();
    

    auto applyBC = [&](int face, double& cell, double interior, int sign, float weightBc ){
             if (types_[face]==BCType::Dirichlet){cell+=weightBc*values_[face];}
            else if (types_[face]==BCType::Neumann){ cell+= (weightBc*sign*2*dx*values_[face]/cond)+interior;}};

    for (std::array<size_type,3> idx:boundaryIdxs){
        auto [i,j,k] = idx;
        int faceNum = static_cast<int>(grid.faceType(i,j,k));
        const std::vector<std::array<size_type,3>> solidNeighbors = grid.findSolidNeigbour(i, j, k);
        int numSolidNeigbours = solidNeighbors.size();
        if(numSolidNeigbours==0){
            std::cout<<"Cell centered at"<<i<<", "<<j<<", "<<k<<"is boundary but no solid neigbour found!"<<std::endl;
            continue;}
        float weightBc = 1.0/numSolidNeigbours;
        std::array<size_t,3> neighborOffsets[6] = {{i-1, j, k}, {i+1, j, k}, {i, j-1, k},
                                                   {i, j+1, k}, {i, j, k-1}, {i, j, k+1}};
        
        grid(i,j,k) = 0.0; //reset so we dont accumulate previous values
        for (std::array<size_type,3> neighbor :solidNeighbors){
            //auto [in,jn,kn] = neighbor; 
            
            if(neighbor==neighborOffsets[0]){applyBC(faceNum, grid(i,j,k),grid(i+2,j,k),-1,weightBc);}
            if(neighbor==neighborOffsets[1]){applyBC(faceNum, grid(i,j,k),grid(i-3,j,k),1,weightBc);}
            if(neighbor==neighborOffsets[2]){applyBC(faceNum, grid(i,j,k), grid(i,j+2,k),-1,weightBc);}
            if(neighbor==neighborOffsets[3]){ applyBC(faceNum, grid(i,j,k),grid(i,j-3,k),1,weightBc);}
            if(neighbor==neighborOffsets[4]){applyBC(faceNum, grid(i,j,k),grid(i,j,k+2),-1,weightBc);}
            if(neighbor==neighborOffsets[5]){applyBC(faceNum, grid(i,j,k),grid(i,j,k-3),1,weightBc);}
        }
    }

}

void BoundaryConditions::applyBCsToRhsMatrix(size_type nx,
                                              size_type ny,
                                               size_type nz, 
                                                double dx,
                                                 double coeff,
                                                  double cond,
                                                   std::vector<double>& b) const
{   
 
    auto applyBC = [&](size_type faceInx, size_type row, int factor){
                    if(types_[faceInx]==BCType::Dirichlet)b[row] += 2.0*coeff*values_[faceInx];
                    else if(types_[faceInx]==BCType::Neumann)b[row] += (factor*2.0*coeff*dx/cond)*values_[faceInx]; };

    // for (size_type k :{size_type(0), nz-1})
    // for (size_type j :{size_type(0), ny-1})
    // for (size_type i :{size_type(0), nx-1})
    for (size_type k =0 ; k< nz; ++k)
    for (size_type j =0 ; j< ny; ++j)
    for (size_type i =0 ; i< nx; ++i)
    {
        size_type row = i + nx*(j+ny*k);

        if(i==0) applyBC(0, row, -1);// ----- X- -----
        if(i==nx-1) applyBC(1,row,1); // ----- X+ -----
        if(j==0) applyBC(2,row,-1); // ----- X- -----       
        if(j==ny-1) applyBC(3,row,1);// ----- Y+ -----        
        if(k==0)applyBC(4,row,-1);// ----- Z- -----       
        if(k==nz-1)applyBC(5,row,1);// ----- Z+ -----

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