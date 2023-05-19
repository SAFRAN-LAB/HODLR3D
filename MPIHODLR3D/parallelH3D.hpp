#ifndef _parallelH3D_HPP
#define _parallelH3D_HPP

#include "FMM3DTreeRAMeff2.hpp"
#include <cmath>
#include <mpi.h>

// Function to check if x is power of 2
bool inline isPowerOfTwo(int n){
    if (n == 0)
        return false;
    return (ceil(log2(n)) == floor(log2(n)));
}
/*
This header files create a parallel construct for the
HODLR3D low-rank structure.
*/

int logn(double x,int n=10){
    return ceil(log(x)/double(log(n)));
}

class parallelH3Dtree{
    FMM3DTree* K;
    int nLevels;
    int np; // Number of MPI process (Expected to be power 2)
    size_t N;
    int TARGET_LEVEL=1;
    int numprocs,myid;
    double Init_time = 0.0;
    double timeHODLR3DRepres = 0.0, timeMatVecProduct = 0.0, timematVecComm = 0.0;
    public:
    // Constructor 
        parallelH3Dtree(int np,kernel*& mykernel, int cubeRootN, int nLevels, int nParticlesInLeafAlong1D, double L, int TOL_POW){
            
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD,&myid);
        
        // Update process information to the TREE class (possible update would be to inherit)
        this->np = np;
        if (numprocs > 1)
            TARGET_LEVEL = ceil(logn(np,8));
        if(myid == 0){
            std::cout << "MPI Process Information set to tree" << std::endl;
            std::cout << "Target Level" << TARGET_LEVEL << std::endl;
        }
        K = new FMM3DTree(mykernel, cubeRootN, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);
        this->N = K->N;
        K->set_process_info(numprocs, myid, TARGET_LEVEL);
        K->set_Uniform_Nodes();
        K->createTree();
        K->assign_Tree_Interactions();
        K->assign_Center_Location();
        K->assignChargeLocations();
        K->assignNonLeafChargeLocations();
        this->nLevels = nLevels;
        schedule_process();   
        }
        void schedule_process();
        void Initialize_parallelH3D(){
            double start	=	MPI_Wtime();
            K->getUVtTree();
            Init_time	=	MPI_Wtime() - start;
            MPI_Barrier(MPI_COMM_WORLD);
            double avg = getAverage(Init_time, MPI_COMM_WORLD);
            double mx = getMax(Init_time, MPI_COMM_WORLD);
            if(myid == 0){
                std::cout << " ++++ Time to find Low-rank basis ++++ " << std::endl;
                std::cout << "(Avg,Max) = " << avg << "," << mx << std::endl; 
            }
        }
        Eigen::VectorXd mat_vec(Eigen::VectorXd& x){
                Eigen::VectorXd bt,b;
                b = Eigen::VectorXd::Zero(N);
                K->assignCharges(x);
                K->evaluateFarField();
                K->evaluate_NearField();
                if (numprocs != 1 && TARGET_LEVEL > 1)
                {
                    for (int i = 1; i < TARGET_LEVEL; i++)
                    {
                        MPI_Barrier(MPI_COMM_WORLD);
                        // Create a new communicator group
                        MPI_Comm MPI_COMM_GROUP;
                        int groupid, groupsize;
                        groupsize = np / pow(8, i);
                        int groupcolor = myid / groupsize;
                        MPI_Comm_split(MPI_COMM_WORLD, groupcolor, myid, &MPI_COMM_GROUP);
                        MPI_Comm_size(MPI_COMM_GROUP, &groupsize);
                        MPI_Comm_rank(MPI_COMM_GROUP, &groupid);
                        std::cout << groupcolor << " "
                                  << "Group RANK/SIZE:" << groupid << "/" << groupsize << std::endl;
                        MPI_Barrier(MPI_COMM_WORLD);
                        K->evaluateFarFieldSuper(i, groupcolor, MPI_COMM_GROUP);
                        MPI_Barrier(MPI_COMM_GROUP);
                        // Free Communicator
                        MPI_Comm_free(&MPI_COMM_GROUP);
                        MPI_Barrier(MPI_COMM_WORLD);
                    }
                }
                K->collectPotential(b);
                bt = Eigen::VectorXd::Zero(N);
                // Measure the Communication time
                double start	=	MPI_Wtime();
                MPI_Allreduce(b.data(),bt.data(),N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                timematVecComm	=	MPI_Wtime() - start;
                timematVecComm += K->timematVecComm;
                K->reorder(bt);
                double avg = getAverage(timematVecComm, MPI_COMM_WORLD);
                double mx = getMax(timematVecComm, MPI_COMM_WORLD);
                if (myid == 0)
                {
                    std::cout << " ++++ Time to Communicate among process ++++ " << std::endl;
                    std::cout << "(Avg,Max) = " << avg << "," << mx << std::endl;
                }
                timeHODLR3DRepres += K->assTime;
                avg = getAverage(timeHODLR3DRepres, MPI_COMM_WORLD);
                mx = getMax(timeHODLR3DRepres, MPI_COMM_WORLD);
                if (myid == 0)
                {
                    std::cout << " ++++ Time to generate entries for HODLR3D ++++ " << std::endl;
                    std::cout << "(Avg,Max) = " << avg << "," << mx << std::endl;
                }
                timeMatVecProduct = K->matVecTime;
                avg = getAverage(timeMatVecProduct, MPI_COMM_WORLD);
                mx = getMax(timeMatVecProduct, MPI_COMM_WORLD);
                if (myid == 0)
                {
                    std::cout << " ++++ Time to matrix-vector product ++++ " << std::endl;
                    std::cout << "(Avg,Max) = " << avg << "," << mx << std::endl;
                }
                return bt;
        }
        size_t system_size(){
            return N;
        }
        ~parallelH3Dtree(){
            delete K;
        }
};
void parallelH3Dtree::schedule_process(){
    size_t strt_node,en_node;
    size_t Chunk = 1;
    //std::cout << " Scheduling ..."<< nLevels << std::endl;
    if(TARGET_LEVEL > 0){
        // This routine sets supernode processor id 
        for (int l = 1; l < TARGET_LEVEL; l++){
            int boxid = myid * pow(8, l) / np;
            K->tree[l][boxid].my_mpi_process = myid;

        }
        Chunk = pow(8, TARGET_LEVEL) / np;
        for (int l = TARGET_LEVEL; l <= nLevels; ++l)
        {
            strt_node = Chunk * myid;
            en_node = strt_node + Chunk;
            //std::cout << "Schedule--" << strt_node << " ," << en_node << " ," << myid << std::endl;
            for (size_t j = strt_node; j < en_node; ++j)
                K->tree[l][j].my_mpi_process = myid;
            Chunk *= 8;
        }       
    }
    else{// Without any parallelization
        for (int l = 1; l <= nLevels; ++l)
            for (size_t k = 0; k < K->nBoxesPerLevel[l]; k++)
                K->tree[l][k].my_mpi_process = 0;
    }
    std::cout << " Scheduled ..." << std::endl;
}
#endif