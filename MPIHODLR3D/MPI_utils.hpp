#ifndef _MPI_utils_HPP
#define _MPI_utils_HPP

#include <mpi.h>

void inline set_parallel_work(const size_t N, const int numproc, const int myid, size_t &strt_id, size_t &chunk)
{
    int rem = N % numproc;
    strt_id = 0;
    if (myid + 1 <= rem)
    {
        chunk = ceil(float(N) / float(numproc));
        strt_id = chunk * myid;
    }
    else
    {
        chunk = floor(float(N) / float(numproc));
        strt_id = rem + chunk * myid;
    }
}

void set_for_gatherv(int *&rcv, int *&dspls, int veclen, int nProcs)
{
    int Q, R;
    rcv = new int[nProcs];
    dspls = new int[nProcs];
    Q = veclen / nProcs;
    R = veclen % nProcs;
    for (int i = 0; i < nProcs; i++)
    {
        if (R != 0)
        {
            rcv[i] = Q + 1;
            R--;
        }
        else
            rcv[i] = Q;
    }
    dspls[0] = 0;
    for (int i = 1; i < nProcs; i++)
        dspls[i] = dspls[i - 1] + rcv[i - 1];
}

double getAverage(double var, MPI_Comm MPI_GROUP = MPI_COMM_WORLD)
{
    double AVG = 0.0;
    int np;
    MPI_Allreduce(&var, &AVG, 1, MPI_DOUBLE, MPI_SUM, MPI_GROUP);
    MPI_Comm_size(MPI_GROUP, &np);
    AVG /= np;
    return AVG;
}

double getMax(double var, MPI_Comm MPI_GROUP = MPI_COMM_WORLD)
{
    double MAX;
    MPI_Allreduce(&var, &MAX, 1, MPI_DOUBLE, MPI_MAX, MPI_GROUP);
    return MAX;
}

#endif