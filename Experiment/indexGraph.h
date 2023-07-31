#ifndef INDEX_GRAPH_H
#define INDEX_GRAPH_H
#include "hashGraph.h"
// #define CNNBV 2
class indexGraph : public hashGraph
{
public:
    indexGraph(const vec *OG, int *hashes, const int N, const int M)
        : hashEntries(new vec[N]), hashIds(new vec[N]), hashGraph(hashes, N, M)
    {
        processHGIndex(OG);
    }
    indexGraph(const vec *OG, const int N, const string hfilename, const string idfilename, const int M, bool ismc = false)
        : hashEntries(new vec[N]), hashIds(new vec[N]), hashGraph(N, hfilename, idfilename, M)
    {
        if (ismc)
        {
            processHGIndexMC(OG);
        }
        else
            processHGIndex(OG);
    }
    ~indexGraph()
    {
        delete[] hashEntries;
        delete[] hashIds;
    }
    // double calTri(bool reportNumber=false);
    virtual void commonNeighbor(const int i, const int n, unsigned long &result);
    virtual void commonNeighborNB(const int i, const int n, unsigned long &result);
    void reportRatio(ostream &out = cout);
    virtual double calF(const pvec &nodePairs, const bool NB);
    double calF2(const pvec &nodePairs, const bool NB);
    int *hashEntriesCA;
    int *hashIdsCA;
    int *hashOA;
    int *degreesh;
    void calRate();
    double MC(const vec &nodes, const int *Vrank);
    void BKP(vec &P, vec &X, unsigned long &result);
    void BKP(vec &PN, vec &XN, vec &hashEntsP, vec &hashEntsX, vec &hashIdP, vec &hashIdX, unsigned long &result);
    void BKP2(vec &PN, vec &XN, vec &hashEntsP, vec &hashEntsX, vec &hashIdP, vec &hashIdX, unsigned long &result);
    int selectPivot(const vec &P, const vec &X, const int mp, const int Mp);
    int selectPivot(const vec &PN, const int ps, const vec &XN, const int mp, const int Mp, vec &hashEntsP, vec &hashEntsX, vec &hashIdP, vec &hashIdX);
    void commonNeighborNB(const int *n1, const int d1, const int *n2, const int d2, unsigned long &result, vec &hashEnts1, int *hashEnts2, vec &hashId1, int *hashId2);
    void commonNeighborNB(const int *n1, const int d1, const int *n2, const int d2, vec &result, vec &hashEnts1, int *hashEnts2, vec &hashId1, int *hashId2, vec &hashEntsR, vec &hashIdR);

protected:
    vec *hashEntries;
    vec *hashIds;
    int **hashEntriesA;
    int **hashIdsA;

private:
    void processHGIndex(const vec *OG);
    void processHGIndexMC(const vec *OG);
    pair<int, int> calCP2(const int i, const int adj2, const int m, const int M);
};

#endif