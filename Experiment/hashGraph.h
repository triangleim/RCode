#ifndef HASH_GRAPH_H
#define HASH_GRAPH_H
#include "Graph.h"

class hashGraph : public Graph
{
public:
    hashGraph(int *hashes, const int N, int n_edge)
        : ids(new int[N]), hashes(new int[N]), hashesMin(new int[N]), hashesMax(new int[N]), Graph(N), M(n_edge)
    {
        utils::sort_indexes(hashes, this->hashes, this->ids, N);
    }
    hashGraph(const int N, const string hfilename, const string idfilename, int n_edge)
        : ids(new int[N]), hashes(new int[N]), hashesMin(new int[N]), hashesMax(new int[N]), Graph(N), M(n_edge)
    {
        utils::readCSV(hfilename, hashes);
        utils::readId(idfilename, ids);
    }
    hashGraph(const vec *OG, int *hashes, const int N, int n_edge)
        : hashGraph(hashes, N, n_edge)
    {
        processHG(OG);
    }
    hashGraph(const vec *OG, const int N, const string hfilename, const string idfilename, int n_edge)
        : hashGraph(N, hfilename, idfilename, n_edge)
    {
        processHG(OG);
    }
    hashGraph(const vec *OG, const int N, const string grofilename, int n_edge)
        : ids(new int[N]), hashes(new int[N]), hashesMin(new int[N]), hashesMax(new int[N]), Graph(N), M(n_edge)
    {
        int *currHashes = new int[N];
        utils::readCSV(grofilename, currHashes, false, ' ');
        utils::sort_indexes(currHashes, this->hashes, this->ids, N);
        delete[] currHashes;
        processHG(OG);
    }
    ~hashGraph()
    {
        delete[] G;
        delete[] ids;
        delete[] hashes;
        delete[] hashesMin;
        delete[] hashesMax;
    }
    int getId(const int i) const
    {
        return ids[i];
    }
    void reportRatio(ostream &);
    double calTri(bool reportNumber = false);

    virtual void commonNeighbor(const int i, const int n, unsigned long &result);
    virtual void commonNeighborNB(const int i, const int n, unsigned long &result);
    virtual void commonNeighborBS(const int i, const int j, unsigned long &cn);
    virtual void commonNeighborNBBS(const int i, const int j, unsigned long &cn);
    virtual void commonNeighborCP(const int i, const int n, unsigned long &result);
    virtual void commonNeighborNBCP(const int i, const int n, unsigned long &result);

    virtual double calCN(const pvec &nodePairs, const bool NB);
    virtual double calCNBS(const pvec &nodePairs, const bool NB);
    virtual unsigned long calCNCP(const pvec &nodePairs, const bool NB);
    virtual double calF(const pvec &nodePairs, const bool NB);
    virtual double calR(const pvec &nodePairs);

    virtual double MC(const vec &nodes, const int *Vrank);

    void writeHashId(const string hfilename, const string idfilename);

    int *oa;
    int *ca;

    int M;

    int *ids;
    int *hashes;
    int *hashesMin;
    int *hashesMax;
#ifndef CSR
    void commonNeighborNB(const vec &vec1, const vec &vec2, unsigned long &result, const int l1, const int l2, const int m, const int M);
    void commonNeighborNB(const vec &vec1, const vec &vec2, vec &result, const int l1, const int l2, const int m, const int M);
#else
    void commonNeighborNB(const int *vec1, int *vec2, unsigned long &result, const int l1, const int l2, const int m, const int M);
    void commonNeighborNB(const int *vec1, int *vec2, vec &result, const int l1, const int l2, const int m, const int M);
#endif
protected:
    // #ifndef CSR
    int selectPivot(const vec &P, const vec &X, const int mp, const int Mp);
    virtual void BKP(vec &P, vec &X, unsigned long &result);
    // #else
    //     int selectPivot(int* P, int* X, int ps,int xs,const int mp, const int Mp);
    //     virtual void BKP(int* P, int* X, int ps,int xs,unsigned long &result);
    // #endif
private:
    void processHG(const vec *OG);
};

#endif