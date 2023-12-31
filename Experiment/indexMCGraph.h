#ifndef INDEX_VEC_GRAPH_H
#define INDEX_VEC_GRAPH_H
#include "utils.h"

typedef struct
{
    int first;
    vec *second;
} entry;
// typedef pair<int, vec*> entry;
typedef vector<entry> adjList;

class indexVecGraph
{
public:
    indexVecGraph(const int VN, const int num_edge)
        : M(num_edge), N(VN), ids(new int[VN]), hashes(new int[VN]), hashesMin(new int[VN]), hashesMax(new int[VN]), degrees(new int[VN]), degreesh(new int[VN]), G(new adjList[VN]), GO(new vec[VN]) {}
    indexVecGraph(const vec *OG, int *hashes, const int VN, const int num_edge)
        : indexVecGraph(VN,num_edge)
    {
        utils::sort_indexes(hashes, this->hashes, this->ids, N);
        processHGIndex(OG);
    }
    indexVecGraph(const vec *OG, const int VN, const string hfilename, const string idfilename, const int num_edge,bool ismc=true)
        : indexVecGraph(VN,num_edge)
    {
        utils::readCSV(hfilename, hashes);
        utils::readId(idfilename, ids);
        processHGIndex(OG);
    }
    void reportRatio(ostream &);
    int *oa;
    entry *ca;
    ~indexVecGraph()
    {
        delete[] ids;
        delete[] hashes;
        delete[] hashesMin;
        delete[] hashesMax;
        delete[] degrees;
        delete[] degreesh;
        delete[] oa;
        delete[] ca;
        for (int i = 0; i != N; i++)
        {
            for (entry &adj : G[i])
                delete adj.second;
        }
        delete[] G;
        delete[] GO;
    }
    int getId(const int i) const
    {
        return ids[i];
    }
    // int* hashEntriesCA;
    // int* hashIdsCA;
    // int* hashOA;

    int M;
    double MC(const vec &nodes, const int *Vrank);

private:
    const int N;
    int *ids;
    int *hashes;
    int *hashesMin;
    int *hashesMax;
    int *degrees;
    int *degreesh;
    adjList *G;
    vec *GO;
    void processHGIndex(const vec *OG);
    int selectPivot(const adjList &PN, const int ps, const adjList &XN, const int mp, const int Mp);
    void BKP(adjList &PN, adjList &XN, unsigned long &result);
    int selectPivot(const vec &P, const vec &X, const int mp, const int Mp);
    void BKP(vec &P, vec &X, unsigned long &result);
};

#endif