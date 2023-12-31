#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <climits>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <boost/filesystem.hpp>

using namespace std;

typedef vector<int> vec;
typedef vector<string> stringvec;

template <typename T>
inline void sort_indexes(T *&v, int *&ids, int N)
{
    int *idx = new int[N];
    iota(idx, idx + N, 0);
    stable_sort(idx, idx + N, [v](int i1, int i2)
                { return v[i1] < v[i2]; });

    for (int i = 0; i < N; i++)
        ids[idx[i]] = i;
    delete[] idx;
}

void readGraph(string file, vec *&G, int N)
{
    ifstream hash_File(file.c_str());
    string currLine;
    while (getline(hash_File, currLine))
    {
        stringstream linestream(currLine);
        string ids;
        getline(linestream, ids, ' ');
        int idx = stoi(ids);
        string ns;
        getline(linestream, ns, ' ');
        int nid = stoi(ns);
        G[idx].push_back(nid);
        G[nid].push_back(idx);
    }
    for (int i = 0; i < N; i++)
    {
        sort(G[i].begin(), G[i].end());
        vec::iterator uniqEnd = unique(G[i].begin(), G[i].end());
        G[i].erase(uniqEnd, G[i].end());
    }
}

void readCSV(const string file, int *&hashes, bool withHeader = true, char delim = ',')
{
    ifstream hash_File(file.c_str());
    string currLine;
    if (withHeader)
        getline(hash_File, currLine);
    while (getline(hash_File, currLine))
    {
        stringstream linestream(currLine);
        string ids;
        getline(linestream, ids, delim);
        int idx = stoi(ids);
        string hs;
        getline(linestream, hs, delim);
        hashes[idx] = stoi(hs);
    }
}

void readId(string file, int *&ids)
{
    ifstream hash_File(file.c_str());
    string currLine;
    getline(hash_File, currLine);
    int i = 0;
    while (getline(hash_File, currLine))
    {
        stringstream linestream(currLine);
        string id;
        getline(linestream, id, ',');
        ids[i++] = stoi(id);
    }
}

static void read_directory(const string &name, stringvec &v)
{
    boost::filesystem::path p(name);
    boost::filesystem::directory_iterator start(p);
    boost::filesystem::directory_iterator end;
    while (start != end)
    {
        if (!boost::filesystem::is_directory((*start).path()))
            v.push_back((*start).path().string());
        start++;
    }
}

void processHGSimple(vec *G, int *ids, int N, vec *indexG)
{
    for (int i = 0; i < N; i++)
    {
        int newid = ids[i];
        indexG[newid] = vec();
        vec &adjnodes = indexG[newid];
        adjnodes.reserve(G[i].size());
        for (const auto &adj : G[i])
            adjnodes.push_back(ids[adj]);
        sort(adjnodes.begin(), adjnodes.end());
    }
}

void writeG(const vec *G, string outfile, const int N)
{
    ofstream nodefile;
    nodefile.open(outfile + ".g");
    for (int i = 0; i != N; i++)
    {
        for (const auto &n : G[i])
            nodefile << i << " " << n << endl;
    }
}

int main(int argc, char **argv)
{
    string directname = "../data/";
    string outdirect = "../maxnode_hash_data/test_";
    string groDirect = "../gro_data/";
    string name = argv[1];
    int N = stoi(argv[2]);
    int dense = 0;
    if (argc > 3)
        dense = stoi(argv[3]);
    string filename = directname + name + ".txt";
    // string finaldirect = "./data/" + name + "/" + name;
    string finaldirect = "../sm_data/" + name;
    ofstream labelfile;
    labelfile.open("../sm_data/" + name + ".l");
    std::random_device e;
    std::uniform_real_distribution<double> u(0, 1); //随机数分布对象
    // vector<int> tmps(6, 0);
    for (int i = 0; i != N; i++)
    {
        int la;
        if (dense == 1)
            la = int(pow((pow(N / 1000, 2)) * u(e), (double)1 / (2)));
        else if (dense == 2)
            la = int(pow((pow(N / 100, 2)) * u(e), (double)1 / (2)));
        else la = int(pow((pow(N / 5000, 2)) * u(e), (double)1 / (2)));
        // auto la = 0;
        labelfile << i << " " << la << endl; //
        // tmps[la]++;
    }
    // for (int i = 0; i != 6; i++){
    //     cout<<tmps[i]<<endl;
    // }

    vec *G = new vec[N];
    readGraph(filename, G, N);
    writeG(G, finaldirect, N);
    // groDirect += name;
    stringvec groFiles;
    read_directory(groDirect, groFiles);
    string finaldirect2 = "../sm_data/gro/" + name;
    for (const string &f : groFiles)
    {
        int idx = f.find(name);       //
        if (idx != std::string::npos) //不存在。
        {
            cout << f << endl;
            int *hashes = new int[N];
            int *ids = new int[N];
            readCSV(f, hashes, false, ' ');
            sort_indexes(hashes, ids, N);
            vec *IG = new vec[N];
            processHGSimple(G, ids, N, IG);
            size_t b = f.find_last_of('/');
            size_t first = f.find_first_of('_', b);
            size_t second = f.find_last_of('_');
            cout << f.substr(first, second - first) << endl;
            writeG(IG, finaldirect2 + f.substr(first, second - first), N);
            delete[] hashes;
            delete[] ids;
            delete[] IG;
        }
    }
    {
        string hdirect = outdirect + name + "_range_nb_";
        string hid = hdirect + "id.csv";
        int *ids = new int[N];
        readId(hid, ids);
        vec *IG = new vec[N];
        processHGSimple(G, ids, N, IG);
        writeG(IG, finaldirect + "_range", N);
        delete[] ids;
        delete[] IG;
    }
    {
        string idirect = outdirect + name + "_index_nb_";
        string iid = idirect + "id.csv";
        int *ids = new int[N];
        readId(iid, ids);
        vec *IG = new vec[N];
        processHGSimple(G, ids, N, IG);
        writeG(IG, finaldirect + "_index", N);
        delete[] ids;
        delete[] IG;
    }
    delete[] G;
}