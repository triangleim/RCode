#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <climits>
#include <assert.h>
#include <time.h>
#include <vector>
#include <queue>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <unordered_set>
#include <numeric>
#include <set>
#include <array>
#include <chrono>
#include <boost/dynamic_bitset.hpp>
#include "../include/tqdm.h"

using namespace std;

typedef vector<int> vec;

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

template <typename T>
void sort_indexes(T *&v, int *&ids, int N)
{
    int *idx = new int[N];
    iota(idx, idx + N, 0);
    stable_sort(idx, idx + N, [v](int i1, int i2)
                { return v[i1] < v[i2]; });

    T *newV = new T[N];
    for (int i = 0; i < N; i++)
    {
        ids[idx[i]] = i;
        newV[i] = v[idx[i]];
    }
    delete[] v;
    v = newV;
}

inline int binarySearch(int _key, const vec &_list, int _left, int _right)
{
    int _mid;
    while (_left < _right)
    {
        _mid = (_right + _left) / 2;
        if (_list[_mid] < _key)
            _left = _mid + 1;
        else
            _right = _mid;
    }
    return _left;
}

inline void commonNeighbor(vec &vec1, vec &vec2, int &cn, int id1, int ed1, int id2, int ed2)
{
    while (id1 != ed1 && id2 != ed2)
    {
        if (vec1[id1] == vec2[id2])
        {
            cn++;
            id1++;
            id2++;
        }
        else if (vec1[id1] < vec2[id2])
            id1++;
        else
            id2++;
    }
}

inline void updateMinMaxNB(vec *&G, int node, int h, unordered_map<int, unsigned> *&hashesNums, queue<int> &candidates, unordered_set<int> &candidatesUniq, int *&hashesMin, int *&hashesMax, int d1_threshold, int d2_threshold)
{
    //更新某个节点赋值label之后，其周围邻居的RangeCode
    // candidates用来遍历邻居加入候选的节点，candidatesUniq用来防止加入重复的节点
    for (const auto &n : G[node])
    {
        if (h < hashesMin[n])
            hashesMin[n] = h;
        if (h > hashesMax[n])
            hashesMax[n] = h;
        if (candidatesUniq.find(n) == candidatesUniq.end())
        {
            candidates.push(n);
            candidatesUniq.insert(n);
        }
    }
    for (const auto &n : G[node])
    {
        //挑选部分节点进行计算
        int d2_count = 0;
        for (const auto &n2 : G[n])
        {
            if (G[node].size() > d1_threshold && d2_count++ > d2_threshold)    break;
            if (h < hashesMin[n2] || h > hashesMax[n2])
            {
                unordered_map<int, unsigned> &hn2 = hashesNums[n2];
                if (hn2.find(h) != hn2.end())
                    hn2[h]++;
                else
                    hn2[h] = 1;
            }
        }
    }
}

void hashMergeDiffNBEx(vec *&G, int *&except, int *&hashes, int *&hashesMin, int *&hashesMax, int maxDegreeNode, int N, int d1_threshold, int d2_threshold)
{
    // hashesNum是T(v,h)->hashesNums[v][h]
    unordered_map<int, unsigned> *hashesNums = new unordered_map<int, unsigned>[N]();
    boost::dynamic_bitset<> adjVec(N); //记录节点n邻居的二进制set
    boost::dynamic_bitset<> marked(N); //记录节点n是否赋值label的二进制set
    marked.set(maxDegreeNode);
    hashes[maxDegreeNode] = 1;
    int upper = 2;
    int lower = 0;
    long long R_count = 0;
    //每个点的lower upper需要检索更新，不再是全局的lower upper
    // 1、点v两跳邻居的范围 最小l 最大r
    // 2、点v一跳邻居的范围 最大l 最小r
    // vector<unsigned> currNeighbors = adjlist[maxDegreeNode];
    queue<int> candidates;             //候选的节点
    unordered_set<int> candidatesUniq; //用来去重已经赋值过的节点
    int markedNum = 1;
    updateMinMaxNB(G, maxDegreeNode, 1, hashesNums, candidates, candidatesUniq, hashesMin, hashesMax, d1_threshold, d2_threshold);
    tqdm bar;
    bar.progress(1, N);
    maxDegreeNode = 1;
    while (markedNum != N)
    {
        int n;
        if (!candidates.empty())
        {
            n = candidates.front();
            candidates.pop();
            candidatesUniq.erase(n);
            if (n == maxDegreeNode)
                maxDegreeNode++;
        }
        else
        {
            n = maxDegreeNode;
            maxDegreeNode++;
        }
        if (!marked.test(n))
        {
            if (except[n] == 0)
            {
                marked.set(n);
                markedNum++;
                bar.progress(markedNum, N);
                continue;
            }
            vec &gn = G[n];
            for (const auto &adj : gn)
                adjVec.set(adj); // adjVec 用来判断是否为节点n的邻居
            int prev = 0;
            int minus = 0;
            int hnum = upper - lower + 1;
            R_count += hnum;
            int *currNum = new int[hnum]; // P1
            fill_n(currNum, hnum, 0);
            unsigned triNum = 0; // CN
            //计算score
            for (const int adj : gn)
            {
                if (hashesMin[adj] != INT_MAX)
                {
                    unordered_map<int, unsigned> &hn = hashesNums[adj]; // hn->h(v,.)
                    for (int uh = lower; uh <= hashesMin[adj] - 1; uh++)
                    {
                        if (hn.find(uh) != hn.end())
                            prev += hn[uh]; // prev = p2(u, lower-1) = U(u,h)
                    }
                }
                //计算p1
                //挑选部分节点进行计算
                int d2_count = 0;
                for (const int adj2 : G[adj])
                {
                    if (gn.size() > d1_threshold && d2_count++ > d2_threshold)    break;
                    if (adjVec.test(adj2))
                        triNum++;
                    else if (hashesMin[adj2] != INT_MAX)
                    {
                        currNum[hashesMin[adj2] - lower]++;
                        currNum[hashesMax[adj2] + 1 - lower]--;
                    }
                }
            }
            for (int i = 1; i < hnum; i++)
                currNum[i] = currNum[i - 1] + currNum[i]; // P1(u,h−1) + P1(u,h)
            int minScore = INT_MAX;
            int hstar = upper;
            for (int h = lower; h <= upper; h++)
            {
                int score = triNum + currNum[h - lower] - minus;
                for (const auto &adj : G[n])
                {
                    unordered_map<int, unsigned> &hn = hashesNums[adj];
                    if (hn.find(h) != hn.end())
                    {
                        if (h < hashesMin[adj])
                            minus = minus + hn[h]; // minus = L(u,h-1)
                        if (h > hashesMax[adj])
                            prev = prev + hn[h];
                    }
                }
                score = prev + score;
                if (h == lower || (h != upper && score <= minScore) || (h == upper && score < minScore))
                {
                    minScore = score;
                    hstar = h;
                }
            }
            // cout << hstar << endl;
            hashes[n] = hstar;
            if (hstar == upper)
                upper++;
            else if (hstar == lower)
                lower--;
            updateMinMaxNB(G, n, hstar, hashesNums, candidates, candidatesUniq, hashesMin, hashesMax, d1_threshold, d2_threshold);
            delete[] currNum;
            adjVec.reset();
            marked.set(n);
            markedNum++;
            bar.progress(markedNum, N);
        }
    }
    double R = (double)R_count / (double)N;
    cout << "R = " << R << endl;
    delete[] hashesNums;
}

//更新某个点的range范围
inline void updateMinMaxIndexNB(vec *&G, int node, int h, unordered_map<int, pair<unsigned, unsigned>> *&hashesNums, unordered_set<int> *&hashEntries, queue<int> &candidates, unordered_set<int> &candidatesUniq, int *&hashesMin, int *&hashesMax, int d1_threshold, int d2_threshold, double& part3_time)
{
    auto start_time_part3 = chrono::high_resolution_clock::now();
    for (const auto &n : G[node])
    {
        unordered_set<int> &hashEnts = hashEntries[n];
        bool found = (hashEnts.find(h) != hashEnts.end());//Nh(n)是否存在
        int d2_count = 0;
        for (const auto &n2 : G[n])
        {
            if (G[node].size() > d1_threshold && d2_count++ > d2_threshold) {
                //cout << "d2_count = " << d2_count << endl;
                break;
            }    
            unordered_map<int, pair<unsigned, unsigned>> &hn = hashesNums[n2];
            if (hn.find(h) != hn.end())
            {
                hn[h].first++;
                if (!found)
                    hn[h].second++;
            }
            else
                hn[h] = {1, 1};
        }
        if (h < hashesMin[n])
            hashesMin[n] = h;
        if (h > hashesMax[n])
            hashesMax[n] = h;
        if (!found)
            hashEnts.insert(h);
        if (candidatesUniq.find(n) == candidatesUniq.end())
        {
            candidates.push(n);
            candidatesUniq.insert(n);
        }
    }
    auto end_time_part3 = chrono::high_resolution_clock::now();
    chrono::duration<double> diff_part3 = end_time_part3 - start_time_part3;
    part3_time += diff_part3.count();
}

void hashMergeDiffIndexNBEx(vec *&G, int *&except, int *&hashes, int maxDegreeNode, int N, int d1_threshold, int d2_threshold, double& part1_time, double& part2_time, double& part3_time, int lamda)
{
    unordered_map<int, pair<unsigned, unsigned>> *hashesNums = new unordered_map<int, pair<unsigned, unsigned>>[N]();
    unordered_set<int> *hashEntries = new unordered_set<int>[N]();
    boost::dynamic_bitset<> marked(N);
    marked.set(maxDegreeNode);
    int *hashesMin = new int[N];
    int *hashesMax = new int[N];
    // fill_n函数的作用是：给你一个起始点，然后再给你一个数值count和val。把从起始点开始依次赋予count个元素val的值。
    //将N个点的hashesMin 设置为INT_MAX
    fill_n(hashesMin, N, INT_MAX);
    fill_n(hashesMax, N, INT_MIN);
    boost::dynamic_bitset<> adjVec(N);
    hashes[maxDegreeNode] = 1;
    int upper = 2;
    int lower = 0;
    // vector<unsigned> currNeighbors = adjlist[maxDegreeNode];
    queue<int> candidates;
    unordered_set<int> candidatesUniq;
    updateMinMaxIndexNB(G, maxDegreeNode, 1, hashesNums, hashEntries, candidates, candidatesUniq, hashesMin, hashesMax, d1_threshold, d2_threshold, part3_time);
    tqdm bar;
    bar.progress(1, N);
    int markedNum = 1;
    maxDegreeNode = 1;
    int R_count = 0;
    while (markedNum != N)
    {
        int n;
        if (!candidates.empty())
        {
            n = candidates.front();
            candidates.pop();
            candidatesUniq.erase(n);
            if (n == maxDegreeNode)
                maxDegreeNode++;
        }
        else
            n = maxDegreeNode++;
        if (!marked.test(n))
        {
            if (except[n] == 0)
            {
                marked.set(n);
                markedNum++;
                bar.progress(markedNum, N);
                continue;
            }
            const vec &gn = G[n];
            for (const auto &adj : gn)
                adjVec.set(adj);
            int d = gn.size();
            int *prev = new int[d];
            fill_n(prev, d, 0);
            int *minus = new int[d];
            fill_n(minus, d, 0);
            int hnum = upper - lower + 1;
            R_count += hnum;
            unsigned *adj2num = new unsigned[d];
            fill_n(adj2num, d, 0);
            vector<unsigned *> currNum(d);
            auto start_time_part1 = chrono::high_resolution_clock::now();
            for (int i = 0; i != d; i++)
            {
                int adj = gn[i];
                unordered_map<int, pair<unsigned, unsigned>> &hn = hashesNums[adj];
                if (hashesMin[adj] != INT_MAX)
                {
                    for (int uh = lower; uh <= hashesMin[adj] - 1; uh++)
                    {
                        if (hn.find(uh) != hn.end())
                            prev[i] += hn[uh].second;
                    }
                }
                unsigned *currNumVec = new unsigned[hnum];
                fill_n(currNumVec, hnum, 0);
                int d2_count = 0;
                for (const auto &adj2 : G[adj])
                {
                    if (d > d1_threshold && d2_count++ > d2_threshold)    break;
                    if (adjVec.test(adj2))
                    {
                        for (const int &h : hashEntries[adj2])
                            currNumVec[h - lower]++;
                        adj2num[i]++;
                    }
                }
                currNum[i] = currNumVec;
            }
            auto end_time_part1 = chrono::high_resolution_clock::now();
            chrono::duration<double> diff_part1 = end_time_part1 - start_time_part1;
            part1_time += diff_part1.count();
            int minScore = INT_MAX;
            int hstar = upper;


            auto start_time_part2 = chrono::high_resolution_clock::now();
            for (int h = lower; h <= upper; h++)
            {
                //分成段去计算 
                int score = 0;
                for (int i = 0; i < d; i++)
                {
                    //if (hnum > R_threshold && i > d3_threshold)   break;
                    int adj = gn[i];
                    int currMinus = minus[i];
                    unsigned a2n = adj2num[i];
                    unordered_map<int, pair<unsigned, unsigned>> &hn = hashesNums[adj];
                    bool fd = hn.find(h) != hn.end();
                    int hnf = 0;
                    int hns = 0;
                    if (fd)
                    {
                        hnf = hn[h].first;
                        hns = hn[h].second;
                    }
                    if (fd)
                    {
                        if (hashesMin[adj] == INT_MAX || (hashesMin[adj] <= h && h <= hashesMax[adj]))
                        {
                            minus[i] += hns;
                            prev[i] += hns;
                        }
                        else if (h < hashesMin[adj])
                            minus[i] += hns;
                        else
                            prev[i] += hns;
                    }
                    if (hashEntries[adj].find(h) != hashEntries[adj].end())//Nh(v)存在
                        score += (a2n - currNum[i][h - lower]) + hns;
                    else//Nh(v)不存在
                        score = score + a2n + (a2n - currNum[i][h - lower]) - currMinus + prev[i] + hnf;
                }
                double random = rand() % 100 / (double)100;
                if (h == lower || (h != upper && score <= minScore) || (h == upper && score < minScore && random < (double)1 / (lamda * d)) ||
                    (hstar == lower && random >= (double)1 / (lamda * d)))
                {
                    minScore = score;
                    hstar = h;
                }
                
            }
            auto end_time_part2 = chrono::high_resolution_clock::now();
            chrono::duration<double> diff_part2 = end_time_part2 - start_time_part2;
            part2_time += diff_part2.count();
            hashes[n] = hstar;
            if (hstar == upper)
                upper++;
            else if (hstar == lower)
                lower--;
            updateMinMaxIndexNB(G, n, hstar, hashesNums, hashEntries, candidates, candidatesUniq, hashesMin, hashesMax, d1_threshold, d2_threshold, part3_time);
            adjVec.reset();
            delete[] prev;
            delete[] minus;
            delete[] adj2num;
            for (const auto &hn : currNum)
                delete[] hn;
            marked.set(n);
            markedNum++;
            bar.progress(markedNum, N);
        }
    }
    double R = (double)R_count / (double)N;
    cout << "R = " << R << endl;
    delete[] hashesMin;
    delete[] hashesMax;
    delete[] hashesNums;
    delete[] hashEntries;
}

inline int binarySearch(int _key, const vec &_list, int *&_hashes)
{
    int _left = 0;
    int _right = _list.size();
    int _mid;
    while (_left < _right)
    {
        _mid = (_right + _left) / 2;
        if (_hashes[_list[_mid]] < _key)
            _left = _mid + 1;
        else
            _right = _mid;
    }
    return _left;
}

int main(int argc, char **argv)
{
    string directname = "../data/";
    string outdirect = "../hash_data/d1select_d2limit_Rlimit_diff_part_time/";
    string name = argv[1];  //第一个参数 文件名
    int N = stoi(argv[2]);  //第二个参数 顶点数目
    bool R = stoi(argv[3]); //第三个参数 0 or 1
    int d1_threshold = stoi(argv[4]);
    int d2_threshold = stoi(argv[5]);
    int lamda = stoi(argv[6]);
    int maxDegreeNode = 0;
    int maxDegree=0;
    string filename = directname + name + ".txt";
    outdirect += argv[4];
    outdirect += '_';
    outdirect += argv[5];
    outdirect += '_';
    outdirect += argv[6];
    
    vec *G = new vec[N];
    // G[i][j]表示顶点i对应的邻接顶点j
    readGraph(filename, G, N);
    int *hashes = new int[N];
    fill_n(hashes, N, INT_MAX);
    int *degrees = new int[N];
    fill_n(degrees, N, 0);
    // 从最大度节点开始coding
    for (int i = 0; i < N; i++)
    {
        if(G[i].size()>maxDegree){
            maxDegree=G[i].size();
            maxDegreeNode=i;
        }
    }
    //记录度
    for (int i = 0; i < N; i++)
    {
        if (G[i].size() == 1)
        {
            degrees[G[i][0]] = degrees[G[i][0]] - 1;
            degrees[i] = 0;
        }
        else
            degrees[i] = degrees[i] + G[i].size();
    }
    int *except = new int[N];

    for (int i = 0; i < N; i++)
    {
        if (degrees[i] <= 1)
            except[i] = 0;
        else
        {
            int CN = 0;
            for (int j = 0; j != G[i].size(); j++)
            {
                int n = G[i][j];
                int id2 = binarySearch(n, G[n], 0, G[n].size());
                commonNeighbor(G[i], G[n], CN, j + 1, G[i].size(), id2, G[n].size());
            }
            except[i] = CN; // CN是i点和其所有邻居的公共点数目总和
        }
    }
    if (R)
    {
        int *hashesMin = new int[N];
        int *hashesMax = new int[N];
        fill_n(hashesMin, N, INT_MAX);
        fill_n(hashesMax, N, INT_MIN);
        auto start_time = chrono::high_resolution_clock::now();
        hashMergeDiffNBEx(G, except, hashes, hashesMin, hashesMax, maxDegreeNode, N, d1_threshold, d2_threshold);
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = end_time - start_time;
        cout << diff.count() << " s" << endl;
        delete[] hashesMin;
        delete[] hashesMax;
    }
    else
    {
        double part1_time = 0;
        double part2_time = 0;
        double part3_time = 0;
        auto start_time = chrono::high_resolution_clock::now();
        hashMergeDiffIndexNBEx(G, except, hashes, maxDegreeNode, N, d1_threshold, d2_threshold, part1_time, part2_time, part3_time, lamda);
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = end_time - start_time;
        cout << diff.count() << " s" << endl;
        cout << "part1_time = " << part1_time << endl;
        cout << "part2_time = " << part2_time << endl;
        cout << "part3_time = " << part3_time << endl;
    }
    int *ids = new int[N];
    sort_indexes(hashes, ids, N);
    string outname = name;

    if (R)
        outname += "_range_nb";
    else
        outname += "_index_nb";
    ofstream idfile;

    idfile.open(outdirect + outname + "_id.csv");
    cout << outdirect + outname + "_id.csv" << endl;
    idfile << ":ID" << endl;
    for (int i = 0; i < N; i++)
        idfile << ids[i] << endl;
    ofstream nodefile;
    nodefile.open(outdirect + outname + "_node.csv");
    nodefile << ":ID,h:int" << endl;
    for (int i = 0; i < N; i++)
        nodefile << i << "," << hashes[i] << endl;
    delete[] ids;
    delete[] hashes;
    delete[] except;
    delete[] G;
}