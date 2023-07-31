#include "indexGraph.h"
#include <fstream>
#include <vector>
#include <set>
#include <unordered_map>
;
std::string get_FileBaseName(std::string path)
{
    std::string name;
    for (int i = path.size() - 1; i > 0; i--)
    {
        if (path[i] == '\\' || path[i] == '/')
        {
            name = path.substr(i + 1);
            return name;
        }
    }
    name = path;
    return name;
}
std::string get_GroMethod(std::string path)
{
    std::string name;
    for (int i = 0; i < path.size(); i++)
    {
        if (path[i] == '_')
        {
            name = path.substr(i + 1);
            return name;
        }
    }
    name = path;
    return name;
}
void commonNeighborIf(Graph &g, const int i, const int j, vector<int> &cn)
{

    int id1 = 0, id2 = 0, ed1 = g.degrees[i], ed2 = g.degrees[j];
    vec &vec1 = g.G[i];
    vec &vec2 = g.G[j];
#ifdef TIMECOUNT
    std::chrono::_V2::system_clock::time_point merge_start_time;
    merge_start_time = std::chrono::high_resolution_clock::now();
#endif
    // cout<<"here"<<endl;
    while (id1 != ed1 && id2 != ed2)
    {
        if (vec1[id1] == vec2[id2])
        {
            // if(vec1[id1]>j)
            cn.push_back(vec1[id1]);
            // (*orgcmpr)--;
            id1++;
            id2++;
        }
        else if (vec1[id1] < vec2[id2])
            id1++;
        else
            id2++;
    }
#ifdef TIMECOUNT
    mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
}
double calGen(Graph &g, bool reportNumber, vector<set<int>> &NeoG,vector<unordered_map<int,int>>& NumG,int thres)
{

    tqdm bar;
    auto start_time = chrono::high_resolution_clock::now();
    unsigned long result = 0;
    for (int i = 0; i != g.N; i++)
    {

        vec &adj = g.G[i];
        int s = adj.size();
        vector<int> res;
        for (int j = 0; j != s; j++)
        {
            int n = adj[j];
            if(i>n){
                continue;
            }
            res.clear();
            // if(n>i)
            commonNeighborIf(g, i, n, res);
            for (int k = 1; k < res.size(); k++)
            {

                for (int t = 0; t < k; t++)
                {
                    // cout << t << "\t" << k << endl;
                    auto& tmp =NumG[res[t]][res[k]];
                    if(tmp!=thres)
                        tmp+=1;
                    else
                        NeoG[res[t]].insert(res[k]);
                }
            }
        }

        bar.progress(i + 1, g.N);
    }
    auto end_time = chrono::high_resolution_clock::now();
    // if (reportNumber)
    //     cout << endl
    //          << "Total Triangles:" << result << endl;
    chrono::duration<double> diff = end_time - start_time;
    return diff.count();
}
int main(int argc, char **argv)
{
    string directname = "../part_data/";
    string outdirect = "../neo_data/";
    string groDirect = "../gro_data/";
    string name = argv[1];
    int N = stoi(argv[2]);
    int thres = stoi(argv[3]);
    int times = stoi(argv[4]);
    int M = 0;
    vec *G = new vec[N];
    // ofstream out(string("tri.log.csro0.dyn")+to_string(BINRATIO), ios::app);
    ofstream out(string(outdirect + name + "_local.txt"), ios::app);
    utils::readGraph(directname + name + ".txt", G, N, M);
    Graph OG(G, N);
    double origtime = 0, rangetime = 0, indextime = 0;
    double OGT = 0;
    vector<unordered_map<int,int>> NumG(N);
    vector<set<int>> NeoG(N);
    OGT += calGen(OG, times == 1, NeoG,NumG,thres);
    for (int i = 0; i < N; i++)
    {
        for (auto j : NeoG[i])
            out << i << "\t" << j << "\n";
    }
    delete[] G;
}