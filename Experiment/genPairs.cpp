#include "indexGraph.h"
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <sys/io.h>
#include <sys/stat.h>
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
int main(int argc, char **argv)
{
    string directname = "../data/";
    string pairdirect = "/root/setintersect/RCode/pair_data/";
    string name = argv[1];
    int N = stoi(argv[2]);
    int M = 0;
    vec *G = new vec[N];
    utils::readGraph(directname + name + ".txt", G, N, M);
    Graph OG(G, N);
    vec node1s = vec(10000);
    vec node2s = vec(10000);
    // Generate local pairs
    srand(time(nullptr));
    auto cnt = 0;
    while (cnt < 10000)
    {
        auto i = rand() % N;
        vec &gi = G[i];
        if (gi.size() != 0)
        {
            auto ind_n = rand() % gi.size();
            auto n = gi[ind_n];
            node1s[cnt] = i;
            node2s[cnt] = n;
            cnt++;
        }
    }
    if (access((pairdirect + name ).c_str(), 0) == -1)
        int re = mkdir((pairdirect + name).c_str(), 0777);
    if (access((pairdirect + name + "/local/").c_str(), 0) == -1)
        int re = mkdir((pairdirect + name + "/local/").c_str(), 0777);
    ofstream out(pairdirect + name + "/local/pair.txt", ios::app);
    auto outcnt = 0;
    while (outcnt < 10000)
    {
        out << node1s[outcnt] << " " << node2s[outcnt] << "\n";
        outcnt++;
    }
    // Generate global pairs
    cnt = 0;
    while (cnt < 10000)
    {
        auto i = rand() % N;
        vec &gi = G[i];
        if (gi.size() != 0)
        {
            auto ind_n = rand() % gi.size();
            auto n = gi[ind_n];
            vec &gn = G[n];
            if (gn.size() != 0)
            {
                auto ind_n2 = rand() % gn.size();
                auto n2 = gn[ind_n2];
                node1s[cnt] = i;
                node2s[cnt] = n2;
                cnt++;
            }
        }
    }
    if (access((pairdirect + name + "/global/").c_str(), 0) == -1)
        int re = mkdir((pairdirect + name + "/global/").c_str(), 0777);
    ofstream out2(pairdirect + name + "/global/pair.txt", ios::app);
    outcnt = 0;
    while (outcnt < 10000)
    {
        out2 << node1s[outcnt] << " " << node2s[outcnt] << "\n";
        outcnt++;
    }
    delete[] G;
}