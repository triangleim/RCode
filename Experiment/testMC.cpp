#include "indexGraph.h"
#include "hashGraph.h"
#include "indexMCGraph.h"
#include <fstream>
typedef unordered_set<int> intset;
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
void getDegOrder(vec *&G, int *&Vrank, int* &order, int N)
{
    order = new int[N];
    // tqdm bar;
    int *degrees = new int[N];
    map<int, intset> D;
    int d;
    for (int i = 0; i < N; i++)
    {
        d = G[i].size();
        degrees[i] = d;
        D[d].insert(i);
    }
    boost::dynamic_bitset<> mark(N);
    int marked = 0;
    int n;
    // int prevd=0;
    // for (map<int, intset>::iterator it=D.begin(); it!=D.end(); it++) {
    //     // if(prevd>it->first)
    //         cout<<it->first<<endl<<flush;
    //     // prevd=it->first;
    //     }
    // cout<<"OK!"<<endl<<flush;
    while (marked < N)
    {
        for (map<int, intset>::iterator it = D.begin(); it != D.end(); it++)
        {
            if (!it->second.empty())
            {
                n = *(it->second.begin());
                break;
            }
        }
        mark.set(n);
        d = degrees[n];
        D[d].erase(n);
        for (const auto &adj : G[n])
        {
            if (!mark.test(n))
            {
                d = degrees[adj];
                D[d].erase(adj);
                D[d - 1].insert(adj);
                degrees[adj] = d - 1;
            }
        }
        Vrank[n] = marked;
        order[marked++] = n;
        // bar.progress(marked, N);
    }
}
int* VrankGlobal;
int main(int argc, char **argv)
{
    string directname = "../data/";

    string outdirect = "../hash_data/";
    string groDirect = "../gro_data/";
    string name = argv[1];
    int N = stoi(argv[2]);
    int PN = stoi(argv[3]);
    int test = stoi(argv[4]);
    if (test == 7)
    {
        directname = "../ordered_data/";
    }
    int times = stoi(argv[5]);
    vec nodes(PN);
    for (int i = 0; i < PN; i++)
        nodes[i] = (N / PN) * i;
    vec *G = new vec[N];
    int M = 0;
    utils::readGraph(directname + name + ".txt", G, N, M);
    int *Vrank = new int[N];
    
    // iota(Vrank,Vrank+N,0);
    int* neworder;
    getDegOrder(G, Vrank, neworder, N);
    Graph OG(G, N);
    double OGT = 0;
    ofstream out("mc.log.csro0.intv", ios::app);
    double origtime = 0, rangetime = 0, indextime = 0;
    for (int i = 0; i < times; i++){
        if(test==-1) break;
        OGT += OG.MC(nodes, Vrank);
    }
    OG.reportRatio(out);
    

    if (test == -1)
    {
        out << "\n>" << name << "<\n"
            << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
        out << "Original Time: " << OGT / times << endl;
        out << "Method\tAvgTime" << endl;
        stringvec groFiles;
        utils::read_directory(groDirect, groFiles);
        for (const string &f : groFiles)
        {
            int idx = f.find(name);       //在aa中查找bb.
            if (idx != std::string::npos) //不存在。
            {
                hashGraph HG(G, N, f, M);
                int *VrankN = new int[N];
                // iota(VrankN,VrankN+N,0);
                for (int i = 0; i != N; i++)
                    VrankN[HG.getId(i)] = Vrank[i];
                vec newNodes(nodes.size());
                for (int i = 0; i != nodes.size(); i++)
                    newNodes[i] = HG.getId(nodes[i]);
                double GROT = 0;
                for (int i = 0; i < times; i++)
                {
                    GROT += HG.MC(newNodes, VrankN);
                    cout << "\r";
                }

                // out << f << " Rate: " << OGT / GROT << endl;
                delete[] VrankN;
                out << get_GroMethod(get_FileBaseName(f)) << "\t" << GROT / times << endl;
            }
        }
        out << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    }
    if (test == -1)
    {
        delete[] G;
        return 0;
    }
    else if (test == 1)
    {
        name = "test_" + name;
        outdirect = "../maxnode_hash_data/";
    }
    else if (test == 2)
    {
        outdirect = "../degree_hash_data/0/";
    }
    else if (test == 3)
    {
        outdirect = "../degree_hash_data/1/";
    }
    else if (test == 4)
    {
        outdirect = "../bfsdegree_hash_data/0/";
    }
    else if (test == 5)
    {
        outdirect = "../bfsdegree_hash_data/1/";
    }
    else if (test == 6)
    {
        outdirect = "../onlyp2_hash_data/";
    }
    else if (test == 7)
    {
        outdirect = "../ordered_hash_data/";
    }
    {
        int *hashes = new int[N];
        iota(hashes, hashes + N, 0);
        // cout<<"Building!\n"<<flush;
        hashGraph HG(G, hashes, N, M);
        // cout<<"Calculation start!\n"<<flush;
        for (int i = 0; i < times; i++)
            origtime += HG.MC(nodes, Vrank);
        // cout << "Orig Order Rate: " << OGT/HGT << endl;
        delete[] hashes;
        HG.reportRatio(out);
        // cout<<"Orignal Order time:"<<HGT<<"s\n";
    }
    string hdirect = outdirect + name + "_range_nb_";
    string hnode = hdirect + "node.csv";
    string hid = hdirect + "id.csv";
    // cout << hnode << endl;
    // cout << hid << endl;
    {
        hashGraph HG(G, N, hnode, hid, M);
        int *VrankN = new int[N];
        for (int i = 0; i != N; i++)
            VrankN[HG.getId(i)] = Vrank[i];
        // iota(VrankN,VrankN+N,0);
        vec newNodes(nodes.size());
        for (int i = 0; i != nodes.size(); i++)
            newNodes[i] = HG.getId(nodes[i]);
        for (int i = 0; i < times; i++)
            rangetime += HG.MC(newNodes, VrankN);
        // cout << "HG Rate: " << OGT / HGT << endl;
        delete[] VrankN;
        HG.reportRatio(out);
        // cout << "Range time:" << HGT << "s\n";
    }
    string idirect = outdirect + name + "_index_nb_";
    string inode = idirect + "node.csv";
    string iid = idirect + "id.csv";
    // cout << inode << endl;
    // cout << iid << endl;
    {
        indexVecGraph IG(G, N, inode, iid, M,true);
        // cout << "Build Index Graph Done" << endl<<flush;
        int *VrankN = new int[N];
        for (int i = 0; i != N; i++)
            VrankN[IG.getId(i)] = Vrank[i];
        // iota(VrankN,VrankN+N,0);
        vec newNodes(nodes.size());
        for (int i = 0; i != nodes.size(); i++)
            newNodes[i] = IG.getId(nodes[i]);
        // cout << "Start" << endl<<flush;
        for (int i = 0; i < times; i++)
            indextime += IG.MC(newNodes, VrankN);
        // cout << "IG Rate: " << OGT / IGT << endl;
        delete[] VrankN;
        IG.reportRatio(out);
        // cout << "Index time:" << IGT << "s\n";
    }

    out << "\n>" << name << "<\n"
        << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    out << OGT / times << "\t" << origtime / times << "\t" << rangetime / times << "\t" << indextime / times << endl;
    // out << OGT / origtime << "\t" << OGT / rangetime << "\t" << OGT / indextime << endl;
    out << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    delete[] Vrank;
    delete[] G;
}