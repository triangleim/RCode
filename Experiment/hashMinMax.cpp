#include "indexGraph.h"
#include <fstream>
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
    string outdirect = "../hash_data/";
    string groDirect = "../gro_data/";
    string name = argv[1];
    int N = stoi(argv[2]);
    int test = stoi(argv[3]);
    int times = stoi(argv[4]);
    int M = 0;
    vec *G = new vec[N];
    // ofstream out("tri.log.csro0", ios::app);
    utils::readGraph(directname + name + ".txt", G, N, M);
    Graph OG(G, N);
    double OGT = 0;
    // for (int i = 0; i < times; i++)
    // {
    //     OGT += OG.calTri(times == 1);
    //     cout << "\r";
    // }
    // OG.reportRatio(out);
    double origtime = 0, rangetime = 0, indextime = 0;
    // cout << "Origin Time: " << OGT << endl;
    if (test == -1)
    {
        // out << "\n>" << name << "<\n"
        //     << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
        // out << "Original Time: " << OGT / times << endl;
        // out << "Method\tRate\tAvgTime" << endl;
        stringvec groFiles;
        utils::read_directory(groDirect, groFiles);
        for (const string &f : groFiles)
        {
            int idx = f.find(name);       //
            if (idx != std::string::npos) //不存在。
            {
                hashGraph HG(G, N, f, M);
                double GROT = 0;
                ofstream out(outdirect + name + "_"+ get_GroMethod(get_FileBaseName(f)) + "_nb_minmax.txt", ios::app);
                for (int i = 0; i < N; i++)
                {
                    out << HG.hashesMin[i] << "," << HG.hashesMax[i] << endl;
                }
                // out << get_GroMethod(get_FileBaseName(f)) << "\t" << OGT / GROT << "\t" << GROT / times << endl;
            }
        }
        // out << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
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
    // {
    //     int *hashes = new int[N];
    //     iota(hashes, hashes + N, 0);
    //     hashGraph HG(G, hashes, N, M);
    //     ofstream out(outdirect + name + "_orig_nb_minmax.txt", ios::app);
    //     for (int i = 0; i < N; i++)
    //     {
    //         out << HG.hashesMin[i] << "," << HG.hashesMax[i] << endl;
    //     }
    //     // HG.reportRatio(out);
    //     // cout << "Orig Order Rate: " << OGT / HGT << endl;
    //     // cout << "Orig Order Time: " << HGT << endl;
    //     delete[] hashes;
    // }
    {
        string hdirect = outdirect + name + "_range_nb_";
        string hnode = hdirect + "node.csv";
        string hid = hdirect + "id.csv";
        // cout << hnode << endl;
        // cout << hid << endl;
        hashGraph HG(G, N, hnode, hid, M);
        ofstream out(outdirect + name+"_range_nb_minmax.txt", ios::app);
        for (int i = 0; i < N; i++)
        {
            out << HG.hashesMin[i] << "," << HG.hashesMax[i] << endl;
        }
        // HG.reportRatio(out);
        // cout << "HG Rate: " << OGT / HGT << endl;
        // cout << "Range Time: " << HGT << endl;
    }
    if (test != 6)
    {
        string idirect = outdirect + name + "_index_nb_";
        string inode = idirect + "node.csv";
        string iid = idirect + "id.csv";
        // cout << inode << endl;
        // cout << iid << endl;
        indexGraph IG(G, N, inode, iid, M);
        ofstream out(outdirect + name+ "_index_nb_minmax.txt", ios::app);
        for (int i = 0; i < N; i++)
        {
            out << IG.hashesMin[i] << "," << IG.hashesMax[i] << endl;
        }
        // IG.reportRatio(out);
        // cout << "IG Rate: " << OGT / IGT << endl;
        // cout << "Index Time: " << IGT << endl;
    }
    else
    {
        indextime = 1;
    }

    // out << "\n>" << name << "<\n"
    //     << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    // out << OGT / times << "\t" << origtime / times << "\t" << rangetime / times << "\t" << indextime / times << endl;
    // out << OGT / origtime << "\t" << OGT / rangetime << "\t" << OGT / indextime << endl;
    // out << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    delete[] G;
}