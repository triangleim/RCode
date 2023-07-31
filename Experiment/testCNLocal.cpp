#include "indexGraph.h"

int main(int argc, char **argv)
{
    string directname = "../data/";
    string outdirect = "../hash_data/";
    string pairdirect = "../pair_data/";
    string name = argv[1];
    int N = stoi(argv[2]);
    int M = 0;
    // string PN = argv[3];
    int test = stoi(argv[3]);

    vec *G = new vec[N];
    // // cout << directname + name + ".txt" << endl;
    utils::readGraph(directname + name + ".txt", G, N, M);
    Graph OG(G, N);
    // // cout << "Build Original Graph Done" << endl;
    if (test == 1)
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
    string hdirect = outdirect + name + "_range_nb_";
    string hnode = hdirect + "node.csv";
    string hid = hdirect + "id.csv";
    // cout << hnode << endl;
    // cout << hid << endl;
    hashGraph HG(G, N, hnode, hid, M);
    // cout << "Build Hash Graph Done" << endl;
    string idirect = outdirect + name + "_index_nb_";
    string inode = idirect + "node.csv";
    string iid = idirect + "id.csv";
    // cout << inode << endl;
    // cout << iid << endl;
    indexGraph IG(G, N, inode, iid, M);
    // cout << "Build Index Graph Done" << endl;
    // cout << "--------------------------" << endl;
    if (test == 1)
    {
        name = name.substr(5);
    }
    pairdirect += name + "/local/";
    stringvec pairFiles;
    utils::read_directory(pairdirect, pairFiles);
    double OGF = 0.0;
    double HGF = 0.0;
    double IGF = 0.0;
    double OGT = 0.0;
    double HGT = 0.0;
    double IGT = 0.0;
    for (const string &f : pairFiles)
    {
        // cout << "Start Processing " << f << endl;
        pvec pairs;
        utils::readPair(f, pairs);
        // cout << pairs.size() << endl;
        OGF += OG.calF(pairs, true);
        HGF += HG.calF(pairs, true);
        IGF += IG.calF(pairs, true);
        OGT += OG.calCN(pairs, true);
        HGT += HG.calCN(pairs, true);
        IGT += IG.calCN(pairs, true);
        // cout << "Finish Processing " << f << endl;
    }
    int fnum = pairFiles.size();
    // cout << "Average Orig F: " << OGF / fnum << endl;
    // cout << "Average Hash F: " << HGF / fnum << endl;
    // cout << "Average Index F: " << IGF / fnum << endl;
    // cout << "Average Hash Rate: " << HGR / fnum << endl;
    // cout << "Average Index Rate: " << IGR / fnum << endl;
    ofstream out(string("localcn.log.csro0"), ios::app);
    out<< OGF/ fnum <<"\t"<<HGF/ fnum <<"\t"<<IGF/ fnum <<"\t"<<OGT/ fnum <<"\t"<<HGT/ fnum <<"\t"<<IGT/ fnum<<endl;
    delete[] G;
}