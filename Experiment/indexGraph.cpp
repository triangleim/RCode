#include "indexGraph.h"

void indexGraph::processHGIndex(const vec *OG)
{
    // cout << "\nHere\n"
    //      << flush;
    tqdm bar;
    oa = new int[N + 1];
    hashOA = new int[N + 1];
    ca = new int[2 * M];
    degreesh = new int[N + 1];
    for (int i = 0; i < N + 1; i++)
    {
        degreesh[i] = -1;
    }
    // hashIdsA = new int *[N];
    // hashEntriesA = new int *[N];
    int hashcount = 0;
    for (int i = 0; i < N; i++)
    {
        int newid = ids[i];
        G[newid].reserve(OG[i].size());
        for (const int &adj : OG[i])
            G[newid].push_back(ids[adj]);
        sort(G[newid].begin(), G[newid].end());
        int prev = INT_MIN;
        int j = 0;
        while (j != G[newid].size() && hashes[newid] != INT_MAX)
        {
            int adj = G[newid][j];
            if (hashes[adj] != prev)
            {
                prev = hashes[adj];
                if (prev == INT_MAX)
                {
                    degreesh[newid] = j;
                    break;
                }
                hashEntries[newid].push_back(prev);
                hashIds[newid].push_back(j);
                hashcount++;
            }
            j++;
        }
        if (degreesh[newid] == -1)
        {
            degreesh[newid] = j;
        }
        if (!hashEntries[newid].empty())
        {
            hashesMin[newid] = hashEntries[newid].front();
            hashesMax[newid] = hashEntries[newid].back();
        }
        else
        {
            hashesMin[newid] = INT_MAX;
            hashesMax[newid] = INT_MIN;
        }
        degrees[newid] = hashEntries[newid].size();
        hashEntries[newid].push_back(INT_MAX);

        hashIds[newid].push_back(j);
        hashcount++;
        // hashIdsA[newid] = new int[hashIds[newid].size()];
        // // hashEntriesA[newid]=new int[hashEntries[newid].size()];
        // for (int i = 0; i < hashIds[newid].size(); i++)
        // {
        //     hashIdsA[newid][i] = hashIds[newid][i];
        // }
        // for(int i = 0;i<hashEntries[newid].size();i++){
        //     hashEntriesA[newid][i]=hashEntries[newid][i];
        // }
        bar.progress(i + 1, N);
    }
    int cur_start = 0;
    for (int i = 0; i < N; i++)
    {
        oa[i] = cur_start;
        for (auto nb : G[i])
        {
            ca[cur_start++] = nb;
        }
        // cout<<"\r"<<marked<<"/"<<N<<flush;
        // bar.progress(++marked, N);
    }
    oa[N] = cur_start;

    int hash_cur_start = 0;
    hashEntriesCA = new int[hashcount];
    hashIdsCA = new int[hashcount];
    for (int i = 0; i < N; i++)
    {
        hashOA[i] = hash_cur_start;
        for (int j = 0; j < hashEntries[i].size(); j++)
        {
            hashIdsCA[hash_cur_start] = hashIds[i][j];
            hashEntriesCA[hash_cur_start++] = hashEntries[i][j];
        }
        // cout<<"\r"<<marked<<"/"<<N<<flush;
        // bar.progress(++marked, N);
    }
    hashOA[N] = hash_cur_start;
}
void indexGraph::processHGIndexMC(const vec *OG)
{
    // cout << "\nHere\n"
    //      << flush;
    tqdm bar;
    oa = new int[N + 1];
    hashOA = new int[N + 1];
    ca = new int[2 * M];
    degreesh = new int[N + 1];
    for (int i = 0; i < N + 1; i++)
    {
        degreesh[i] = -1;
    }
    // hashIdsA = new int *[N];
    // hashEntriesA = new int *[N];
    int hashcount = 0;
    for (int i = 0; i < N; i++)
    {
        int newid = ids[i];
        G[newid].reserve(OG[i].size());
        for (const int &adj : OG[i])
            G[newid].push_back(ids[adj]);
        sort(G[newid].begin(), G[newid].end());
        int prev = INT_MIN;
        int j = 0;
        while (j != G[newid].size())
        {
            int adj = G[newid][j];
            if (hashes[adj] != prev)
            {
                prev = hashes[adj];
                if (prev == INT_MAX && degreesh[newid] == -1)
                {
                    degreesh[newid] = j;
                }
                hashEntries[newid].push_back(prev);
                hashIds[newid].push_back(j);
                hashcount++;
            }
            j++;
        }
        if (degreesh[newid] == -1)
        {
            degreesh[newid] = j;
        }
        if (!hashEntries[newid].empty())
        {
            hashesMin[newid] = hashEntries[newid].front();
            hashesMax[newid] = hashEntries[newid].back();
        }
        else
        {
            hashesMin[newid] = INT_MAX;
            hashesMax[newid] = INT_MIN;
        }
        if (hashesMax[newid] == INT_MAX)
        {
            if (hashEntries[newid].size() > 1)
            {
                hashesMax[newid] = hashEntries[newid][hashEntries[newid].size() - 2];
            }
            else if (hashEntries[newid].size() == 1)
            {
                hashesMin[newid] = INT_MAX;
                hashesMax[newid] = INT_MIN;
            }
            degrees[newid] = hashEntries[newid].size() - 1;
        }
        else
        {
            degrees[newid] = hashEntries[newid].size();
            hashEntries[newid].push_back(INT_MAX);

            hashIds[newid].push_back(j);
            hashcount++;
        }
        // hashIdsA[newid] = new int[hashIds[newid].size()];
        // // hashEntriesA[newid]=new int[hashEntries[newid].size()];
        // for (int i = 0; i < hashIds[newid].size(); i++)
        // {
        //     hashIdsA[newid][i] = hashIds[newid][i];
        // }
        // for(int i = 0;i<hashEntries[newid].size();i++){
        //     hashEntriesA[newid][i]=hashEntries[newid][i];
        // }
        bar.progress(i + 1, N);
    }
    int cur_start = 0;
    for (int i = 0; i < N; i++)
    {
        oa[i] = cur_start;
        for (auto nb : G[i])
        {
            ca[cur_start++] = nb;
        }
        // cout<<"\r"<<marked<<"/"<<N<<flush;
        // bar.progress(++marked, N);
    }
    oa[N] = cur_start;

    int hash_cur_start = 0;
    hashEntriesCA = new int[hashcount];
    hashIdsCA = new int[hashcount];
    for (int i = 0; i < N; i++)
    {
        hashOA[i] = hash_cur_start;
        for (int j = 0; j < hashEntries[i].size(); j++)
        {
            hashIdsCA[hash_cur_start] = hashIds[i][j];
            hashEntriesCA[hash_cur_start++] = hashEntries[i][j];
        }
        // cout<<"\r"<<marked<<"/"<<N<<flush;
        // bar.progress(++marked, N);
    }
    hashOA[N] = hash_cur_start;
}
extern unsigned long long mins;
extern unsigned long long maxs;
extern double avgs;
double avgs2 = 0;
extern unsigned long long times;
extern unsigned long long actcmpr;
extern unsigned long long orgcmpr;
extern chrono::duration<double> mergetime;
extern chrono::duration<double> bintime;
void indexGraph::reportRatio(ostream &out)
{

    cout<<"Reduced ratio:"<<(float)mins/maxs<<endl<<flush;
    // cout<<"Reduced ratio(AVG):"<<avgs/times<<endl<<flush;
    // if (orgcmpr)
    //     cout << "Reduced ratio(Actual):" << (float)actcmpr / orgcmpr << endl
    //          << flush;
    // cout<<"Binary Search Time:"<<bintime.count()<<endl<<flush;
    // out << "Index Merge Time:" << mergetime.count() << endl;
    // //     << flush;
    // out << "#P-Hash/#P-Entry: " << mins / (double)maxs << "\t#X-Hash/#X-Entry: " << actcmpr / (double)orgcmpr << endl;
    // out << "AVG #P-Hash/#P-Entry: " << avgs / times << "\tAVG #X-Hash/#X-Entry: " << avgs2 / times << endl;
    mins = 0;
    maxs = 0;
    avgs = 0;
    times = 0;
    actcmpr = 0;
    mergetime = chrono::duration<double>(0);
    bintime = chrono::duration<double>(0);
}
#if CNV == 1
void indexGraph::commonNeighbor(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    int m1 = hashesMin[i];
    int m2 = hashesMin[n];
    int M1 = hashesMax[i];
    int M2 = hashesMax[n];
#ifndef CSR
    vec &vec1 = G[i];
    vec &vec2 = G[n];
    int s1 = vec1.size();
    int s2 = vec2.size();
#else
    int *vec1 = ca + oa[i];
    int *vec2 = ca + oa[n];
    int s1 = oa[i + 1] - oa[i];
    int s2 = oa[n + 1] - oa[n];
#endif
    if (s1 == 1 && l1 == 0)
    {
        l2 = hashIds[n][l2];
        int v = vec1[0];
        if (s2 == 1 && v == vec2[0])
            result++;
        else if (l2 < s2 && v >= vec2[l2] && v <= vec2[s2 - 1])
        {
            int lb = utils::binarySearch(v, vec2, l2, s2);
            if (lb != s2 && vec2[lb] == v)
                result++;
        }
    }
    else if (s2 == 1 && l2 == 0 && l1 < s1)
    {
        l1 = hashIds[i][l1];
        int v = vec2[0];
        if (v >= vec1[l1] && v <= vec1[s1 - 1])
        {
            int lb = utils::binarySearch(v, vec1, l1, s1);
            if (lb != s1 && vec1[lb] == v)
                result++;
        }
    }
    else if (l1 != 0 && l2 != 0)
    {
#ifndef CSR
        vec &hashEnts1 = hashEntries[i];
        vec &hashEnts2 = hashEntries[n];
        vec &hashId1 = hashIds[i];
        vec &hashId2 = hashIds[n];
#else
        int *hashEnts1 = hashEntriesCA + hashOA[i];
        int *hashEnts2 = hashEntriesCA + hashOA[n];
        int *hashId1 = hashIdsCA + hashOA[i];
        int *hashId2 = hashIdsCA + hashOA[n];
#endif
        if (l1 == 1)
        {
            if (l2 == 1 && vec1[0] == vec2[0])
                result++;
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[l2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, l2);
                    if (lb != l2 && vec2[lb] == v)
                        result++;
                }
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[s1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (max(m1, m2) <= min(M1, M2))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            maxs+=s1+s2;
            if (m1 < m2)
                b1 = utils::binarySearch(m2, hashEnts1, b1, f1);
            else if (m1 > m2)
                b2 = utils::binarySearch(m1, hashEnts2, b2, f2);
            if (M1 > M2)
                f1 = utils::binarySearch(M2 + 1, hashEnts1, b1, f1);
            else if (M1 < M2)
                f2 = utils::binarySearch(M1 + 1, hashEnts2, b2, f2);
            mins+=f2-b2+f1-b1;
            while (b1 != f1 && b2 != f2)
            {
                if (hashEnts1[b1] == hashEnts2[b2])
                {
                    int vit1 = hashId1[b1];
                    int vit2 = hashId2[b2];
                    int ved1 = hashId1[b1 + 1];
                    int ved2 = hashId2[b2 + 1];
                    mins+=ved2-vit2+ved1-vit1;
                    while (vit1 != ved1 && vit2 != ved2)
                    {
                        if (vec1[vit1] == vec2[vit2])
                        {
                            result++;
                            vit1++;
                            vit2++;
                        }
                        else if (vec1[vit1] < vec2[vit2])
                            vit1++;
                        else
                            vit2++;
                    }
                    b1++;
                    b2++;
                }
                else if (hashEnts1[b1] < hashEnts2[b2])
                    b1++;
                else
                    b2++;
            }
        }
    }
}
#elif CNV == 2
void indexGraph::commonNeighbor(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    int m1 = hashesMin[i];
    int m2 = hashesMin[n];
    int M1 = hashesMax[i];
    int M2 = hashesMax[n];
    if (l1 != 0 && l2 != 0)
    {
        vec &vec1 = G[i];
        vec &vec2 = G[n];
        vec &hashEnts1 = hashEntries[i];
        vec &hashEnts2 = hashEntries[n];
        vec &hashId1 = hashIds[i];
        vec &hashId2 = hashIds[n];
        int s1 = hashId1[l1];
        int s2 = hashId2[l2];
        if (s1 == 1)
        {
            if (s2 == 1 && vec1[0] == vec2[0])
                result++;
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[s2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, s2);
                    if (lb != s2 && vec2[lb] == v)
                        result++;
                }
            }
        }
        else if (s2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[s1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, s1);
                if (lb != s1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (max(m1, m2) <= min(M1, M2))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (m1 < m2)
                b1 = utils::binarySearch(m2, hashEnts1, b1, f1);
            else if (m1 > m2)
                b2 = utils::binarySearch(m1, hashEnts2, b2, f2);
            if (M1 > M2)
                f1 = utils::binarySearch(M2 + 1, hashEnts1, b1, f1);
            else if (M1 < M2)
                f2 = utils::binarySearch(M1 + 1, hashEnts2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (hashEnts1[b1] == hashEnts2[b2])
                {
                    int vit1 = hashId1[b1];
                    int vit2 = hashId2[b2];
                    int ved1 = hashId1[b1 + 1];
                    int ved2 = hashId2[b2 + 1];
                    while (vit1 != ved1 && vit2 != ved2)
                    {
                        if (vec1[vit1] == vec2[vit2])
                        {
                            result++;
                            vit1++;
                            vit2++;
                        }
                        else if (vec1[vit1] < vec2[vit2])
                            vit1++;
                        else
                            vit2++;
                    }
                    b1++;
                    b2++;
                }
                else if (hashEnts1[b1] < hashEnts2[b2])
                    b1++;
                else
                    b2++;
            }
        }
    }
}
#elif CNV == 3
void indexGraph::commonNeighbor(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    vec &vec1 = G[i];
    vec &vec2 = G[n];
    int s1 = vec1.size();
    int s2 = vec2.size();
    int m1 = hashesMin[i];
    int m2 = hashesMin[n];
    int M1 = hashesMax[i];
    int M2 = hashesMax[n];
    if (s1 == 1 && l1 == 0)
    {
        l2 = hashIds[n][l2];
        int v = vec1[0];
        if (s2 == 1 && v == vec2[0])
            result++;
        else if (l2 < s2 && v >= vec2[l2] && v <= vec2[s2 - 1])
        {
            int lb = utils::binarySearch(v, vec2, l2, s2);
            if (lb != s2 && vec2[lb] == v)
                result++;
        }
    }
    else if (s2 == 1 && l2 == 0 && l1 < s1)
    {
        l1 = hashIds[i][l1];
        int v = vec2[0];
        if (v >= vec1[l1] && v <= vec1[s1 - 1])
        {
            int lb = utils::binarySearch(v, vec1, l1, s1);
            if (lb != s1 && vec1[lb] == v)
                result++;
        }
    }
    else if (l1 != 0 && l2 != 0 && max(m1, m2) <= min(M1, M2))
    {
        vec &hashEnts1 = hashEntries[i];
        vec &hashEnts2 = hashEntries[n];
        vec &hashId1 = hashIds[i];
        vec &hashId2 = hashIds[n];
        int b1 = 0;
        int f1 = l1;
        int b2 = 0;
        int f2 = l2;
        if (m1 < m2)
            b1 = utils::binarySearch(m2, hashEnts1, b1, f1);
        else if (m1 > m2)
            b2 = utils::binarySearch(m1, hashEnts2, b2, f2);
        if (M1 > M2)
            f1 = utils::binarySearch(M2 + 1, hashEnts1, b1, f1);
        else if (M1 < M2)
            f2 = utils::binarySearch(M1 + 1, hashEnts2, b2, f2);
        while (b1 != f1 && b2 != f2)
        {
            if (hashEnts1[b1] == hashEnts2[b2])
            {
                int vit1 = hashId1[b1];
                int vit2 = hashId2[b2];
                int ved1 = hashId1[b1 + 1];
                int ved2 = hashId2[b2 + 1];
                while (vit1 != ved1 && vit2 != ved2)
                {
                    if (vec1[vit1] == vec2[vit2])
                    {
                        result++;
                        vit1++;
                        vit2++;
                    }
                    else if (vec1[vit1] < vec2[vit2])
                        vit1++;
                    else
                        vit2++;
                }
                b1++;
                b2++;
            }
            else if (hashEnts1[b1] < hashEnts2[b2])
                b1++;
            else
                b2++;
        }
    }
}
#endif

#if CNNBV == 3
void indexGraph::commonNeighborNB(const int i1, const int n1, unsigned long &result)
{
    int i = i1, n = n1;
    // times+=1;
    // if(i1>n1){
    //     i=n1;
    //     n=i1;
    // }
    int l1 = degrees[i];
    int l2 = degrees[n];
    // maxs+=G[i].size()+G[n].size();
    // int cmprs=0;
    if (l1 != 0 && l2 != 0)
    {
        int m1 = hashesMin[i];
        int m2 = hashesMin[n];
        int M1 = hashesMax[i];
        int M2 = hashesMax[n];
#ifndef CSR
        vec &vec1 = G[i];
        vec &vec2 = G[n];
        vec &hashEnts1 = hashEntries[i];
        vec &hashEnts2 = hashEntries[n];
        vec &hashId1 = hashIds[i];
        vec &hashId2 = hashIds[n];
        int s1 = vec1.size();
        int s2 = vec2.size();
#else
        int *vec1 = ca + oa[i];
        int *vec2 = ca + oa[n];
        int *hashEnts1 = hashEntriesCA + hashOA[i];
        int *hashEnts2 = hashEntriesCA + hashOA[n];
        int *hashId1 = hashIdsCA + hashOA[i];
        int *hashId2 = hashIdsCA + hashOA[n];
        int s1 = oa[i + 1] - oa[i];
        int s2 = oa[n + 1] - oa[n];
#endif
        if (s1 == 1)
        {
            if (s2 == 1 && vec1[0] == vec2[0])
            {
                result++;
                cout << "Something strange happened!" << endl
                     << flush;
            }
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[s2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, s2);
                    if (lb != s2 && vec2[lb] == v)
                    {
                        result++;
                        cout << "Something strange happened!" << endl
                             << flush;
                    }
                }
            }
        }
        else if (s2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[s1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, s1);
                if (lb != s1 && vec1[lb] == v)
                {
                    result++;
                    cout << "Something strange happened!" << endl
                         << flush;
                }
            }
        }
        else if (max(m1, m2) <= min(M1, M2))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            // auto bin_start_time = chrono::high_resolution_clock::now();
            if (m1 < m2)
                b1 = utils::binarySearch(m2, hashEnts1, b1, f1);
            else if (m1 > m2)
                b2 = utils::binarySearch(m1, hashEnts2, b2, f2);
            if (M1 > M2)
                f1 = utils::binarySearch(M2 + 1, hashEnts1, b1, f1);
            else if (M1 < M2)
                f2 = utils::binarySearch(M1 + 1, hashEnts2, b2, f2);
            // bintime += chrono::high_resolution_clock::now()-bin_start_time;
            // mins+=f1-b1+f2-b2;
            // actcmpr+=f1-b1+f2-b2;
            // cmprs+=f1-b1+f2-b2;
            // auto merge_start_time = std::chrono::high_resolution_clock::now();
            // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
            while (b1 != f1 && b2 != f2)
            {
                if (hashEnts1[b1] == hashEnts2[b2])
                {

                    // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
                    int vit1 = hashId1[b1];
                    int vit2 = hashId2[b2];
                    int ved1 = hashId1[b1 + 1];
                    int ved2 = hashId2[b2 + 1];

                    // actcmpr++;
                    // mins+=ved1-vit1+ved2-vit2;
                    // actcmpr+=ved1-vit1+ved2-vit2-1;
                    // cmprs+=ved1-vit1+ved2-vit2;
                    while (vit1 != ved1 && vit2 != ved2)
                    {
                        if (vec1[vit1] == vec2[vit2])
                        {
                            result++;
                            // actcmpr--;
                            // actcmpr++;
                            vit1++;
                            vit2++;
                        }
                        else if (vec1[vit1] < vec2[vit2])
                        {
                            // if ((ved1 - vit1) > 32 * (ved2 - vit2))
                            vit1 = utils::binarySearch(vec2[vit2], vec1, vit1, ved1);
                            // else
                            // {
                            //     while (vec1[vit1] < vec2[vit2])
                            //         vit1++;
                            // }
                            // actcmpr++;
                        }
                        else
                        {
                            // if ((ved2 - vit2) > 32 * (ved1 - vit1))
                            vit2 = utils::binarySearch(vec1[vit1], vec2, vit2, ved2);
                            // else
                            // {
                            //     while (vec1[vit1] > vec2[vit2])
                            //         vit2++;
                            // }
                            // actcmpr++;
                        }
                        if (vit1 == ved1 || vit2 == ved2)
                        {
                            break;
                        }
                        else if (vec1[ved1 - 1] == vec2[ved2 - 1])
                        {
                            result++;
                            ved1--;
                            ved2--;
                        }
                        else if (vec1[ved1 - 1] > vec2[ved2 - 1])
                        {
                            // if ((ved1 - vit1) > 32 * (ved2 - vit2))
                            ved1 = utils::binarySearch(vec2[ved2 - 1] + 1, vec1, vit1, ved1);
                            // else
                            // {
                            //     while (vec1[ved1 - 1] > vec2[ved2 - 1])
                            //     {
                            //         ved1--;
                            //     }
                            // }
                        }
                        else
                        {
                            // if ((ved2 - vit2) > 32 * (ved1 - vit1))
                            ved2 = utils::binarySearch(vec1[ved1 - 1] + 1, vec2, vit2, ved2);
                            // else
                            // {
                            //     while (vec1[ved1 - 1] < vec2[ved2 - 1])
                            //     {
                            //         ved2--;
                            //     }
                            // }
                        }
                    }
                    // actcmpr-=ved1-vit1+ved2-vit2;
                    b1++;
                    b2++;
                }
                else if (hashEnts1[b1] < hashEnts2[b2])
                {
                    // if ((f1 - b1) > 32 * (f2 - b2))
                    b1 = utils::binarySearch(hashEnts2[b2], hashEnts1, b1, f1);
                    // else
                    // {
                    //     while (hashEnts1[b1] < hashEnts2[b2])
                    //         b1++;
                    // }

                    // actcmpr++;
                }
                else
                {
                    // if ((f2 - b2) > 32 * (f1 - b1))
                    b2 = utils::binarySearch(hashEnts1[b1], hashEnts2, b2, f2);
                    // else
                    // {
                    //     while (hashEnts1[b1] > hashEnts2[b2])
                    //         b2++;
                    // }
                    // actcmpr++;
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                if (hashEnts1[f1 - 1] == hashEnts2[f2 - 1])
                {

                    // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
                    int vit1 = hashId1[f1 - 1];
                    int vit2 = hashId2[f2 - 1];
                    int ved1 = hashId1[f1];
                    int ved2 = hashId2[f2];

                    // actcmpr++;
                    // mins+=ved1-vit1+ved2-vit2;
                    // actcmpr+=ved1-vit1+ved2-vit2-1;
                    // cmprs+=ved1-vit1+ved2-vit2;
                    while (vit1 != ved1 && vit2 != ved2)
                    {
                        if (vec1[vit1] == vec2[vit2])
                        {
                            result++;
                            // actcmpr--;
                            // actcmpr++;
                            vit1++;
                            vit2++;
                        }
                        else if (vec1[vit1] < vec2[vit2])
                        {
                            vit1 = utils::binarySearch(vec2[vit2], vec1, vit1, ved1);
                            // actcmpr++;
                        }
                        else
                        {
                            vit2 = utils::binarySearch(vec1[vit1], vec2, vit2, ved2);
                            // actcmpr++;
                        }
                        if (vit1 == ved1 || vit2 == ved2)
                        {
                            break;
                        }
                        else if (vec1[ved1 - 1] == vec2[ved2 - 1])
                        {
                            result++;
                            ved1--;
                            ved2--;
                        }
                        else if (vec1[ved1 - 1] > vec2[ved2 - 1])
                            ved1 = utils::binarySearch(vec2[ved2 - 1] + 1, vec1, vit1, ved1);
                        else
                            ved2 = utils::binarySearch(vec1[ved1 - 1] + 1, vec2, vit2, ved2);
                    }
                    // actcmpr-=ved1-vit1+ved2-vit2;
                    f1--;
                    f2--;
                }
                else if (hashEnts1[f1 - 1] > hashEnts2[f2 - 1])
                {
                    // if ((f1 - b1) > 32 * (f2 - b2))
                    f1 = utils::binarySearch(hashEnts2[f2 - 1] + 1, hashEnts1, b1, f1);
                    // else
                    // {
                    //     while (hashEnts1[f1 - 1] > hashEnts2[f2 - 1])
                    //     {
                    //         f1--;
                    //     }
                    // }

                    // actcmpr++;
                }
                else
                {
                    // if ((f2 - b2) > 32 * (f1 - b1))
                    f2 = utils::binarySearch(hashEnts1[f1 - 1] + 1, hashEnts2, b2, f2);
                    // else
                    // {
                    //     while (hashEnts1[f1 - 1] < hashEnts2[f2 - 1])
                    //     {
                    //         f2--;
                    //     }
                    // }

                    // actcmpr++;
                }
            }

            // actcmpr-=f1-b1+f2-b2;
        }
    }
    // avgs+=(float)cmprs/(G[i].size()+G[n].size());
}
#elif CNNBV == 4
void indexGraph::commonNeighborNB(const int i1, const int n1, unsigned long &result)
{
    int i = i1, n = n1;
    // times+=1;
    // if(i1>n1){
    //     i=n1;
    //     n=i1;
    // }
    int l1 = degrees[i];
    int l2 = degrees[n];
    // maxs+=G[i].size()+G[n].size();
    // int cmprs=0;
    if (l1 != 0 && l2 != 0)
    {
        int m1 = hashesMin[i];
        int m2 = hashesMin[n];
        int M1 = hashesMax[i];
        int M2 = hashesMax[n];
#ifndef CSR
        vec &vec1 = G[i];
        vec &vec2 = G[n];
        vec &hashEnts1 = hashEntries[i];
        vec &hashEnts2 = hashEntries[n];
        vec &hashId1 = hashIds[i];
        vec &hashId2 = hashIds[n];
        int s1 = vec1.size();
        int s2 = vec2.size();
#else
        int *vec1 = ca + oa[i];
        int *vec2 = ca + oa[n];
        int *hashEnts1 = hashEntriesCA + hashOA[i];
        int *hashEnts2 = hashEntriesCA + hashOA[n];
        int *hashId1 = hashIdsCA + hashOA[i];
        int *hashId2 = hashIdsCA + hashOA[n];
        int s1 = oa[i + 1] - oa[i];
        int s2 = oa[n + 1] - oa[n];
#endif
        if (s1 == 1)
        {
            if (s2 == 1 && vec1[0] == vec2[0])
            {
                result++;
                cout << "Something strange happened!" << endl
                     << flush;
            }
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[s2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, s2);
                    if (lb != s2 && vec2[lb] == v)
                    {
                        result++;
                        cout << "Something strange happened!" << endl
                             << flush;
                    }
                }
            }
        }
        else if (s2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[s1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, s1);
                if (lb != s1 && vec1[lb] == v)
                {
                    result++;
                    cout << "Something strange happened!" << endl
                         << flush;
                }
            }
        }
        else if (max(m1, m2) <= min(M1, M2))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            // auto bin_start_time = chrono::high_resolution_clock::now();
            if (m1 < m2)
                b1 = utils::binarySearch(m2, hashEnts1, b1, f1);
            else if (m1 > m2)
                b2 = utils::binarySearch(m1, hashEnts2, b2, f2);
            if (M1 > M2)
                f1 = utils::binarySearch(M2 + 1, hashEnts1, b1, f1);
            else if (M1 < M2)
                f2 = utils::binarySearch(M1 + 1, hashEnts2, b2, f2);
            // bintime += chrono::high_resolution_clock::now()-bin_start_time;
            // mins+=f1-b1+f2-b2;
            // actcmpr+=f1-b1+f2-b2;
            // cmprs+=f1-b1+f2-b2;
            // auto merge_start_time = std::chrono::high_resolution_clock::now();
            // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
            while (b1 != f1 && b2 != f2)
            {
                if (hashEnts1[b1] == hashEnts2[b2])
                {

                    // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
                    int vit1 = hashId1[b1];
                    int vit2 = hashId2[b2];
                    int ved1 = hashId1[b1 + 1];
                    int ved2 = hashId2[b2 + 1];

                    // actcmpr++;
                    // mins+=ved1-vit1+ved2-vit2;
                    // actcmpr+=ved1-vit1+ved2-vit2-1;
                    // cmprs+=ved1-vit1+ved2-vit2;
                    while (vit1 != ved1 && vit2 != ved2)
                    {
                        if (vec1[vit1] == vec2[vit2])
                        {
                            result++;
                            // actcmpr--;
                            // actcmpr++;
                            vit1++;
                            vit2++;
                        }
                        else if (vec1[vit1] < vec2[vit2])
                        {
                            if ((ved1 - vit1) > BINRATIO * (ved2 - vit2))
                                vit1 = utils::binarySearch(vec2[vit2], vec1, vit1, ved1);
                            else
                            {
                                while (vec1[vit1] < vec2[vit2])
                                    vit1++;
                            }
                            // actcmpr++;
                        }
                        else
                        {
                            if ((ved2 - vit2) > BINRATIO * (ved1 - vit1))
                                vit2 = utils::binarySearch(vec1[vit1], vec2, vit2, ved2);
                            else
                            {
                                while (vec1[vit1] > vec2[vit2])
                                    vit2++;
                            }
                            // actcmpr++;
                        }
                        if (vit1 == ved1 || vit2 == ved2)
                        {
                            break;
                        }
                        else if (vec1[ved1 - 1] == vec2[ved2 - 1])
                        {
                            result++;
                            ved1--;
                            ved2--;
                        }
                        else if (vec1[ved1 - 1] > vec2[ved2 - 1])
                        {
                            if ((ved1 - vit1) > BINRATIO * (ved2 - vit2))
                                ved1 = utils::binarySearch(vec2[ved2 - 1] + 1, vec1, vit1, ved1);
                            else
                            {
                                while (vec1[ved1 - 1] > vec2[ved2 - 1])
                                {
                                    ved1--;
                                }
                            }
                        }
                        else
                        {
                            if ((ved2 - vit2) > BINRATIO * (ved1 - vit1))
                                ved2 = utils::binarySearch(vec1[ved1 - 1] + 1, vec2, vit2, ved2);
                            else
                            {
                                while (vec1[ved1 - 1] < vec2[ved2 - 1])
                                {
                                    ved2--;
                                }
                            }
                        }
                    }
                    // actcmpr-=ved1-vit1+ved2-vit2;
                    b1++;
                    b2++;
                }
                else if (hashEnts1[b1] < hashEnts2[b2])
                {
                    if ((f1 - b1) > BINRATIO * (f2 - b2))
                        b1 = utils::binarySearch(hashEnts2[b2], hashEnts1, b1, f1);
                    else
                    {
                        while (hashEnts1[b1] < hashEnts2[b2])
                            b1++;
                    }

                    // actcmpr++;
                }
                else
                {
                    if ((f2 - b2) > BINRATIO * (f1 - b1))
                        b2 = utils::binarySearch(hashEnts1[b1], hashEnts2, b2, f2);
                    else
                    {
                        while (hashEnts1[b1] > hashEnts2[b2])
                            b2++;
                    }
                    // actcmpr++;
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                if (hashEnts1[f1 - 1] == hashEnts2[f2 - 1])
                {

                    // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
                    int vit1 = hashId1[f1 - 1];
                    int vit2 = hashId2[f2 - 1];
                    int ved1 = hashId1[f1];
                    int ved2 = hashId2[f2];

                    // actcmpr++;
                    // mins+=ved1-vit1+ved2-vit2;
                    // actcmpr+=ved1-vit1+ved2-vit2-1;
                    // cmprs+=ved1-vit1+ved2-vit2;
                    while (vit1 != ved1 && vit2 != ved2)
                    {
                        if (vec1[vit1] == vec2[vit2])
                        {
                            result++;
                            // actcmpr--;
                            // actcmpr++;
                            vit1++;
                            vit2++;
                        }
                        else if (vec1[vit1] < vec2[vit2])
                        {
                            if ((ved1 - vit1) > BINRATIO * (ved2 - vit2))
                                vit1 = utils::binarySearch(vec2[vit2], vec1, vit1, ved1);
                            else
                            {
                                while (vec1[vit1] < vec2[vit2])
                                    vit1++;
                            }
                            // actcmpr++;
                        }
                        else
                        {
                            if ((ved2 - vit2) > BINRATIO * (ved1 - vit1))
                                vit2 = utils::binarySearch(vec1[vit1], vec2, vit2, ved2);
                            else
                            {
                                while (vec1[vit1] > vec2[vit2])
                                    vit2++;
                            }
                            // actcmpr++;
                        }
                        if (vit1 == ved1 || vit2 == ved2)
                        {
                            break;
                        }
                        else if (vec1[ved1 - 1] == vec2[ved2 - 1])
                        {
                            result++;
                            ved1--;
                            ved2--;
                        }
                        else if (vec1[ved1 - 1] > vec2[ved2 - 1])
                        {
                            if ((ved1 - vit1) > BINRATIO * (ved2 - vit2))
                                ved1 = utils::binarySearch(vec2[ved2 - 1] + 1, vec1, vit1, ved1);
                            else
                            {
                                while (vec1[ved1 - 1] > vec2[ved2 - 1])
                                {
                                    ved1--;
                                }
                            }
                        }
                        else
                        {
                            if ((ved2 - vit2) > BINRATIO * (ved1 - vit1))
                                ved2 = utils::binarySearch(vec1[ved1 - 1] + 1, vec2, vit2, ved2);
                            else
                            {
                                while (vec1[ved1 - 1] < vec2[ved2 - 1])
                                {
                                    ved2--;
                                }
                            }
                        }
                    }
                    // actcmpr-=ved1-vit1+ved2-vit2;
                    f1--;
                    f2--;
                }
                else if (hashEnts1[f1 - 1] > hashEnts2[f2 - 1])
                {
                    if ((f1 - b1) > BINRATIO * (f2 - b2))
                        f1 = utils::binarySearch(hashEnts2[f2 - 1] + 1, hashEnts1, b1, f1);
                    else
                    {
                        while (hashEnts1[f1 - 1] > hashEnts2[f2 - 1])
                        {
                            f1--;
                        }
                    }

                    // actcmpr++;
                }
                else
                {
                    if ((f2 - b2) > BINRATIO * (f1 - b1))
                        f2 = utils::binarySearch(hashEnts1[f1 - 1] + 1, hashEnts2, b2, f2);
                    else
                    {
                        while (hashEnts1[f1 - 1] < hashEnts2[f2 - 1])
                        {
                            f2--;
                        }
                    }

                    // actcmpr++;
                }
            }

            // actcmpr-=f1-b1+f2-b2;
        }
    }
    // avgs+=(float)cmprs/(G[i].size()+G[n].size());
}
#elif CNNBV == 1
void indexGraph::commonNeighborNB(const int i1, const int n1, unsigned long &result)
{
    int i = i1, n = n1;
    // times+=1;
    // if(i1>n1){
    //     i=n1;
    //     n=i1;
    // }
    int l1 = degrees[i];
    int l2 = degrees[n];
    maxs+=G[i].size()+G[n].size();
    // int cmprs=0;
    if (l1 != 0 && l2 != 0)
    {
        int m1 = hashesMin[i];
        int m2 = hashesMin[n];
        int M1 = hashesMax[i];
        int M2 = hashesMax[n];
#ifndef CSR
        vec &vec1 = G[i];
        vec &vec2 = G[n];
        vec &hashEnts1 = hashEntries[i];
        vec &hashEnts2 = hashEntries[n];
        vec &hashId1 = hashIds[i];
        vec &hashId2 = hashIds[n];
        int s1 = vec1.size();
        int s2 = vec2.size();
#else
        int *vec1 = ca + oa[i];
        int *vec2 = ca + oa[n];
        int *hashEnts1 = hashEntriesCA + hashOA[i];
        int *hashEnts2 = hashEntriesCA + hashOA[n];
        int *hashId1 = hashIdsCA + hashOA[i];
        int *hashId2 = hashIdsCA + hashOA[n];
        int s1 = oa[i + 1] - oa[i];
        int s2 = oa[n + 1] - oa[n];
#endif
        times++;
        mins += hashOA[i + 1] - hashOA[i]+hashOA[n + 1] - hashOA[n];
        // maxs += s1+s2;
        actcmpr += hashOA[n + 1] - hashOA[n];
        orgcmpr += s2;
        avgs += (double)(hashOA[i + 1] - hashOA[i]) / s1;
        avgs2 += (double)(hashOA[n + 1] - hashOA[n]) / s2;
        if (s1 == 1)
        {
            if (s2 == 1 && vec1[0] == vec2[0])
            {
                result++;
                cout << "Something strange happened!" << endl
                     << flush;
            }
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[s2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, s2);
                    if (lb != s2 && vec2[lb] == v)
                    {
                        result++;
                        cout << "Something strange happened!" << endl
                             << flush;
                    }
                }
            }
        }
        else if (s2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[s1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, s1);
                if (lb != s1 && vec1[lb] == v)
                {
                    result++;
                    cout << "Something strange happened!" << endl
                         << flush;
                }
            }
        }
        else if (max(m1, m2) <= min(M1, M2))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
// auto bin_start_time = chrono::high_resolution_clock::now();
#ifdef TIMECOUNT
            std::chrono::_V2::system_clock::time_point merge_start_time;
            merge_start_time = std::chrono::high_resolution_clock::now();
#endif
            if (m1 < m2)
                b1 = utils::binarySearch(m2, hashEnts1, b1, f1);
            else if (m1 > m2)
                b2 = utils::binarySearch(m1, hashEnts2, b2, f2);
            if (M1 > M2)
                f1 = utils::binarySearch(M2 + 1, hashEnts1, b1, f1);
            else if (M1 < M2)
                f2 = utils::binarySearch(M1 + 1, hashEnts2, b2, f2);
            // bintime += chrono::high_resolution_clock::now()-bin_start_time;
            // mins+=f1-b1+f2-b2;
            // actcmpr+=f1-b1+f2-b2;
            // cmprs+=f1-b1+f2-b2;
            // auto merge_start_time = std::chrono::high_resolution_clock::now();
            // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
            while (b1 != f1 && b2 != f2)
            {
                if (hashEnts1[b1] == hashEnts2[b2])
                {

                    // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
                    int vit1 = hashId1[b1];
                    int vit2 = hashId2[b2];
                    int ved1 = hashId1[b1 + 1];
                    int ved2 = hashId2[b2 + 1];

                    // actcmpr++;
                    mins+=ved1-vit1+ved2-vit2;
                    // actcmpr+=ved1-vit1+ved2-vit2-1;
                    // cmprs+=ved1-vit1+ved2-vit2;
                    while (vit1 != ved1 && vit2 != ved2)
                    {
                        if (vec1[vit1] == vec2[vit2])
                        {
                            // if(vec1[vit1]>n)
                            result++;
                            // actcmpr--;
                            // actcmpr++;
                            vit1++;
                            vit2++;
                        }
                        else if (vec1[vit1] < vec2[vit2])
                        {
                            vit1++;
                            // actcmpr++;
                        }
                        else
                        {
                            vit2++;
                            // actcmpr++;
                        }
                    }
                    // actcmpr-=ved1-vit1+ved2-vit2;
                    b1++;
                    b2++;
                }
                else if (hashEnts1[b1] < hashEnts2[b2])
                {
                    b1++;
                    // actcmpr++;
                }
                else
                {
                    b2++;
                    // actcmpr++;
                }
            }
#ifdef TIMECOUNT
            mergetime += chrono::high_resolution_clock::now() - merge_start_time;
#endif

            // actcmpr-=f1-b1+f2-b2;
        }
    }
    // avgs+=(float)cmprs/(G[i].size()+G[n].size());
}
#elif CNNBV == 2
void indexGraph::commonNeighborNB(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    int m1 = hashesMin[i];
    int m2 = hashesMin[n];
    int M1 = hashesMax[i];
    int M2 = hashesMax[n];
    if (l1 != 0 && l2 != 0 && max(m1, m2) <= min(M1, M2))
    {
        vec &vec1 = G[i];
        vec &vec2 = G[n];
        vec &hashEnts1 = hashEntries[i];
        vec &hashEnts2 = hashEntries[n];
        vec &hashId1 = hashIds[i];
        vec &hashId2 = hashIds[n];
        int b1 = 0;
        int f1 = l1;
        int b2 = 0;
        int f2 = l2;
        if (m1 < m2)
            b1 = utils::binarySearch(m2, hashEnts1, b1, f1);
        else if (m1 > m2)
            b2 = utils::binarySearch(m1, hashEnts2, b2, f2);
        if (M1 > M2)
            f1 = utils::binarySearch(M2 + 1, hashEnts1, b1, f1);
        else if (M1 < M2)
            f2 = utils::binarySearch(M1 + 1, hashEnts2, b2, f2);
        while (b1 != f1 && b2 != f2)
        {
            if (hashEnts1[b1] == hashEnts2[b2])
            {
                int vit1 = hashId1[b1];
                int vit2 = hashId2[b2];
                int ved1 = hashId1[b1 + 1];
                int ved2 = hashId2[b2 + 1];
                while (vit1 != ved1 && vit2 != ved2)
                {
                    if (vec1[vit1] == vec2[vit2])
                    {
                        result++;
                        vit1++;
                        vit2++;
                    }
                    else if (vec1[vit1] < vec2[vit2])
                        vit1++;
                    else
                        vit2++;
                }
                b1++;
                b2++;
            }
            else if (hashEnts1[b1] < hashEnts2[b2])
                b1++;
            else
                b2++;
        }
    }
}
#endif
void indexGraph::commonNeighborNB(const int *n1, const int d1, const int *n2, const int d2, unsigned long &result, vec &hashEnts1, int *hashEnts2, vec &hashId1, int *hashId2)
{
    if (d1 == 1)
    {
        int h = hashEnts1[0];
        if (d2 == 1)
        {
            if (h == hashEnts2[0])
                utils::commonNeighbor(n1, n2, result, 0, hashId1[1], 0, hashId2[1]);
        }
        else
        {
            int lb = utils::binarySearch(h, hashEnts2, 0, d2);
            if (lb != d2 && hashEnts2[lb] == h)
                utils::commonNeighbor(n1, n2, result, 0, hashId1[1], hashId2[lb], hashId2[lb + 1]);
        }
    }
    else if (d2 == 1)
    {
        int h = hashEnts2[0];
        int lb = utils::binarySearch(h, hashEnts1, 0, d1);
        if (lb != d1 && hashEnts1[lb] == h)
            utils::commonNeighbor(n1, n2, result, hashId1[lb], hashId1[lb + 1], 0, hashId2[1]);
    }
    else
    {
        int b1 = 0;
        int f1 = d1;
        int b2 = 0;
        int f2 = d2;
        int m1 = hashEnts1[0];
        int m2 = hashEnts2[0];
        int M1 = hashEnts1[d1 - 1];
        int M2 = hashEnts2[d2 - 1];
        if (max(m1, m2) <= min(M1, M2))
        {
            if (m1 < m2)
                b1 = utils::binarySearch(m2, hashEnts1, b1, f1);
            else if (m2 < m1)
                b2 = utils::binarySearch(m1, hashEnts2, b2, f2);
            if (M1 > M2)
                f1 = utils::binarySearch(M2 + 1, hashEnts1, b1, f1);
            else if (M1 < M2)
                f2 = utils::binarySearch(M1 + 1, hashEnts2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (hashEnts1[b1] == hashEnts2[b2])
                {
                    utils::commonNeighbor(n1, n2, result, hashId1[b1], hashId1[1 + b1], hashId2[b2], hashId2[1 + b2]);
                    b1++;
                    b2++;
                }
                else if (hashEnts1[b1] < hashEnts2[b2])
                    b1++;
                else
                    b2++;
            }
        }
    }
}
void indexGraph::commonNeighborNB(const int *n1, const int d1, const int *n2, const int d2, vec &result, vec &hashEnts1, int *hashEnts2, vec &hashId1, int *hashId2, vec &hashEntsR, vec &hashIdR)
{
    if (d1 == 1)
    {
        int h = hashEnts1[0];
        if (d2 == 1)
        {
            if (h == hashEnts2[0])
            {
                utils::commonNeighbor(n1, n2, result, 0, hashId1[1], 0, hashId2[1]);
                if (result.size())
                {
                    hashEntsR.push_back(h);

                    hashIdR.push_back(0);
                }
            }
        }
        else
        {
            int lb = utils::binarySearch(h, hashEnts2, 0, d2);
            if (lb != d2 && hashEnts2[lb] == h)
            {
                utils::commonNeighbor(n1, n2, result, 0, hashId1[1], hashId2[lb], hashId2[lb + 1]);
                if (result.size())
                {
                    hashEntsR.push_back(h);
                    hashIdR.push_back(0);
                }
            }
        }
    }
    else if (d2 == 1)
    {
        int h = hashEnts2[0];
        int lb = utils::binarySearch(h, hashEnts1, 0, d1);
        if (lb != d1 && hashEnts1[lb] == h)
        {
            utils::commonNeighbor(n1, n2, result, hashId1[lb], hashId1[lb + 1], 0, hashId2[1]);
            if (result.size())
            {
                hashEntsR.push_back(h);

                hashIdR.push_back(0);
            }
        }
    }
    else
    {
        int b1 = 0;
        int f1 = d1;
        int b2 = 0;
        int f2 = d2;
        int m1 = hashEnts1[0];
        int m2 = hashEnts2[0];
        int M1 = hashEnts1[d1 - 1];
        int M2 = hashEnts2[d2 - 1];
        if (max(m1, m2) <= min(M1, M2))
        {

            if (m1 < m2)
                b1 = utils::binarySearch(m2, hashEnts1, b1, f1);
            else if (m2 < m1)
                b2 = utils::binarySearch(m1, hashEnts2, b2, f2);
            if (M1 > M2)
                f1 = utils::binarySearch(M2 + 1, hashEnts1, b1, f1);
            else if (M1 < M2)
                f2 = utils::binarySearch(M1 + 1, hashEnts2, b2, f2);

            while (b1 != f1 && b2 != f2)
            {
                // cout << m1 << " " << m2 << " " << M1 << " " << M2 << " " << endl;
                int h = hashEnts1[b1];
                if (h == hashEnts2[b2])
                {
                    int start = result.size();
                    utils::commonNeighbor(n1, n2, result, hashId1[b1], hashId1[1 + b1], hashId2[b2], hashId2[1 + b2]);
                    if (result.size() != start)
                    {

                        hashEntsR.push_back(h);
                        hashIdR.push_back(start);
                    }
                    b1++;
                    b2++;
                }
                else if (h < hashEnts2[b2])
                    b1++;
                else
                    b2++;
            }
        }
    }
    hashEntsR.push_back(INT_MAX);
    hashIdR.push_back(result.size());
}
int indexGraph::selectPivot(const vec &PN, const int ps, const vec &XN, const int mp, const int Mp, vec &hashEntsP, vec &hashEntsX, vec &hashIdP, vec &hashIdX)
{
    int p = PN[0];
    int si = 0;
    int pss = PN.size();
    // for (const entry &adjl : PN)
    // {
    for (const int v : PN) // v是Node
    {

        int gvs = degreesh[v];
        if (gvs != 0)
        {
            unsigned long cn = 0;
            if (max(mp, hashesMin[v]) <= min(Mp, hashesMax[v]))
            {
#ifdef TIMECOUNT
                std::chrono::_V2::system_clock::time_point merge_start_time;
                merge_start_time = std::chrono::high_resolution_clock::now();
#endif
                //                 if (ps == 1)
                // #ifdef CSR
                //                     commonNeighborNBP(*(PN[0].second), PN[0].first, ca + oa[v], degrees[v], cn);
                // #else
                //                     commonNeighborNBP(*(PN[0].second), PN[0].first, G[v], degrees[v], cn);
                // #endif
                //                 else
                // hashGraph::commonNeighborNB(PN.data(), ca + oa[v], cn, pss, degreesh[v], max(mp, hashes[*(ca + oa[v])]), min(Mp, hashes[*(ca + oa[v] + degreesh[v] - 1)]));
                commonNeighborNB(PN.data(), ps, ca + oa[v], degrees[v], cn, hashEntsP, hashEntriesCA + hashOA[v], hashIdP, hashIdsCA + hashOA[v]);
#ifdef TIMECOUNT
                mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
                if (cn > si)
                {
                    si = cn;
                    p = v;
                }
            }
        }
        // commonNeighborNBP(PN, ps, G[v], degrees[v], cn);
    }
    // }
    // for (const entry &adjl : XN)
    // {
    for (const int v : XN)
    {

        int gvs = degreesh[v];
        if (gvs != 0)
        {
            unsigned long cn = 0;
            if (max(mp, hashesMin[v]) <= min(Mp, hashesMax[v]))
            {
#ifdef TIMECOUNT
                std::chrono::_V2::system_clock::time_point merge_start_time;
                merge_start_time = std::chrono::high_resolution_clock::now();
#endif
                //                 if (ps == 1)
                // #ifdef CSR
                //                     commonNeighborNBP(*(PN[0].second), PN[0].first, ca + oa[v], degrees[v], cn);
                // #else
                //                     commonNeighborNBP(*(PN[0].second), PN[0].first, G[v], degrees[v], cn);
                // #endif
                //                 else
                // hashGraph::commonNeighborNB(PN.data(), ca + oa[v], cn, pss, degreesh[v], max(mp, hashes[*(ca + oa[v])]), min(Mp, hashes[*(ca + oa[v] + degreesh[v] - 1)]));
                commonNeighborNB(PN.data(), ps, ca + oa[v], degrees[v], cn, hashEntsP, hashEntriesCA + hashOA[v], hashIdP, hashIdsCA + hashOA[v]);
#ifdef TIMECOUNT
                mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
                if (cn > si)
                {
                    si = cn;
                    p = v;
                }
            }
        }
        // commonNeighborNBP(PN, ps, G[v], degrees[v], cn);
    }
    // }
    return p;
}
int indexGraph::selectPivot(const vec &P, const vec &X, const int mp, const int Mp)
{
    // return Graph::selectPivot(P,X);
    int p = P[0];
    int si = 0;
    int ps = P.size();
    for (const int &v : P)
    {
        int gvs = degreesh[v];
        if (gvs != 0)
        {
            unsigned long cn = 0;
#ifdef TIMECOUNT
            std::chrono::_V2::system_clock::time_point merge_start_time;
            merge_start_time = std::chrono::high_resolution_clock::now();
#endif
#ifndef CSR
            commonNeighborNB(P, G[v], cn, ps, gvs, max(mp, hashes[G[v][0]]), min(Mp, hashes[G[v][gvs - 1]]));
#else
            hashGraph::commonNeighborNB(P.data(), ca + oa[v], cn, ps, gvs, max(mp, hashes[*(ca + oa[v])]), min(Mp, hashes[*(ca + oa[v] + gvs - 1)]));
#endif

#ifdef TIMECOUNT
            mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
            if (cn > si)
            {
                si = cn;
                p = v;
            }
        }
    }
    for (const int &v : X)
    {
        int gvs = degreesh[v];
        if (gvs != 0)
        {
            unsigned long cn = 0;
#ifdef TIMECOUNT
            std::chrono::_V2::system_clock::time_point merge_start_time;
            merge_start_time = std::chrono::high_resolution_clock::now();
#endif
#ifndef CSR
            commonNeighborNB(P, G[v], cn, ps, gvs, max(mp, hashes[G[v][0]]), min(Mp, hashes[G[v][gvs - 1]]));
#else
            hashGraph::commonNeighborNB(P.data(), ca + oa[v], cn, ps, gvs, max(mp, hashes[*(ca + oa[v])]), min(Mp, hashes[*(ca + oa[v] + gvs - 1)]));
#endif
#ifdef TIMECOUNT
            mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
            if (cn > si)
            {
                si = cn;
                p = v;
            }
        }
    }
    return p;
}
void indexGraph::BKP(vec &P, vec &X, unsigned long &result)
{
    // Graph::BKP(P,X,result);
    // return;
    std::chrono::_V2::system_clock::time_point merge_start_time;
    merge_start_time = std::chrono::high_resolution_clock::now();
    if (P.empty())
    {
        if (X.empty())
            result++;
        return;
    }
    int ps = P.size();
    int xs = X.size();
    int u = selectPivot(P, X, hashes[P[0]], hashes[P[ps - 1]]);
    int vit1 = 0;
    int vit2 = 0;
    int ud = degreesh[u];
    mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
    while (vit1 != ps)
    {
        if (vit2 == ud || P[vit1] < G[u][vit2])
        {
            merge_start_time = std::chrono::high_resolution_clock::now();
            int v = P[vit1];
            int gvs = degreesh[v];
            vec NP;
            vec NX;
            if (gvs != 0)
            {
                NP.reserve(min(ps, gvs));
                NX.reserve(min(xs, gvs));
#ifdef TIMECOUNT
                std::chrono::_V2::system_clock::time_point merge_start_time;
                merge_start_time = std::chrono::high_resolution_clock::now();
#endif
                hashGraph::commonNeighborNB(P.data(), ca + oa[v], NP, ps, gvs, max(hashes[P[0]], hashes[*(ca + oa[v])]), min(hashes[P[ps - 1]], hashes[*(ca + oa[v] + gvs - 1)]));
                if (xs != 0)
                    hashGraph::commonNeighborNB(X.data(), ca + oa[v], NX, xs, gvs, max(hashes[X[0]], hashes[*(ca + oa[v])]), min(hashes[X[xs - 1]], hashes[*(ca + oa[v] + gvs - 1)]));
#ifdef TIMECOUNT
                mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
            }
            mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
            BKP(NP, NX, result);
            merge_start_time = std::chrono::high_resolution_clock::now();
            P.erase(P.begin() + vit1);
            X.insert(X.begin() + utils::binarySearch(v, X, 0, xs), v);
            ps--;
            xs++;
            mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
        }
        else if (P[vit1] == G[u][vit2])
        {
            vit1++;
            vit2++;
        }
        else
            vit2++;
    }
}
void indexGraph::BKP2(vec &PN, vec &XN, vec &hashEntsP, vec &hashEntsX, vec &hashIdP, vec &hashIdX, unsigned long &result)
{

    return;
}
void indexGraph::BKP(vec &PN, vec &XN, vec &hashEntsP, vec &hashEntsX, vec &hashIdP, vec &hashIdX, unsigned long &result)
{

    if (PN.empty())
    {
        if (XN.empty())
            result++;
        // else
        // {
        //     for (entry &x : XN)
        //         delete x.second;
        // }
        return;
    }
    // std::chrono::_V2::system_clock::time_point merge_start_time;
    // merge_start_time = std::chrono::high_resolution_clock::now();
    int ps = hashEntsP.size() - 1;
    int xs = hashEntsX.size() - 1;
    // int bxs = xs;
    int pns = PN.size();
    if (xs < 0)
    {
        xs = 0;
        hashEntsX.push_back(INT_MAX);
        hashIdX.push_back(0);
    }
    mins += ps;
    maxs += PN.size();
    actcmpr += xs;
    orgcmpr += XN.size();
    times++;
    avgs += (double)ps / PN.size();
    avgs2 += XN.empty() ? 0 : (double)xs / XN.size();
    // if(xs<0){
    //     cout<<"wrongg"<<endl;
    // }
    int u = selectPivot(PN, ps, XN, hashes[PN[0]], hashes[PN[pns - 1]], hashEntsP, hashEntsX, hashIdP, hashIdX);
    // std::chrono::_V2::system_clock::time_point merge_start_time;
    // merge_start_time = std::chrono::high_resolution_clock::now();
    // int u1 = selectPivot(PN, XN, hashes[PN[0]], hashes[PN[pns - 1]]);
    // if (u1 != u)
    // {
    //     cout << u << "\t" << u1 << endl;
    // }
    int ud = degreesh[u]; // u的邻居hash数
    // vec PMu;
    // PMu.reserve(min(pns, ud));
    int vit1 = 0;
    int vit2 = 0;
    // int ud = oa[u + 1] - oa[u];
    // utils::setMinus(PN,G[u],PMu);
    // setMinus(PN.data(), ps, ca+oa[u], degrees[u], PMu,hashEntsP, hashEntriesCA + hashOA[u], hashIdP, hashIdsCA + hashOA[u]);
    // mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
    // for (const int v : PMu)
    // {

    while (vit1 != pns) // PN-N(u)
    {
        if (vit2 == ud || PN[vit1] < G[u][vit2])
        {
            // std::chrono::_V2::system_clock::time_point merge_start_time;
            // merge_start_time = std::chrono::high_resolution_clock::now();
            int v = PN[vit1];
            int *gv = ca + oa[v];
            int gvs = degrees[v];
            // adjList NP;
            // adjList NX;
            vec NP;
            vec NX;
            vec hashEntsNP, hashEntsNX, hashIdNP, hashIdNX;
            if (gvs != 0)
            {

                if (max(hashEntsP[0], hashesMin[v]) <= min(hashEntsP[ps - 1], hashesMax[v]))
                {
                    hashEntsNP.reserve(min(ps, gvs));
                    hashIdNP.reserve(min(ps, gvs));
                    NP.reserve(min(int(PN.size()), oa[v + 1] - oa[v]));
#ifdef TIMECOUNT
                    std::chrono::_V2::system_clock::time_point merge_start_time;
                    merge_start_time = std::chrono::high_resolution_clock::now();
#endif
                    //                 if (ps == 1)
                    // #ifndef CSR
                    //                     commonNeighborNBP(*(PN[0].second), PN[0].first, G[v], degrees[v], NP);
                    // #else
                    //                     commonNeighborNBP(*(PN[0].second), PN[0].first, ca + oa[v], degrees[v], NP);
                    // #endif
                    //                 else
                    commonNeighborNB(PN.data(), ps, gv, gvs, NP, hashEntsP, hashEntriesCA + hashOA[v], hashIdP, hashIdsCA + hashOA[v], hashEntsNP, hashIdNP); // NP 是PN和G[v]的所有交集
#ifdef TIMECOUNT
                    mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
                    // commonNeighborNBP(PN, ps, gv, gvs, NP);
                }

                if (xs != 0 && max(hashEntsX[0], hashesMin[v]) <= min(hashEntsX[xs - 1], hashesMax[v]))
                {

                    hashEntsNX.reserve(min(xs, gvs));
                    hashIdNX.reserve(min(xs, gvs));
                    NX.reserve(min(int(XN.size()), oa[v + 1] - oa[v]));
#ifdef TIMECOUNT
                    std::chrono::_V2::system_clock::time_point merge_start_time;
                    merge_start_time = std::chrono::high_resolution_clock::now();
#endif
                    //                 if (xs == 1)
                    // #ifdef CSR
                    //                     commonNeighborNBP(*(XN[0].second), XN[0].first, ca + oa[v], degrees[v], NX);
                    // #else
                    //                     commonNeighborNBP(*(XN[0].second), XN[0].first, G[v], degrees[v], NX);
                    // #endif
                    //                 else

                    commonNeighborNB(XN.data(), xs, gv, gvs, NX, hashEntsX, hashEntriesCA + hashOA[v], hashIdX, hashIdsCA + hashOA[v], hashEntsNX, hashIdNX); // NX 是XN和G[v]的所有交集
                                                                                                                                                              // if (hashEntsNX.size() == 0)
                                                                                                                                                              // {
                                                                                                                                                              //     cout << "wrong!" << endl;
                                                                                                                                                              // }
#ifdef TIMECOUNT
                    mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
                }
            }
            // mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
            if ((NP.size() / (ps + 1) > 4 || NX.size() / (xs + 1) > 4) || (NP.size() / (ps + 1) > 2 && NX.size() / (xs + 1) > 2))
                BKP(NP, NX, hashEntsNP, hashEntsNX, hashIdNP, hashIdNX, result); // 反复调用回落成一层的
            else
                BKP(NP, NX, result);
            // merge_start_time = std::chrono::high_resolution_clock::now();
            int h = hashes[v];
            int hpit = utils::binarySearch(h, hashEntsP, 0, ps);
            //
            // utils::binarySearch(v, PN, hashIdP[hpit], hashIdP[hpit + 1]);
            //
            PN.erase(PN.begin() + utils::binarySearch(v, PN, hashIdP[hpit], hashIdP[hpit + 1]));
            pns--;

            if (hashIdP[hpit + 1] - hashIdP[hpit] == 1)
            {

                for (int i = hpit; i < hashIdP.size() - 1; i++)
                {
                    hashIdP[i] = hashIdP[i + 1] - 1;
                    hashEntsP[i] = hashEntsP[i + 1];
                }
                hashIdP.pop_back();
                hashEntsP.pop_back();
                ps--;
                // cout << hashEntsP.size() << "\t" << ps<<endl;
            }
            else
            {

                for (int i = hpit + 1; i < hashIdP.size(); i++)
                {
                    hashIdP[i]--;
                }
            }

            int hxit = utils::binarySearch(h, hashEntsX, 0, xs);

            if (hxit != xs && hashEntsX[hxit] == h)
            {

                XN.insert(XN.begin() + utils::binarySearch(v, XN, hashIdX[hxit], hashIdX[hxit + 1]), v);
                for (int i = hxit + 1; i < hashIdX.size(); i++)
                {
                    hashIdX[i]++;
                }
            }
            else
            {
                int newpos = hashIdX[hxit];
                hashIdX.insert(hashIdX.begin() + hxit, newpos);
                hashEntsX.insert(hashEntsX.begin() + hxit, h);
                XN.insert(XN.begin() + newpos, v);
                for (int i = hxit + 1; i < hashIdX.size(); i++)
                {
                    hashIdX[i]++;
                }
                xs++;
            }
            // mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
        }
        else if (PN[vit1] == G[u][vit2])
        {
            // merge_start_time = std::chrono::high_resolution_clock::now();
            vit1++;
            vit2++;
            // mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
        }
        else
            vit2++;
    }

    // for (entry &p : PN)
    //     delete p.second;
    // for (entry &x : XN)
    //     delete x.second;
}

double indexGraph::MC(const vec &nodes, const int *Vrank)
{
    mergetime = std::chrono::duration<double>(0);
    tqdm bar;
    unsigned long result = 0;
    auto start_time = chrono::high_resolution_clock::now();
    int count = 0;
    for (const int &n : nodes)
    {
        bar.progress(count++, nodes.size());
        int r = Vrank[n];
        int gs = degrees[n]; // 邻居Hash个数
        vec PN;
        vec XN;
        vec hashEntsP, hashEntsX, hashIdP, hashIdX;
        hashEntsP.reserve(gs);
        hashEntsX.reserve(gs);
        hashIdP.reserve(gs);
        hashIdX.reserve(gs);
        int *hashEntsN = hashEntriesCA + hashOA[n];
        int *hashIdN = hashIdsCA + hashOA[n];
        int *GN = ca + oa[n];
        PN.reserve(oa[n + 1] - oa[n]);
        XN.reserve(oa[n + 1] - oa[n]);
        for (int i = 0; i != hashOA[n + 1] - hashOA[n]; i++)
        {
            int h = hashEntsN[i]; // 当前h
            int *gn = GN + hashIdN[i];
            bool hM = h == INT_MAX;
            int gns = hashIdN[i + 1] - hashIdN[i];
            if (hM)
            {
                gns = oa[n + 1] - oa[n] - hashIdN[i];
            }
            bool pf = true;
            bool xf = true;

            for (int adjcount = 0; adjcount < gns; adjcount++)
            {
                const int adj = gn[adjcount];
                if (Vrank[adj] > r)
                {
                    // if (hashes[adj] != h)
                    // {
                    //     cout << "Wrong!" << endl;
                    // }
                    if (hM)
                    {
                        result++;
                    }
                    else if (pf)
                    {
                        hashEntsP.push_back(h);
                        hashIdP.push_back(PN.size());
                        PN.push_back(adj);
                        pf = false;
                    }
                    else
                        PN.push_back(adj);
                }
                else if (!hM)
                {
                    // if (hashes[adj] != h)
                    // {
                    //     cout << "Wrong!" << endl;
                    // }
                    if (xf)
                    {
                        hashEntsX.push_back(h);
                        hashIdX.push_back(XN.size());
                        // XN.push_back(adj);
                        xf = false;
                    }
                    XN.push_back(adj);
                }
            }
        }
        hashEntsP.push_back(INT_MAX);
        hashEntsX.push_back(INT_MAX);
        hashIdP.push_back(PN.size());
        hashIdX.push_back(XN.size());
        if (PN.empty() && XN.empty())
            continue;

        // mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
        if ((PN.size() / hashIdP.size() > 4 || XN.size() / hashIdX.size() > 4) || (PN.size() / hashIdP.size() > 2 && XN.size() / hashIdX.size() > 2))
            BKP(PN, XN, hashEntsP, hashEntsX, hashIdP, hashIdX, result);
        else
            BKP(PN, XN, result);
    }
    auto end_time = chrono::high_resolution_clock::now();
    cout << endl
         << "Total MC:" << result << endl;
    chrono::duration<double> diff = end_time - start_time;
    return diff.count();
}

pair<int, int> indexGraph::calCP2(const int i, const int n, const int m, const int M)
{
    int result = 0;
    // vec &hashEnts1 = hashEntries[i];
    // vec &hashEnts2 = hashEntries[n];
    // vec &hashId1 = hashIds[i];
    // vec &hashId2 = hashIds[n];
    // int l1 = degrees[i];
    // int l2 = degrees[n];
    // int b1 = 0;
    // int f1 = l1;
    // int b2 = 0;
    // int f2 = l2;
    // if (hashesMin[i] < hashesMin[n])
    //     b1 = utils::binarySearch(hashesMin[n], hashEnts1, b1, f1);
    // else if (hashesMin[i] > hashesMin[n])
    //     b2 = utils::binarySearch(hashesMin[i], hashEnts2, b2, f2);
    // if (hashesMax[i] > hashesMax[n])
    //     f1 = utils::binarySearch(hashesMax[n] + 1, hashEnts1, b1, f1);
    // else if (hashesMax[i] < hashesMax[n])
    //     f2 = utils::binarySearch(hashesMax[i] + 1, hashEnts2, b2, f2);
    // int add = f1 - b1 + f2 - b2;
    // while (b1 != f1 && b2 != f2)
    // {
    //     if (hashEnts1[b1] == hashEnts2[b2])
    //     {
    //         result += hashId1[b1 + 1] - hashId1[b1] + hashId2[b2 + 1] - hashId2[b2];
    //         b1++;
    //         b2++;
    //     }
    //     else if (hashEnts1[b1] < hashEnts2[b2])
    //         b1++;
    //     else
    //         b2++;
    // }

    // times+=1;
    // if(i1>n1){
    //     i=n1;
    //     n=i1;
    // }
    int l1 = degrees[i];
    int l2 = degrees[n];
    // maxs+=G[i].size()+G[n].size();
    // int cmprs=0;
    if (l1 != 0 && l2 != 0)
    {
        int m1 = hashesMin[i];
        int m2 = hashesMin[n];
        int M1 = hashesMax[i];
        int M2 = hashesMax[n];
#ifndef CSR
        vec &vec1 = G[i];
        vec &vec2 = G[n];
        vec &hashEnts1 = hashEntries[i];
        vec &hashEnts2 = hashEntries[n];
        vec &hashId1 = hashIds[i];
        vec &hashId2 = hashIds[n];
        int s1 = vec1.size();
        int s2 = vec2.size();
#else
        int *vec1 = ca + oa[i];
        int *vec2 = ca + oa[n];
        int *hashEnts1 = hashEntriesCA + hashOA[i];
        int *hashEnts2 = hashEntriesCA + hashOA[n];
        int *hashId1 = hashIdsCA + hashOA[i];
        int *hashId2 = hashIdsCA + hashOA[n];
        int s1 = oa[i + 1] - oa[i];
        int s2 = oa[n + 1] - oa[n];
#endif
        times++;
        mins += hashOA[i + 1] - hashOA[i];
        maxs += s1;
        actcmpr += hashOA[n + 1] - hashOA[n];
        orgcmpr += s2;
        avgs += (double)(hashOA[i + 1] - hashOA[i]) / s1;
        avgs2 += (double)(hashOA[n + 1] - hashOA[n]) / s2;
        if (s1 == 1)
        {
            if (s2 == 1 && vec1[0] == vec2[0])
            {
                result++;
                cout << "Something strange happened!" << endl
                     << flush;
            }
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[s2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, s2);
                    if (lb != s2 && vec2[lb] == v)
                    {
                        result++;
                        cout << "Something strange happened!" << endl
                             << flush;
                    }
                }
            }
        }
        else if (s2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[s1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, s1);
                if (lb != s1 && vec1[lb] == v)
                {
                    result++;
                    cout << "Something strange happened!" << endl
                         << flush;
                }
            }
        }
        else if (max(m1, m2) <= min(M1, M2))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
// auto bin_start_time = chrono::high_resolution_clock::now();
#ifdef TIMECOUNT
            std::chrono::_V2::system_clock::time_point merge_start_time;
            merge_start_time = std::chrono::high_resolution_clock::now();
#endif
            if (m1 < m2)
                b1 = utils::binarySearch(m2, hashEnts1, b1, f1);
            else if (m1 > m2)
                b2 = utils::binarySearch(m1, hashEnts2, b2, f2);
            if (M1 > M2)
                f1 = utils::binarySearch(M2 + 1, hashEnts1, b1, f1);
            else if (M1 < M2)
                f2 = utils::binarySearch(M1 + 1, hashEnts2, b2, f2);
            // bintime += chrono::high_resolution_clock::now()-bin_start_time;
            // mins+=f1-b1+f2-b2;
            result+=f1-b1+f2-b2;
            // actcmpr+=f1-b1+f2-b2;
            // cmprs+=f1-b1+f2-b2;
            // auto merge_start_time = std::chrono::high_resolution_clock::now();
            // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
            while (b1 != f1 && b2 != f2)
            {
                if (hashEnts1[b1] == hashEnts2[b2])
                {

                    // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
                    int vit1 = hashId1[b1];
                    int vit2 = hashId2[b2];
                    int ved1 = hashId1[b1 + 1];
                    int ved2 = hashId2[b2 + 1];
                    result+=ved1-vit1+ved2-vit2;
                    // actcmpr++;
                    // mins+=ved1-vit1+ved2-vit2;
                    // actcmpr+=ved1-vit1+ved2-vit2-1;
                    // cmprs+=ved1-vit1+ved2-vit2;
                    // actcmpr-=ved1-vit1+ved2-vit2;
                    b1++;
                    b2++;
                }
                else if (hashEnts1[b1] < hashEnts2[b2])
                {
                    b1++;
                    // actcmpr++;
                }
                else
                {
                    b2++;
                    // actcmpr++;
                }
            }
#ifdef TIMECOUNT
            mergetime += chrono::high_resolution_clock::now() - merge_start_time;
#endif

            // actcmpr-=f1-b1+f2-b2;
        }
    }
    return make_pair(result, result);
}

double indexGraph::calF(const pvec &nodePairs, const bool NB)
{
    tqdm bar;
    int marked = 0;
    long double CNDP = 0.0,CNDP2=0.0;
    
    for (const pair<int, int> &p : nodePairs)
    {
        int i = ids[p.first];
        int n = ids[p.second];
        unsigned long CN = 0;
        int CP = 0;
        commonNeighborNB(i,n,CN);
        if (G[i].size() == 1 || G[n].size() == 1 || hashIds[i][degrees[i]] == 1 || hashIds[n][degrees[n]] == 1)
            CP = 2;
        else if (degrees[i] != 0 && degrees[n] != 0)
        {
            int m = max(hashesMin[i], hashesMin[n]);
            int M = min(hashesMax[i], hashesMax[n]);
            if (m <= M)
                CP = calCP2(i, n, m, M).second;
        }
        CNDP = CNDP + (double)(2 * CN + 1);
        CNDP2 += (CP+1);
        bar.progress(++marked, N);
    }
    CNDP = CNDP / CNDP2;
    return CNDP;
}

double indexGraph::calF2(const pvec &nodePairs, const bool NB)
{
    tqdm bar;
    int marked = 0;
    long double CNDP = 0.0;
    for (const pair<int, int> &p : nodePairs)
    {
        int i = ids[p.first];
        int n = ids[p.second];
        unsigned long CN = 0;
        int CP = 0;
        Graph::commonNeighbor(i, n, CN);
        if (G[i].size() == 1 || G[n].size() == 1 || hashIds[i][degrees[i]] == 1 || hashIds[n][degrees[n]] == 1)
            CP = 2;
        else if (degrees[i] != 0 && degrees[n] != 0)
        {
            int m = max(hashesMin[i], hashesMin[n]);
            int M = min(hashesMax[i], hashesMax[n]);
            if (m <= M)
                CP = calCP2(i, n, m, M).first;
        }
        CNDP = CNDP + (double)(2 * CN + 1) / (CP + 1);
        bar.progress(++marked, N);
    }
    CNDP = CNDP / nodePairs.size();
    return CNDP;
}

void indexGraph::calRate()
{
    double r = 0.0;
    for (int i = 0; i < N; i++)
    {
        r = r + (double)(hashEntries[i].size() - 1) / G[i].size();
    }
    r = r / N;
    cout << "R: " << r << endl;
}