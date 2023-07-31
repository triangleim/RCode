#include "hashGraph.h"

void hashGraph::processHG(const vec *OG)
{
    // tqdm bar;
    int marked = 0;
    oa = new int[N + 1];
    ca = new int[2 * M];

    for (int i = 0; i < N; i++)
    {

        int newid = ids[i];

        G[newid].reserve(OG[i].size());

        if (OG[i].empty()) //|| hashes[newid]==INT_MAX)
        {
            degrees[newid] = 0;
            hashesMin[newid] = INT_MAX;
            hashesMax[newid] = INT_MIN;
        }
        else
        {
            for (const int &adj : OG[i])
            {
                if (ids[adj] >= N)
                {
                    cout << adj << " Wrong!!!\n"
                         << flush;
                }
                G[newid].push_back(ids[adj]);
            }
            sort(G[newid].begin(), G[newid].end());
            if (hashes[G[newid].front()] == INT_MAX)
            {
                degrees[newid] = 0;
                hashesMin[newid] = INT_MAX;
                hashesMax[newid] = INT_MIN;
            }
            else
            {

                hashesMin[newid] = hashes[G[newid].front()];
                int back = utils::binarySearch(INT_MAX, G[newid], hashes);

                hashesMax[newid] = hashes[G[newid][back - 1]];
                degrees[newid] = back;
            }
        }
        ++marked;
        // cout<<"\r"<<marked<<"/"<<N<<flush;
        // bar.progress(++marked, N);
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
}
unsigned long long mins = 0;
unsigned long long maxs = 0;
double avgs = 0;
unsigned long long times = 0;
unsigned long long actcmpr = 0;
extern unsigned long long orgcmpr;
chrono::duration<double> mergetime(0);
chrono::duration<double> bintime(0);
void hashGraph::reportRatio(ostream &out = cout)
{
    cout << "Reduced ratio:" << (float)mins / maxs << endl
         << flush;
    // maxs=1;
    cout << "Reduced ratio(AVG):" << (float)times / maxs << endl
         << flush;
    cout<<maxs<<endl;
    // cout << "Reduced ratio(Actual):" << (float)actcmpr / orgcmpr << endl
    //      << flush;
    // cout << "Binary Search Time:" << bintime.count() << endl
    //      << flush;
    out << "Hash Merge Time:" << mergetime.count() << endl
        << flush;
    mins = 0;
    maxs = 0;
    // avgs = 0;
    times = 0;
    // actcmpr = 0;
    mergetime = chrono::duration<double>(0);
    // bintime = chrono::duration<double>(0);
    // cout << "Reduced ratio:" << (float)mins / maxs << endl
    //      << flush;
    // cout << "Reduced ratio(AVG):" << avgs / times << endl
    //      << flush;
    // cout << "Reduced ratio(Actual):" << (float)actcmpr / orgcmpr << endl
    //      << flush;
    // cout << "Binary Search Time:" << bintime.count() << endl
    //      << flush;
    // out << "Hash Merge Time:" << mergetime.count() << endl
    //     << flush;
    // mins = 0;
    // maxs = 0;
    // // avgs = 0;
    // // times = 0;
    // // actcmpr = 0;
    // mergetime = chrono::duration<double>(0);
    // bintime = chrono::duration<double>(0);
}
#ifdef CSR
double hashGraph::calTri(bool reportNumber)
{
    tqdm bar;
    auto start_time = chrono::high_resolution_clock::now();
    unsigned long result = 0;
    // cout<<"\rstart"<<flush;
    for (int i = 0; i != N; i++)
    {
        for (int j = oa[i]; j != oa[i + 1]; j++)
        {
            int n = ca[j];
            // if(n>i)
            commonNeighborNB(i, n, result);
        }
        bar.progress(i + 1, N);
        // cout<<"\r"<<i+1<<" "<<N<<flush;
    }
    auto end_time = chrono::high_resolution_clock::now();
    if (reportNumber)
        cout << endl
             << "Hash Total Triangles:" << result << endl
             << flush;
    chrono::duration<double> diff = end_time - start_time;
    return diff.count();
}
#else
double hashGraph::calTri(bool reportNumber)
{
    tqdm bar;
    auto start_time = chrono::high_resolution_clock::now();
    unsigned long result = 0;
    for (int i = 0; i != N; i++)
    {

        vec &adj = G[i];
        int s = adj.size();
        for (int j = 0; j != s; j++)
        {
            int n = adj[j];
            commonNeighborNB(i, n, result);
        }
        bar.progress(i + 1, N);
    }
    auto end_time = chrono::high_resolution_clock::now();
    if (reportNumber)
        cout << endl
             << "Total Triangles:" << result << endl;
    chrono::duration<double> diff = end_time - start_time;
    return diff.count();
}
#endif

#if CNV == 1
void hashGraph::commonNeighbor(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
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
    if (s2 == 1 && l2 == 0 && l1 < s1)
    {
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
        if (l1 == 1)
        {
            int v = vec1[0];
            if (l2 == 1 && v == vec2[0])
                result++;
            else if (v >= vec2[0] && v <= vec2[l2 - 1])
            {
                int lb = utils::binarySearch(v, vec2, 0, l2);
                if (lb != l2 && vec2[lb] == v)
                    result++;
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            maxs+=l1+l2;
            mins+=f2+f1-b2-b1;
            auto merge_start_time = std::chrono::high_resolution_clock::now(); 
            while (b1 != f1 && b2 != f2)
            {
                times+=1;
                if (vec1[b1] == vec2[b2])
                {
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                    b1++;
                else
                    b2++;
            }
            mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
        }
    }
}

void hashGraph::commonNeighborBS(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    vec &vec1 = G[i];
    vec &vec2 = G[n];
    int s1 = vec1.size();
    int s2 = vec2.size();
    if (s1 == 1 && l1 == 0)
    {
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
        if (l1 == 1)
        {
            int v = vec1[0];
            if (l2 == 1 && v == vec2[0])
                result++;
            else if (v >= vec2[0] && v <= vec2[l2 - 1])
            {
                int lb = utils::binarySearch(v, vec2, 0, l2);
                if (lb != l2 && vec2[lb] == v)
                    result++;
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            if (f1 - b1 < f2 - b2)
                utils::commonNeighborBS(vec1, vec2, result, b1, f1, b2, f2);
            else
                utils::commonNeighborBS(vec2, vec1, result, b2, f2, b1, f1);
        }
    }
}
void hashGraph::commonNeighborCP(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    vec &vec1 = G[i];
    vec &vec2 = G[n];
    int s1 = vec1.size();
    int s2 = vec2.size();
    if (s1 == 1 && l1 == 0)
    {
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
        if (l1 == 1)
        {
            int v = vec1[0];
            if (l2 == 1 && v == vec2[0])
                result++;
            else if (v >= vec2[0] && v <= vec2[l2 - 1])
            {
                int lb = utils::binarySearch(v, vec2, 0, l2);
                if (lb != l2 && vec2[lb] == v)
                    result++;
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                    b1++;
                else
                    b2++;
                result++;
            }
        }
    }
}

#elif CNV == 2
void hashGraph::commonNeighbor(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    vec &vec1 = G[i];
    vec &vec2 = G[n];
    if (l1 != 0 && l2 != 0)
    {
        if (l1 == 1)
        {
            int v = vec1[0];
            if (l2 == 1 && v == vec2[0])
                result++;
            else if (v >= vec2[0] && v <= vec2[l2 - 1])
            {
                int lb = utils::binarySearch(v, vec2, 0, l2);
                if (lb != l2 && vec2[lb] == v)
                    result++;
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                    b1++;
                else
                    b2++;
            }
        }
    }
}
#elif CNV == 3
void hashGraph::commonNeighbor(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    vec &vec1 = G[i];
    vec &vec2 = G[n];
    int s1 = vec1.size();
    int s2 = vec2.size();
    if (s1 == 1 && l1 == 0)
    {
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
        int v = vec2[0];
        if (v >= vec1[l1] && v <= vec1[s1 - 1])
        {
            int lb = utils::binarySearch(v, vec1, l1, s1);
            if (lb != s1 && vec1[lb] == v)
                result++;
        }
    }
    else if (l1 != 0 && l2 != 0 && max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
    {
        int b1 = 0;
        int f1 = l1;
        int b2 = 0;
        int f2 = l2;
        if (vec1[b1] < vec2[b2])
            b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
        else if (vec1[b1] > vec2[b2])
            b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
        if (vec1[f1 - 1] > vec2[f2 - 1])
            f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
        else if (vec1[f1 - 1] < vec2[f2 - 1])
            f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
        while (b1 != f1 && b2 != f2)
        {
            if (vec1[b1] == vec2[b2])
            {
                result++;
                b1++;
                b2++;
            }
            else if (vec1[b1] < vec2[b2])
                b1++;
            else
                b2++;
        }
    }
}
#endif
#ifndef CSR
#if CNNBV == 1
void hashGraph::commonNeighborNB(const vec &vec1, const vec &vec2, vec &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
        if (l1 == 1)
        {
            if (l2 == 1 && vec1[0] == vec2[0])
                result.push_back(vec1[0]);
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[l2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, l2);
                    if (lb != l2 && vec2[lb] == v)
                        result.push_back(v);
                }
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result.push_back(v);
            }
        }
        else if (m <= M) // maxMin<=minMax
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result.push_back(vec1[b1]);
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                    b1++;
                else
                    b2++;
            }
        }
    }
}
void hashGraph::commonNeighborNB(const vec &vec1, const vec &vec2, unsigned long &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
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
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (m <= M)
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                    b1++;
                else
                    b2++;
            }
        }
    }
}
#elif CNNBV == 3
void hashGraph::commonNeighborNB(const vec &vec1, const vec &vec2, vec &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
        if (l1 == 1)
        {
            if (l2 == 1 && vec1[0] == vec2[0])
                result.push_back(vec1[0]);
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[l2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, l2);
                    if (lb != l2 && vec2[lb] == v)
                        result.push_back(v);
                }
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result.push_back(v);
            }
        }
        else if (m <= M) // maxMin<=minMax
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            std::vector<int> revcands;
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result.push_back(vec1[b1]);
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                {
                    b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
                }
                else
                {
                    b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                else if (vec1[f1 - 1] == vec2[f2 - 1])
                {
                    revcands.push_back(vec1[f1 - 1]);
                    f1--;
                    f2--;
                }
                else if (vec1[f1 - 1] > vec2[f2 - 1])
                {
                    f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
                }
                else
                {
                    f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
                }
            }
            for (int i = revcands.size() - 1; i >= 0; i--)
            {
                result.push_back(revcands[i]);
            }
        }
    }
}
void hashGraph::commonNeighborNB(const vec &vec1, const vec &vec2, unsigned long &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
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
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (m <= M)
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                {
                    b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
                }
                else
                {
                    b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                else if (vec1[f1 - 1] == vec2[f2 - 1])
                {
                    result++;
                    f1--;
                    f2--;
                }
                else if (vec1[f1 - 1] > vec2[f2 - 1])
                {
                    f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
                }
                else
                {
                    f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
                }
            }
        }
    }
}
#elif CNNBV == 4
void hashGraph::commonNeighborNB(const vec &vec1, const vec &vec2, vec &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
        if (l1 == 1)
        {
            if (l2 == 1 && vec1[0] == vec2[0])
                result.push_back(vec1[0]);
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[l2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, l2);
                    if (lb != l2 && vec2[lb] == v)
                        result.push_back(v);
                }
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result.push_back(v);
            }
        }
        else if (m <= M) // maxMin<=minMax
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                {
                    if ((f1 - b1) > BINRATIO * (f2 - b2))
                        b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
                    else
                    {
                        while (vec1[b1] < vec2[b2])
                            b1++;
                    }
                }
                else
                {
                    if ((f2 - b2) > BINRATIO * (f1 - b1))
                        b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
                    else
                    {
                        while (vec1[b1] > vec2[b2])
                            b2++;
                    }
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                else if (vec1[f1 - 1] == vec2[f2 - 1])
                {
                    result++;
                    f1--;
                    f2--;
                }
                else if (vec1[f1 - 1] > vec2[f2 - 1])
                {
                    if ((f1 - b1) > BINRATIO * (f2 - b2))
                        f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
                    else
                    {
                        while (vec1[f1 - 1] > vec2[f2 - 1])
                        {
                            f1--;
                        }
                    }
                }
                else
                {
                    if ((f2 - b2) > BINRATIO * (f1 - b1))
                        f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
                    else
                    {
                        while (vec1[f1 - 1] < vec2[f2 - 1])
                        {
                            f2--;
                        }
                    }
                }
                // if (record)
                // {
                //     record = false;
                //     mins += f1 - b1 + f2 - b2;
                //     maxs += vec1[f1 - 1] - vec1[b1] + vec2[f2 - 1] - vec2[b2];
                // }
            }
        }
    }
}
#endif
#else
#if CNNBV == 1
void hashGraph::commonNeighborNB(const int *vec1, int *vec2, vec &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
        if (l1 == 1)
        {
            if (l2 == 1 && vec1[0] == vec2[0])
                result.push_back(vec1[0]);
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[l2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, l2);
                    if (lb != l2 && vec2[lb] == v)
                        result.push_back(v);
                }
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result.push_back(v);
            }
        }
        else if (m <= M) // maxMin<=minMax
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result.push_back(vec1[b1]);
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                    b1++;
                else
                    b2++;
            }
        }
    }
}
void hashGraph::commonNeighborNB(const int *vec1, int *vec2, unsigned long &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
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
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (m <= M)
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                    b1++;
                else
                    b2++;
            }
        }
    }
}
#elif CNNBV == 3
void hashGraph::commonNeighborNB(const int *vec1, int *vec2, vec &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
        if (l1 == 1)
        {
            if (l2 == 1 && vec1[0] == vec2[0])
                result.push_back(vec1[0]);
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[l2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, l2);
                    if (lb != l2 && vec2[lb] == v)
                        result.push_back(v);
                }
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result.push_back(v);
            }
        }
        else if (m <= M) // maxMin<=minMax
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            std::vector<int> revcands;
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result.push_back(vec1[b1]);
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                {
                    b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
                }
                else
                {
                    b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                else if (vec1[f1 - 1] == vec2[f2 - 1])
                {
                    revcands.push_back(vec1[f1 - 1]);
                    f1--;
                    f2--;
                }
                else if (vec1[f1 - 1] > vec2[f2 - 1])
                {
                    f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
                }
                else
                {
                    f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
                }
            }
            for (int i = revcands.size() - 1; i >= 0; i--)
            {
                result.push_back(revcands[i]);
            }
        }
    }
}
void hashGraph::commonNeighborNB(const int *vec1, int *vec2, unsigned long &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
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
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (m <= M)
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                {
                    b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
                }
                else
                {
                    b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                else if (vec1[f1 - 1] == vec2[f2 - 1])
                {
                    result++;
                    f1--;
                    f2--;
                }
                else if (vec1[f1 - 1] > vec2[f2 - 1])
                {
                    f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
                }
                else
                {
                    f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
                }
            }
        }
    }
}
#elif CNNBV == 4
void hashGraph::commonNeighborNB(const int *vec1, int *vec2, vec &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
        if (l1 == 1)
        {
            if (l2 == 1 && vec1[0] == vec2[0])
                result.push_back(vec1[0]);
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[l2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, l2);
                    if (lb != l2 && vec2[lb] == v)
                        result.push_back(v);
                }
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result.push_back(v);
            }
        }
        else if (m <= M) // maxMin<=minMax
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            std::vector<int> revcands;
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result.push_back(vec1[b1]);
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                {
                    b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
                }
                else
                {
                    b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                else if (vec1[f1 - 1] == vec2[f2 - 1])
                {
                    revcands.push_back(vec1[f1 - 1]);
                    f1--;
                    f2--;
                }
                else if (vec1[f1 - 1] > vec2[f2 - 1])
                {
                    f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
                }
                else
                {
                    f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
                }
            }
            for (int i = revcands.size() - 1; i >= 0; i--)
            {
                result.push_back(revcands[i]);
            }
        }
    }
}
void hashGraph::commonNeighborNB(const int *vec1, int *vec2, unsigned long &result, const int l1, const int l2, const int m, const int M)
{
    if (l1 != 0 && l2 != 0)
    {
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
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (m <= M)
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                {
                    b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
                }
                else
                {
                    b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                else if (vec1[f1 - 1] == vec2[f2 - 1])
                {
                    result++;
                    f1--;
                    f2--;
                }
                else if (vec1[f1 - 1] > vec2[f2 - 1])
                {
                    f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
                }
                else
                {
                    f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
                }
            }
        }
    }
}
#endif
#endif

#if CNNBV == 1
void hashGraph::commonNeighborNB(const int i, const int n, unsigned long &result)
{
    // times+=1;
    int l1 = degrees[i];
    int l2 = degrees[n];
    maxs+=G[i].size()+G[n].size();
    if (l1 != 0 && l2 != 0 && hashes[i] != INT_MAX && hashes[n] != INT_MAX)
    {
#ifndef CSR
        vec &vec1 = G[i];
        vec &vec2 = G[n];
#else
        int *vec1 = ca + oa[i];
        int *vec2 = ca + oa[n];
#endif
        if (l1 == 1)
        {
            if (l2 == 1 && vec1[0] == vec2[0])
            {
                result++;
                cout << "Something strange happened!" << endl
                     << flush;
            }
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[l2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, l2);
                    if (lb != l2 && vec2[lb] == v)
                    {
                        result++;
                        cout << "Something strange happened!" << endl
                             << flush;
                    }
                }
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                {
                    result++;
                    cout << "Something strange happened!" << endl
                         << flush;
                }
            }
        }
        else if (max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
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
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            // bintime += chrono::high_resolution_clock::now()-bin_start_time;
            
            mins+=f1-b1+f2-b2;
            // actcmpr+=f1-b1+f2-b2;
            // avgs+=(float)(f1-b1+f2-b2)/(G[i].size()+G[n].size());
            // auto merge_start_time = std::chrono::high_resolution_clock::now();
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    // if(vec1[b1]>n)
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                    b1++;
                else
                    b2++;
            }
#ifdef TIMECOUNT
            mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
            // actcmpr-=f1-b1+f2-b2;
        }
    }
}
#elif CNNBV == 3
void hashGraph::commonNeighborNB(const int i, const int n, unsigned long &result)
{
    // times+=1;
    int l1 = degrees[i];
    int l2 = degrees[n];
    // maxs+=G[i].size()+G[n].size();
    if (l1 != 0 && l2 != 0 && hashes[i] != INT_MAX && hashes[n] != INT_MAX)
    {
#ifndef CSR
        vec &vec1 = G[i];
        vec &vec2 = G[n];
#else
        int *vec1 = ca + oa[i];
        int *vec2 = ca + oa[n];
#endif
        if (l1 == 1)
        {
            if (l2 == 1 && vec1[0] == vec2[0])
            {
                result++;
                cout << "Something strange happened!" << endl
                     << flush;
            }
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[l2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, l2);
                    if (lb != l2 && vec2[lb] == v)
                    {
                        result++;
                        cout << "Something strange happened!" << endl
                             << flush;
                    }
                }
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                {
                    result++;
                    cout << "Something strange happened!" << endl
                         << flush;
                }
            }
        }
        else if (max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
        {

            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            bool record = true;
            // auto bin_start_time = chrono::high_resolution_clock::now();

            // bintime += chrono::high_resolution_clock::now()-bin_start_time;
            // mins+=f1-b1+f2-b2;
            // actcmpr+=f1-b1+f2-b2;
            // avgs+=(float)(f1-b1+f2-b2)/(G[i].size()+G[n].size());
            // auto merge_start_time = std::chrono::high_resolution_clock::now();

            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                {
                    // if ((f1 - b1) > 32 * (f2 - b2))
                    b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
                    // else
                    // {
                    //     while (vec1[b1] < vec2[b2])
                    //         b1++;
                    // }
                }
                else
                {
                    // if ((f2 - b2) > 32 * (f1 - b1))
                    b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
                    // else
                    // {
                    //     while (vec1[b1] > vec2[b2])
                    //         b2++;
                    // }
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                else if (vec1[f1 - 1] == vec2[f2 - 1])
                {
                    result++;
                    f1--;
                    f2--;
                }
                else if (vec1[f1 - 1] > vec2[f2 - 1])
                {
                    // if ((f1 - b1) > 32 * (f2 - b2))
                    f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
                    // else
                    // {
                    //     while (vec1[f1 - 1] > vec2[f2 - 1])
                    //     {
                    //         f1--;
                    //     }
                    // }
                }
                else
                {
                    // if ((f2 - b2) > 32 * (f1 - b1))
                    f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
                    // else
                    // {
                    //     while (vec1[f1 - 1] < vec2[f2 - 1])
                    //     {
                    //         f2--;
                    //     }
                    // }
                }
                // if (record)
                // {
                //     record = false;
                //     mins += f1 - b1 + f2 - b2;
                //     maxs += vec1[f1 - 1] - vec1[b1] + vec2[f2 - 1] - vec2[b2];
                // }
            }

            // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
            // actcmpr-=f1-b1+f2-b2;
        }
    }
}
#elif CNNBV == 4
void hashGraph::commonNeighborNB(const int i, const int n, unsigned long &result)
{
    // times+=1;
    int l1 = degrees[i];
    int l2 = degrees[n];
    // maxs+=G[i].size()+G[n].size();
    if (l1 != 0 && l2 != 0 && hashes[i] != INT_MAX && hashes[n] != INT_MAX)
    {
#ifndef CSR
        vec &vec1 = G[i];
        vec &vec2 = G[n];
#else
        int *vec1 = ca + oa[i];
        int *vec2 = ca + oa[n];
#endif
        if (l1 == 1)
        {
            if (l2 == 1 && vec1[0] == vec2[0])
            {
                result++;
                cout << "Something strange happened!" << endl
                     << flush;
            }
            else
            {
                int v = vec1[0];
                if (v >= vec2[0] && v <= vec2[l2 - 1])
                {
                    int lb = utils::binarySearch(v, vec2, 0, l2);
                    if (lb != l2 && vec2[lb] == v)
                    {
                        result++;
                        cout << "Something strange happened!" << endl
                             << flush;
                    }
                }
            }
        }
        else if (l2 == 1)
        {
            int v = vec2[0];
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                {
                    result++;
                    cout << "Something strange happened!" << endl
                         << flush;
                }
            }
        }
        else if (max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            bool record = true;
            // auto bin_start_time = chrono::high_resolution_clock::now();

            // bintime += chrono::high_resolution_clock::now()-bin_start_time;
            // mins+=f1-b1+f2-b2;
            // actcmpr+=f1-b1+f2-b2;
            // avgs+=(float)(f1-b1+f2-b2)/(G[i].size()+G[n].size());
            // auto merge_start_time = std::chrono::high_resolution_clock::now();

            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    result++;
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                {
                    if ((f1 - b1) > BINRATIO * (f2 - b2))
                        b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
                    else
                    {
                        while (vec1[b1] < vec2[b2])
                            b1++;
                    }
                }
                else
                {
                    if ((f2 - b2) > BINRATIO * (f1 - b1))
                        b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
                    else
                    {
                        while (vec1[b1] > vec2[b2])
                            b2++;
                    }
                }
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                else if (vec1[f1 - 1] == vec2[f2 - 1])
                {
                    result++;
                    f1--;
                    f2--;
                }
                else if (vec1[f1 - 1] > vec2[f2 - 1])
                {
                    if ((f1 - b1) > BINRATIO * (f2 - b2))
                        f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
                    else
                    {
                        while (vec1[f1 - 1] > vec2[f2 - 1])
                        {
                            f1--;
                        }
                    }
                }
                else
                {
                    if ((f2 - b2) > BINRATIO * (f1 - b1))
                        f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
                    else
                    {
                        while (vec1[f1 - 1] < vec2[f2 - 1])
                        {
                            f2--;
                        }
                    }
                }
                // if (record)
                // {
                //     record = false;
                //     mins += f1 - b1 + f2 - b2;
                //     maxs += vec1[f1 - 1] - vec1[b1] + vec2[f2 - 1] - vec2[b2];
                // }
            }

            // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
            // actcmpr-=f1-b1+f2-b2;
        }
    }
}
#elif CNNBV == 2
void hashGraph::commonNeighborNB(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    if (l1 != 0 && l2 != 0 && max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
    {
        vec &vec1 = G[i];
        vec &vec2 = G[n];
        int b1 = 0;
        int f1 = l1;
        int b2 = 0;
        int f2 = l2;
        if (vec1[b1] < vec2[b2])
            b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
        else if (vec1[b1] > vec2[b2])
            b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
        if (vec1[f1 - 1] > vec2[f2 - 1])
            f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
        else if (vec1[f1 - 1] < vec2[f2 - 1])
            f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
        while (b1 != f1 && b2 != f2)
        {
            if (vec1[b1] == vec2[b2])
            {
                result++;
                b1++;
                b2++;
            }
            else if (vec1[b1] < vec2[b2])
                b1++;
            else
                b2++;
        }
    }
}
#endif
void hashGraph::commonNeighborNBCP(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    if (l1 != 0 && l2 != 0)
    {
        vec &vec1 = G[i];
        vec &vec2 = G[n];
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
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                    b1++;
                else
                    b2++;
                result++;
            }
        }
    }
}

void hashGraph::commonNeighborNBBS(const int i, const int n, unsigned long &result)
{
    int l1 = degrees[i];
    int l2 = degrees[n];
    if (l1 != 0 && l2 != 0)
    {
        vec &vec1 = G[i];
        vec &vec2 = G[n];
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
            if (v >= vec1[0] && v <= vec1[l1 - 1])
            {
                int lb = utils::binarySearch(v, vec1, 0, l1);
                if (lb != l1 && vec1[lb] == v)
                    result++;
            }
        }
        else if (max(hashesMin[i], hashesMin[n]) <= min(hashesMax[i], hashesMax[n]))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
            if (vec1[b1] < vec2[b2])
                b1 = utils::binarySearch(vec2[b2], vec1, b1, f1);
            else if (vec1[b1] > vec2[b2])
                b2 = utils::binarySearch(vec1[b1], vec2, b2, f2);
            if (vec1[f1 - 1] > vec2[f2 - 1])
                f1 = utils::binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
            else if (vec1[f1 - 1] < vec2[f2 - 1])
                f2 = utils::binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
            if (f1 - b1 < f2 - b2)
                utils::commonNeighborBS(vec1, vec2, result, b1, f1, b2, f2);
            else
                utils::commonNeighborBS(vec2, vec1, result, b2, f2, b1, f1);
        }
    }
}

double hashGraph::calCN(const pvec &nodePairs, const bool NB)
{
    auto start_time = chrono::high_resolution_clock::now();
    unsigned long result = 0;
    if (NB)
    {
        for (const pair<int, int> &p : nodePairs)
            commonNeighborNB(ids[p.first], ids[p.second], result);
    }
    else
    {
        for (const pair<int, int> &p : nodePairs)
            commonNeighbor(ids[p.first], ids[p.second], result);
    }
    auto end_time = chrono::high_resolution_clock::now();
    cout << endl
         << "Total CNs:" << result << endl;
    chrono::duration<double> diff = end_time - start_time;
    return diff.count();
}

double hashGraph::calCNBS(const pvec &nodePairs, const bool NB)
{
    auto start_time = chrono::high_resolution_clock::now();
    unsigned long result = 0;
    if (NB)
    {
        for (const pair<int, int> &p : nodePairs)
            commonNeighborNBBS(ids[p.first], ids[p.second], result);
    }
    else
    {
        for (const pair<int, int> &p : nodePairs)
            commonNeighborBS(ids[p.first], ids[p.second], result);
    }
    auto end_time = chrono::high_resolution_clock::now();
    cout << endl
         << "Total CNs:" << result << endl;
    chrono::duration<double> diff = end_time - start_time;
    return diff.count();
}

unsigned long hashGraph::calCNCP(const pvec &nodePairs, const bool NB)
{
    unsigned long result = 0;
    if (NB)
    {
        for (const pair<int, int> &p : nodePairs)
            commonNeighborNBCP(ids[p.first], ids[p.second], result);
    }
    else
    {
        for (const pair<int, int> &p : nodePairs)
            commonNeighborCP(ids[p.first], ids[p.second], result);
    }
    return result;
}

double hashGraph::calF(const pvec &nodePairs, const bool NB)
{
    tqdm bar;
    long double CNDP = 0.0, CNDP2 = 0.0;
    int marked = 0;
    int psize = nodePairs.size();
    for (const pair<int, int> &p : nodePairs)
    {
        int i = ids[p.first];
        int n = ids[p.second];
        unsigned long CN = 0;
        int CP = 0;
        Graph::commonNeighbor(i, n, CN);
        if (G[i].size() == 1 || G[n].size() == 1 || degrees[i] == 1 || degrees[n] == 1)
            CP = 2;
        else if (degrees[i] != 0 && degrees[n] != 0)
        {
            int m = max(hashesMin[i], hashesMin[n]);
            int M = min(hashesMax[i], hashesMax[n]);
            if (m <= M)
            {
                int il = degrees[i];
                int jl = degrees[n];
                CP = lower_bound(G[i].begin(), G[i].begin() + il, G[n][jl - 1] + 1) - lower_bound(G[i].begin(), G[i].begin() + il, G[n].front()) + lower_bound(G[n].begin(), G[n].begin() + jl, G[i][il - 1] + 1) - lower_bound(G[n].begin(), G[n].begin() + jl, G[i].front());
            }
        }
        CNDP = CNDP + (double)(2 * CN + 1);
        CNDP2 += (CP + 1);
        bar.progress(++marked, psize);
    }
    CNDP = CNDP / CNDP2;
    return CNDP;
}

double hashGraph::calR(const pvec &nodePairs)
{
    tqdm bar;
    long double R = 0.0;
    int marked = 0;
    int psize = nodePairs.size();
    for (const pair<int, int> &p : nodePairs)
    {
        int i = p.first;
        int n = p.second;
        int CP1 = 0;
        int CP2 = 0;
        if (G[i].size() != 1 && G[n].size() != 1 && degrees[i] > 1 && degrees[n] > 1)
        {
            int m = max(hashesMin[i], hashesMin[n]);
            int M = min(hashesMax[i], hashesMax[n]);
            if (m <= M)
            {
                int il = degrees[i];
                int jl = degrees[n];
                CP1 = lower_bound(G[i].begin(), G[i].begin() + il, G[n][jl - 1] + 1) - lower_bound(G[i].begin(), G[i].begin() + il, G[n].front());
                CP2 = lower_bound(G[n].begin(), G[n].begin() + jl, G[i][il - 1] + 1) - lower_bound(G[n].begin(), G[n].begin() + jl, G[i].front());
            }
        }
        R += (CP1 != 0 && CP2 != 0) ? (double)min(CP1, CP2) / max(CP1, CP2) : 1;
        bar.progress(++marked, psize);
    }
    R /= psize;
    return R;
}

int hashGraph::selectPivot(const vec &P, const vec &X, const int mp, const int Mp)
{
    // return Graph::selectPivot(P,X);
    int p = P[0];
    int si = 0;
    int ps = P.size();
    for (const int &v : P)
    {
        int gvs = degrees[v];
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
            commonNeighborNB(P.data(), ca + oa[v], cn, ps, gvs, max(mp, hashes[*(ca + oa[v])]), min(Mp, hashes[*(ca + oa[v] + gvs - 1)]));
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
        int gvs = degrees[v];
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
            commonNeighborNB(P.data(), ca + oa[v], cn, ps, gvs, max(mp, hashes[*(ca + oa[v])]), min(Mp, hashes[*(ca + oa[v] + gvs - 1)]));
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
void hashGraph::BKP(vec &P, vec &X, unsigned long &result)
{
    // Graph::BKP(P,X,result);
    // return;
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
    int ud = degrees[u];
    while (vit1 != ps)
    {
        if (vit2 == ud || P[vit1] < G[u][vit2])
        {
            int v = P[vit1];
            int gvs = degrees[v];
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
                commonNeighborNB(P.data(), ca + oa[v], NP, ps, gvs, max(hashes[P[0]], hashes[*(ca + oa[v])]), min(hashes[P[ps - 1]], hashes[*(ca + oa[v] + gvs - 1)]));
                if (xs != 0)
                    commonNeighborNB(X.data(), ca + oa[v], NX, xs, gvs, max(hashes[X[0]], hashes[*(ca + oa[v])]), min(hashes[X[xs - 1]], hashes[*(ca + oa[v] + gvs - 1)]));
#ifdef TIMECOUNT
                mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
            }
            BKP(NP, NX, result);
            P.erase(P.begin() + vit1);
            X.insert(X.begin() + utils::binarySearch(v, X, 0, xs), v);
            ps--;
            xs++;
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
// #else
// int hashGraph::selectPivot(int* P, int* X, int ps, int xs,const int mp, const int Mp)
// {
//     // return Graph::selectPivot(P,X);
//     int p = P[0];
//     int si = 0;
//     for (int count = 0;count<ps;count++)
//     {
//         int& v=P[count];
//         int gvs = degrees[v];
//         if (gvs != 0)
//         {
//             unsigned long cn = 0;
// #ifdef TIMECOUNT
//             std::chrono::_V2::system_clock::time_point merge_start_time;
//             merge_start_time = std::chrono::high_resolution_clock::now();
// #endif
//             commonNeighborNB(P, ca+oa[v], cn, ps, gvs, max(mp, hashes[*(ca+oa[v])]), min(Mp, hashes[*(ca+oa[v]+gvs - 1)]));
// #ifdef TIMECOUNT
//             mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
// #endif
//             if (cn > si)
//             {
//                 si = cn;
//                 p = v;
//             }
//         }
//     }
//     for (int count=0;count<xs;count++)
//     {
//         int& v=X[count];
//         int gvs = degrees[v];
//         if (gvs != 0)
//         {
//             unsigned long cn = 0;
// #ifdef TIMECOUNT
//             std::chrono::_V2::system_clock::time_point merge_start_time;
//             merge_start_time = std::chrono::high_resolution_clock::now();
// #endif
//             commonNeighborNB(P, ca+oa[v], cn, ps, gvs, max(mp, hashes[*(ca+oa[v])]), min(Mp, hashes[*(ca+oa[v]+gvs - 1)]));
// #ifdef TIMECOUNT
//             mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
// #endif
//             if (cn > si)
//             {
//                 si = cn;
//                 p = v;
//             }
//         }
//     }
//     return p;
// }
// void hashGraph::BKP(int* P, int* X, int ps, int xs,unsigned long &result)
// {
//     // Graph::BKP(P,X,result);
//     // return;
//     if (ps==0)
//     {
//         if (xs==0)
//             result++;
//         return;
//     }
//     int u = selectPivot(P, X, ps,xs,hashes[P[0]], hashes[P[ps - 1]]);
//     int vit1 = 0;
//     int vit2 = 0;
//     int ud = degrees[u];
//     while (vit1 != ps)
//     {
//         if (vit2 == ud || P[vit1] < G[u][vit2])
//         {
//             int v = P[vit1];
//             int gvs = degrees[v];
//             int* NP;
//             int* NX;
//             int nps=0;
//             int nxs=0;
//             if (gvs != 0)
//             {
//                 NP=new int[min(ps, gvs)];
//                 NX=new int[min(xs, gvs)];
// #ifdef TIMECOUNT
//                 std::chrono::_V2::system_clock::time_point merge_start_time;
//                 merge_start_time = std::chrono::high_resolution_clock::now();
// #endif
//                 commonNeighborNB(P, ca+oa[v], NP, ps, gvs, max(hashes[P[0]], hashes[*(ca+oa[v])]), min(hashes[P[ps - 1]], hashes[*(ca+oa[v]+gvs - 1)]));
//                 if (xs != 0)
//                     commonNeighborNB(X, ca+oa[v], NX, xs, gvs, max(hashes[X[0]], hashes[*(ca+oa[v])]), min(hashes[X[xs - 1]], hashes[*(ca+oa[v]+gvs - 1)]));
// #ifdef TIMECOUNT
//                 mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
// #endif
//             }
//             BKP(NP, NX, nps,nxs,result);
//             P.erase(P.begin() + vit1);
//             X.insert(X.begin() + utils::binarySearch(v, X, 0, xs), v);
//             ps--;
//             xs++;
//         }
//         else if (P[vit1] == G[u][vit2])
//         {
//             vit1++;
//             vit2++;
//         }
//         else
//             vit2++;
//     }
// }
// #endif

double hashGraph::MC(const vec &nodes, const int *Vrank)
{
    tqdm bar;
    unsigned long result = 0;
    auto start_time = chrono::high_resolution_clock::now();
    int count = 0;
    for (const auto &n : nodes)
    {
        bar.progress(count++, nodes.size());
        int i = Vrank[n];
        int gns = degrees[n];
        vec P;
        P.reserve(gns);
        vec X;
        X.reserve(gns);
        for (const int &adj : G[n])
        {
            if (Vrank[adj] > i && P.size() < PSIZEUB)
            {
                if (hashes[adj] == INT_MAX)
                    result++;
                else
                    P.push_back(adj);
            }
            else if (hashes[adj] != INT_MAX)
                X.push_back(adj);
        }
        if (P.empty() && X.empty())
            continue;
        BKP(P, X, result);
    }
    auto end_time = chrono::high_resolution_clock::now();
    cout << endl
         << "Total MC:" << result << endl;
    chrono::duration<double> diff = end_time - start_time;
    return diff.count();
}

void hashGraph::writeHashId(const string hfilename, const string idfilename)
{
    ofstream idfile;
    idfile.open(idfilename);
    idfile << ":ID" << endl;
    for (int i = 0; i < N; i++)
        idfile << ids[i] << endl;
    ofstream nodefile;
    nodefile.open(hfilename);
    nodefile << ":ID,h:int" << endl;
    for (int i = 0; i < N; i++)
        nodefile << i << "," << hashes[i] << endl;
}