#include <sstream>
#include <set>
#include "subgraph_match_index.hpp"
#include <malloc.h>
#include "../include/tqdm.h"
#include "util.hpp"
// TODO : New method implement
extern std::chrono::duration<double> mergetime;
#if CNNBV == 1
std::pair<int, int> CNRange(int *vec1, int l1, int *vec2, int l2, int *hashEntries1, int *hashIds1, int el1, int *hashEntries2, int *hashIds2, int el2, int *cn, int *cne, int *cni)
{
    int size_c = 0;
    int size_e = 0;
    if (l1 != 0 && l2 != 0)
    {
        int v1s = vec1[0];
        int v2s = vec2[0];
        int v1e = vec1[l1 - 1];
        int v2e = vec2[l2 - 1];
        int h1s = hashEntries1[0];
        int h2s = hashEntries2[0];
        int h1e = hashEntries1[el1 - 1];
        int h2e = hashEntries2[el2 - 1];
        if (l1 == 1)
        {
            if (l2 == 1 && v1s == v2s)
            {
                cn[size_c++] = v1s;
                cne[size_e] = h1s;
                cni[size_e++] = 0;
            }
            else if (v1s >= v2s && v1s <= v2e)
            {
                int lb = binarySearch(v1s, vec2, 0, l2);
                if (lb != l2 && vec2[lb] == v1s)
                {
                    cn[size_c++] = v1s;
                    cne[size_e] = h1s;
                    cni[size_e++] = 0;
                }
            }
        }
        else if (l2 == 1 && v2s >= v1s && v2s <= v1e)
        {
            int lb = binarySearch(v2s, vec1, 0, l1);
            if (lb != l1 && vec1[lb] == v2s)
            {
                cn[size_c++] = v2s;
                cne[size_e] = h2s;
                cni[size_e++] = 0;
            }
        }
        else if (std::max(h1s, h2s) <= std::min(h1e, h2e))
        {
            int b1 = 0;
            int f1 = el1;
            int b2 = 0;
            int f2 = el2;
// std::vector<int> Revcands;
// std::vector<int> Revids;
#ifdef TIMECOUNT
            std::chrono::_V2::system_clock::time_point merge_start_time;
            merge_start_time = std::chrono::high_resolution_clock::now();
#endif
            // std::vector<int> Revcns;
            if (h1s < h2s)
                b1 = binarySearch(h2s, hashEntries1, b1, f1);
            else if (h1s > h2s)
                b2 = binarySearch(h1s, hashEntries2, b2, f2);
            if (h1e > h2e)
                f1 = binarySearch(h2e + 1, hashEntries1, b1, f1);
            else if (h1e < h2e)
                f2 = binarySearch(h1e + 1, hashEntries2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (hashEntries1[b1] == hashEntries2[b2])
                {
                    int vit1 = hashIds1[b1];
                    int vit2 = hashIds2[b2];
                    int ved1 = hashIds1[b1 + 1];
                    int ved2 = hashIds2[b2 + 1];
                    int prevsize = size_c;
                    // std::vector<int> revcands;
                    while (vit1 != ved1 && vit2 != ved2)
                    {
                        if (vec1[vit1] == vec2[vit2])
                        {
                            cn[size_c++] = vec1[vit1];
                            vit1++;
                            vit2++;
                        }
                        else if (vec1[vit1] < vec2[vit2])
                            vit1++;
                        else
                            vit2++;
                    }
                    if (size_c != prevsize)
                    {
                        cne[size_e] = hashEntries1[b1];
                        cni[size_e++] = prevsize;
                    }
                    b1++;
                    b2++;
                }
                else if (hashEntries1[b1] < hashEntries2[b2])
                {
                    b1++;
                }
                else
                {

                    b2++;
                }
                // std::cout<<"Here2"<<std::flush;
            }
#ifdef TIMECOUNT
            mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
        }
        cni[size_e] = size_c;
    }
    return std::make_pair(size_c, size_e);
}
std::pair<int, int> OrgSubGraphMatchIndex::Copy(int *vec1, int l1, int *hashEntries1, int *hashIds1, int el1, int *cn, int *cne, int *cni, int label)
{
    // align_malloc((void **)&cn, 32, sizeof(int) * l1);
    // align_malloc((void **)&cne, 32, sizeof(int) * el1);
    // align_malloc((void **)&cni, 32, sizeof(int) * (el1 + 1));
    // cn = (int*)realloc(cn,sizeof(int) * l1);
    // cne = (int*)realloc(cne,sizeof(int) * el1);
    // cni = (int*)realloc(cni,sizeof(int) * (el1 + 1));
    // memcpy(cn, vec1, l1 * sizeof(int));
    // memcpy(cne, hashEntries1, el1 * sizeof(int));
    // memcpy(cni, hashIds1, (el1 + 1) * sizeof(int));
    int size_e = 0;
    int size_c = 0;

    int b1 = 0;
    int f1 = el1;
    // std::cout<<"\r"<<"Here!"<<std::flush;
#ifdef TIMECOUNT
    std::chrono::_V2::system_clock::time_point merge_start_time;
    merge_start_time = std::chrono::high_resolution_clock::now();
#endif
    while (b1 != f1)
    {
        if (hashEntries1[b1] != INT_MAX)
        {
            // std::cout<<"\r"<<"here~"<<std::flush;
            int vit1 = hashIds1[b1];

            int ved1 = hashIds1[b1 + 1];
            // std::cout<<"here~~"<<std::flush;
            int prevsize = size_c;
            while (vit1 != ved1)
            {
                if (vertex2label[vec1[vit1]] == label)
                {
                    cn[size_c++] = vec1[vit1];
                }
                vit1++;
            }
            // std::cout<<"here~~~"<<std::flush;
            if (size_c != prevsize)
            {
                cne[size_e] = hashEntries1[b1];
                cni[size_e++] = prevsize;
            }
        }
        b1++;
    }
#ifdef TIMECOUNT
    mergetime += std::chrono::high_resolution_clock::now() - merge_start_time;
#endif
    cni[size_e] = size_c;
    // std::cout<<"Here!!"<<std::flush;
    return std::make_pair(size_c, size_e);
}
#elif CNNBV == 3
std::pair<int, int> CNRange(int *vec1, int l1, int *vec2, int l2, int *hashEntries1, int *hashIds1, int el1, int *hashEntries2, int *hashIds2, int el2, int *cn, int *cne, int *cni)
{
    int size_c = 0;
    int size_e = 0;
    if (l1 != 0 && l2 != 0)
    {
        int v1s = vec1[0];
        int v2s = vec2[0];
        int v1e = vec1[l1 - 1];
        int v2e = vec2[l2 - 1];
        int h1s = hashEntries1[0];
        int h2s = hashEntries2[0];
        int h1e = hashEntries1[el1 - 1];
        int h2e = hashEntries2[el2 - 1];
        if (l1 == 1)
        {
            if (l2 == 1 && v1s == v2s)
            {
                cn[size_c++] = v1s;
                cne[size_e] = h1s;
                cni[size_e++] = 0;
            }
            else if (v1s >= v2s && v1s <= v2e)
            {
                int lb = binarySearch(v1s, vec2, 0, l2);
                if (lb != l2 && vec2[lb] == v1s)
                {
                    cn[size_c++] = v1s;
                    cne[size_e] = h1s;
                    cni[size_e++] = 0;
                }
            }
        }
        else if (l2 == 1 && v2s >= v1s && v2s <= v1e)
        {
            int lb = binarySearch(v2s, vec1, 0, l1);
            if (lb != l1 && vec1[lb] == v2s)
            {
                cn[size_c++] = v2s;
                cne[size_e] = h2s;
                cni[size_e++] = 0;
            }
        }
        else if (std::max(h1s, h2s) <= std::min(h1e, h2e))
        {
            int b1 = 0;
            int f1 = el1;
            int b2 = 0;
            int f2 = el2;
            std::vector<int> Revcands;
            std::vector<int> Revids;
            std::vector<int> Revcns;
            if (h1s < h2s)
                b1 = binarySearch(h2s, hashEntries1, b1, f1);
            else if (h1s > h2s)
                b2 = binarySearch(h1s, hashEntries2, b2, f2);
            if (h1e > h2e)
                f1 = binarySearch(h2e + 1, hashEntries1, b1, f1);
            else if (h1e < h2e)
                f2 = binarySearch(h1e + 1, hashEntries2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (hashEntries1[b1] == hashEntries2[b2])
                {
                    int vit1 = hashIds1[b1];
                    int vit2 = hashIds2[b2];
                    int ved1 = hashIds1[b1 + 1];
                    int ved2 = hashIds2[b2 + 1];
                    int prevsize = size_c;
                    std::vector<int> revcands;
                    while (vit1 != ved1 && vit2 != ved2)
                    {
                        if (vec1[vit1] == vec2[vit2])
                        {
                            cn[size_c++] = vec1[vit1];
                            vit1++;
                            vit2++;
                        }
                        else if (vec1[vit1] < vec2[vit2])
                            vit1 = binarySearch(vec2[vit2], vec1, vit1, ved1);
                        else
                            vit2 = binarySearch(vec1[vit1], vec2, vit2, ved2);
                        if (vit1 == ved1 || vit2 == ved2)
                        {
                            break;
                        }
                        else if (vec1[ved1 - 1] == vec2[ved2 - 1])
                        {
                            revcands.push_back(vec1[ved1 - 1]);
                            ved1--;
                            ved2--;
                        }
                        else if (vec1[ved1 - 1] > vec2[ved2 - 1])
                        {

                            ved1 = binarySearch(vec2[ved2 - 1] + 1, vec1, vit1, ved1);
                        }
                        else
                        {

                            ved2 = binarySearch(vec1[ved1 - 1] + 1, vec2, vit2, ved2);
                        }
                    }
                    for (int i = revcands.size() - 1; i >= 0; i--)
                    {
                        cn[size_c++] = revcands[i];
                    }
                    if (size_c != prevsize)
                    {
                        cne[size_e] = hashEntries1[b1];
                        cni[size_e++] = prevsize;
                    }
                    b1++;
                    b2++;
                }
                else if (hashEntries1[b1] < hashEntries2[b2])
                {
                    b1 = binarySearch(hashEntries2[b2], hashEntries1, b1, f1);
                }
                else
                {

                    b2 = binarySearch(hashEntries1[b1], hashEntries2, b2, f2);
                }
                // std::cout<<"Here2"<<std::flush;
                if (b1 == f1 || b2 == f2)
                {
                    break;
                }
                if (hashEntries1[f1 - 1] == hashEntries2[f2 - 1])
                {

                    // mergetime += chrono::high_resolution_clock::now()-merge_start_time;
                    int vit1 = hashIds1[f1 - 1];
                    int vit2 = hashIds2[f2 - 1];
                    int ved1 = hashIds1[f1 - 1 + 1];
                    int ved2 = hashIds2[f2 - 1 + 1];
                    std::vector<int> revcands;
                    int prevsize = Revcns.size();
                    // actcmpr++;
                    // mins+=ved1-vit1+ved2-vit2;
                    // actcmpr+=ved1-vit1+ved2-vit2-1;
                    // cmprs+=ved1-vit1+ved2-vit2;
                    while (vit1 != ved1 && vit2 != ved2)
                    {
                        if (vec1[ved1 - 1] == vec2[ved2 - 1])
                        {
                            Revcns.push_back(vec1[ved1 - 1]);
                            ved1--;
                            ved2--;
                        }
                        else if (vec1[ved1 - 1] > vec2[ved2 - 1])
                            ved1 = binarySearch(vec2[ved2 - 1] + 1, vec1, vit1, ved1);
                        else
                            ved2 = binarySearch(vec1[ved1 - 1] + 1, vec2, vit2, ved2);
                        if (vit1 == ved1 || vit2 == ved2)
                        {
                            break;
                        }
                        else if (vec1[vit1] == vec2[vit2])
                        {
                            // cn[size_c++] = vec1[vit1];
                            revcands.push_back(vec1[vit1]);
                            // actcmpr--;
                            // actcmpr++;
                            vit1++;
                            vit2++;
                        }
                        else if (vec1[vit1] < vec2[vit2])
                        {
                            vit1 = binarySearch(vec2[vit2], vec1, vit1, ved1);
                            // actcmpr++;
                        }
                        else
                        {
                            vit2 = binarySearch(vec1[vit1], vec2, vit2, ved2);
                            // actcmpr++;
                        }
                    }
                    for (int i = revcands.size() - 1; i >= 0; i--)
                    {
                        Revcns.push_back(revcands[i]);
                    }
                    if (Revcns.size() != prevsize)
                    {
                        Revcands.push_back(hashEntries1[f1 - 1]);
                        Revids.push_back(Revcns.size());
                    }
                    // actcmpr-=ved1-vit1+ved2-vit2;
                    f1--;
                    f2--;
                }
                else if (hashEntries1[f1 - 1] > hashEntries2[f2 - 1])
                {

                    f1 = binarySearch(hashEntries2[f2 - 1] + 1, hashEntries1, b1, f1);
                }
                else
                {

                    f2 = binarySearch(hashEntries1[f1 - 1] + 1, hashEntries2, b2, f2);
                }
            }
            for (int i = Revcns.size() - 1; i >= 0; i--)
            {
                cn[size_c++] = Revcns[i];
            }
            for (int i = Revcands.size() - 1; i >= 0; i--)
            {
                cne[size_e] = Revcands[i];
                cni[size_e++] = size_c - Revids[i];
            }
        }
        cni[size_e] = size_c;
    }
    return std::make_pair(size_c, size_e);
}
std::pair<int, int> OrgSubGraphMatchIndex::Copy(int *vec1, int l1, int *hashEntries1, int *hashIds1, int el1, int *cn, int *cne, int *cni, int label)
{
    // align_malloc((void **)&cn, 32, sizeof(int) * l1);
    // align_malloc((void **)&cne, 32, sizeof(int) * el1);
    // align_malloc((void **)&cni, 32, sizeof(int) * (el1 + 1));
    // cn = (int*)realloc(cn,sizeof(int) * l1);
    // cne = (int*)realloc(cne,sizeof(int) * el1);
    // cni = (int*)realloc(cni,sizeof(int) * (el1 + 1));
    // memcpy(cn, vec1, l1 * sizeof(int));
    // memcpy(cne, hashEntries1, el1 * sizeof(int));
    // memcpy(cni, hashIds1, (el1 + 1) * sizeof(int));
    int size_e = 0;
    int size_c = 0;

    int b1 = 0;
    int f1 = el1;
    // std::cout<<"\r"<<"Here!"<<std::flush;
    while (b1 != f1)
    {
        if (hashEntries1[b1] != INT_MAX)
        {
            // std::cout<<"\r"<<"here~"<<std::flush;
            int vit1 = hashIds1[b1];

            int ved1 = hashIds1[b1 + 1];
            // std::cout<<"here~~"<<std::flush;
            int prevsize = size_c;
            while (vit1 != ved1)
            {
                if (vertex2label[vec1[vit1]] == label)
                {
                    cn[size_c++] = vec1[vit1];
                }
                vit1++;
            }
            // std::cout<<"here~~~"<<std::flush;
            if (size_c != prevsize)
            {
                cne[size_e] = hashEntries1[b1];
                cni[size_e++] = prevsize;
            }
        }
        b1++;
    }
    cni[size_e] = size_c;
    // std::cout<<"Here!!"<<std::flush;
    return std::make_pair(size_c, size_e);
}
#endif

OrgSubGraphMatchIndex::OrgSubGraphMatchIndex()
{
    v_num = 0;
    e_num = 0;
    l_num = 0;
}

OrgSubGraphMatchIndex::~OrgSubGraphMatchIndex()
{
    if (pool_edges != NULL)
        free(pool_edges);
    if (pool_vertices != NULL)
        free(pool_vertices);
    if (pool_entries != NULL)
        free(pool_entries);
    if (pool_ids != NULL)
        free(pool_ids);
    for (int i = 0; i < 2; ++i)
    {
        if (temp_join_res[i] != NULL)
            free(temp_join_res[i]);
        if (temp_join_entries[i] != NULL)
            free(temp_join_entries[i]);
        if (temp_join_ids[i] != NULL)
            free(temp_join_ids[i]);
    }
}

void OrgSubGraphMatchIndex::build(const EdgeVectorh &_e_v, const std::vector<int> _v_l)
{
    vertex2label = _v_l;

    for (auto &e : _e_v)
    {
        v_num = std::max(v_num, e[0]);
        v_num = std::max(v_num, e[1]);
    }
    v_num++;
    e_num = (long long)_e_v.size();
    for (auto l : vertex2label)
        l_num = std::max(l_num, l);
    l_num++;

    align_malloc((void **)&pool_edges, 32, sizeof(int) * e_num);
    align_malloc((void **)&pool_vertices, 32, sizeof(int) * v_num);
    align_malloc((void **)&pool_entries, 32, sizeof(int) * e_num);
    align_malloc((void **)&pool_ids, 32, sizeof(int) * e_num);

    for (int i = 0; i < 2; ++i)
    {
        align_malloc((void **)&temp_join_res[i], 32, sizeof(int) * v_num);
        align_malloc((void **)&temp_join_entries[i], 32, sizeof(int) * v_num);
        align_malloc((void **)&temp_join_ids[i], 32, sizeof(int) * v_num);
    }

    graph.resize(v_num);
    graphh.resize(v_num);
    int cur_node_idx = 0;
    int cur_entries_idx = 0;
    int prev_u = -1;
    int prev_h = INT_MIN;
    for (auto &e : _e_v)
    {
        if (e[0] != prev_u)
        {
            if (prev_u != -1 && prev_h != INT_MAX)
            {
                pool_entries[cur_entries_idx] = INT_MAX;
                pool_ids[cur_entries_idx++] = graph[prev_u].deg;
            }
            prev_u = e[0];
            prev_h = INT_MIN;
            graph[e[0]].start = cur_node_idx;
            graphh[e[0]].start = cur_entries_idx;
        }
        if (e[2] != prev_h)
        {
            prev_h = e[2];
            if (e[2] != INT_MAX)
                graphh[e[0]].deg++; // hash数
            pool_entries[cur_entries_idx] = e[2];
            pool_ids[cur_entries_idx++] = graph[e[0]].deg;
        }
        if (e[2] != INT_MAX)
            graph[e[0]].deg++;
        pool_edges[cur_node_idx++] = e[1];
    }

    EdgeVector label2vertices;
    label2vertices.reserve(v_num);
    for (int i = 0; i < v_num; ++i)
        label2vertices.push_back(Edge(vertex2label[i], i));
    std::sort(label2vertices.begin(), label2vertices.end(), edge_idpair_cmp);
    labels.resize(l_num);
    for (auto l : vertex2label)
        labels[l].deg++;
    for (int i = 0; i < v_num; ++i)
        pool_vertices[i] = label2vertices[i].second;
    labels[0].start = 0;
    for (int i = 1; i < l_num; ++i)
        labels[i].start = labels[i - 1].start + labels[i - 1].deg;

    printf("v_num=%d l_num=%d e_num=%lld\n", v_num, l_num, e_num);
}

std::vector<std::vector<int>> OrgSubGraphMatchIndex::subgraph_matching(const LabelSubgraph &q)
{
    std::vector<std::vector<int>> res;

    // build query graph's edge index.
    std::set<Edge> q_edge_set;
    std::vector<std::vector<int>> q_edge_list(q.v_num);
    for (auto e : q.edge_vec)
    {
        if (e.first > e.second)
            std::swap(e.first, e.second);
        q_edge_set.insert(e);
        q_edge_list[e.first].push_back(e.second);
        q_edge_list[e.second].push_back(e.first);
    }
    // sort query graph's edge lists (q_edge_list).
    for (auto &el : q_edge_list)
    {
        std::sort(el.begin(), el.end(),
                  [&](const int &a, const int &b) -> bool
                  {
                      if (labels[q.vertex2label[a]].deg == labels[q.vertex2label[b]].deg)
                          return a < b;
                      return labels[q.vertex2label[a]].deg < labels[q.vertex2label[b]].deg;
                  });
    }

    // // q.print();

    // decide the join order.
    std::vector<int> join_order;
    join_order.reserve(q.v_num);
    int first_vertex = 0;
    for (int i = 1; i < q.v_num; ++i)
    {
        if (labels[q.vertex2label[i]].deg < labels[q.vertex2label[first_vertex]].deg)
            first_vertex = i;
    }
    std::vector<int> added_vertices(q.v_num, 0);
    gen_join_order(first_vertex, q_edge_list, added_vertices, join_order);

    if ((int)join_order.size() < q.v_num)
    {
        printf("query graph is not conntected!\n");
        printf("join_order.size()=%lu\n", join_order.size());
        return res;
    }

    // printf("join_order: ");
    // for (auto u : join_order) printf("%d ", u); printf("\n");

    int cur_label = q.vertex2label[join_order[0]];
    for (int i = 0; i < labels[cur_label].deg; ++i)
    {
        int u = pool_vertices[labels[cur_label].start + i];
        res.push_back(std::vector<int>(1, u));
    }

    for (int i = 1; i < q.v_num; ++i)
    {
        int sg_u = join_order[i];
        int sg_u_label = q.vertex2label[sg_u];

        std::vector<std::vector<int>> next_res;
        std::vector<int> cur_join_idx;
        for (int j = 0; j < i; ++j)
        {
            Edge e(join_order[j], sg_u);
            if (e.first > e.second)
                std::swap(e.first, e.second);
            if (q_edge_set.find(e) != q_edge_set.end())
                cur_join_idx.push_back(j);
        }
        assert(cur_join_idx.size() > 0);
        // printf("sg_u=%d sg_u_label=%d cur_join_idx: ", sg_u, sg_u_label);
        // for (auto idx : cur_join_idx) printf("%d ", idx); printf("\n");

        // int *cand_list = pool_vertices + labels[sg_u_label].start;
        // int cand_size = labels[sg_u_label].deg;
        int counter = 0;
        tqdm bar;
        for (const auto &rec : res)
        {
            bar.progress(++counter, res.size());
            int g_v = rec[cur_join_idx[0]];
            // printf("cand_list: ");
            // for (int j = 0; j < cand_size; ++j) printf("%d ", cand_list[j]); printf("\n");
            // printf("edge_list: ");
            // for (int j = 0; j < graph[g_v].deg; ++j) printf("%d ", pool_edges[graph[g_v].start + j]); printf("\n");
            std::pair<int, int> sce = Copy(pool_edges + graph[g_v].start, graph[g_v].deg,
                                           pool_entries + graphh[g_v].start, pool_ids + graphh[g_v].start, graphh[g_v].deg,
                                           temp_join_res[0], temp_join_entries[0], temp_join_ids[0], sg_u_label);
            temp_join_size[0] = sce.first;          // graph[g_v].deg
            temp_join_entries_size[0] = sce.second; // graphh[g_v].deg
            // temp_join_size[0] = CNRange(cand_list, cand_size,
            //         pool_edges + graph[g_v].start, graph[g_v].deg, temp_join_res[0], vertex2hash);
            // temp_join_size[0] = intersect(cand_list, cand_size,
            //         pool_edges + graph[g_v].start, graph[g_v].deg, temp_join_res[0]);
            // printf("g_v=%d temp_join_size[0]=%d\n", g_v, temp_join_size[0]);
            for (size_t j = 1; j < cur_join_idx.size(); ++j)
            {
                g_v = rec[cur_join_idx[j]];
                sce = CNRange(temp_join_res[(j & 1) ^ 1], temp_join_size[(j & 1) ^ 1], pool_edges + graph[g_v].start, graph[g_v].deg,
                              temp_join_entries[(j & 1) ^ 1], temp_join_ids[(j & 1) ^ 1], temp_join_entries_size[(j & 1) ^ 1],
                              pool_entries + graphh[g_v].start, pool_ids + graphh[g_v].start, graphh[g_v].deg,
                              temp_join_res[j & 1], temp_join_entries[j & 1], temp_join_ids[j & 1]);
                temp_join_size[j & 1] = sce.first;
                temp_join_entries_size[j & 1] = sce.second;
                // temp_join_size[j & 1] = CNRange(temp_join_res[j ^ 1], temp_join_size[j ^ 1],
                //         pool_edges + graph[g_v].start, graph[g_v].deg, temp_join_res[j & 1], vertex2hash);
                // temp_join_size[j & 1] = intersect(temp_join_res[j ^ 1], temp_join_size[j ^ 1],
                //         pool_edges + graph[g_v].start, graph[g_v].deg, temp_join_res[j & 1]);
            }
            size_t last_join_idx = (cur_join_idx.size() - 1) & 1;
            for (int j = 0; j < temp_join_size[last_join_idx]; ++j)
            {
                int join_v = temp_join_res[last_join_idx][j];
                bool is_legal = true;
                for (auto joined_u : rec)
                    if (joined_u == join_v)
                    {
                        is_legal = false;
                        break;
                    }
                if (is_legal)
                {
                    auto new_rec = rec;
                    new_rec.push_back(join_v);
                    next_res.push_back(new_rec);
                }
            }
        }
        // res = next_res;
        res = std::move(next_res);
    }

    // rearrange result records to the origin order.
    for (auto &rec : res)
    {
        auto new_rec = rec;
        for (int i = 0; i < q.v_num; ++i)
            new_rec[join_order[i]] = rec[i];
        rec = std::move(new_rec);
        // printf("rec: ");
        // for (auto u : rec) printf("%d ", u); printf("\n");
    }
    std::cout << "\r" << res.size() << std::endl;
    return res;
}

void OrgSubGraphMatchIndex::gen_join_order(int u,
                                           std::vector<std::vector<int>> &q_edge_list,
                                           std::vector<int> &visited, std::vector<int> &join_order)
{
    visited[u] = 1;
    join_order.push_back(u);
    auto &nbrs = q_edge_list[u];
    for (auto v : nbrs)
        if (visited[v] == 0)
            gen_join_order(v, q_edge_list, visited, join_order);
}

void OrgSubGraphMatchIndex::reportRatio(std::ostream &out = std::cout)
{
    // cout << "Reduced ratio:" << (float)mins / maxs << endl
    //      << flush;
    // cout << "Reduced ratio(AVG):" << avgs / times << endl
    //      << flush;
    // cout << "Reduced ratio(Actual):" << (float)actcmpr / orgcmpr << endl
    //      << flush;
    // cout << "Binary Search Time:" << bintime.count() << endl
    //      << flush;
    out << "Index Merge Time:" << mergetime.count() << std::endl
        << std::flush;

    // avgs = 0;
    // times = 0;
    // actcmpr = 0;
    mergetime = std::chrono::duration<double>(0);
    // bintime = chrono::duration<double>(0);
}