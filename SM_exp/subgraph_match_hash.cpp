#include <sstream>
#include <set>
#include <vector>
#include "subgraph_match_hash.hpp"
#include "../include/tqdm.h"
#include "util.hpp"
extern std::chrono::duration<double> mergetime;
int intersect1(int *set_a, int size_a, int *set_b, int size_b, int *set_c)
{
    int i = 0, j = 0, size_c = 0;
    while (i != size_a && j != size_b)
    {
        if (set_a[i] == set_b[j])
        {
            set_c[size_c++] = set_a[i];
            i++;
            j++;
        }
        else if (set_a[i] < set_b[j])
        {
            i++;
        }
        else
        {
            j++;
        }
    }

    return size_c;
}
// TODO : New method implement
#if CNNBV == 1
int CNRange(int *vec1, int l1, int *vec2, int l2, int *cn, std::vector<int> &hashes)
{
    int size_c = 0;
    if (l1 != 0 && l2 != 0)
    {
        int v1s = vec1[0];
        int v2s = vec2[0];
        int v1e = vec1[l1 - 1];
        int v2e = vec2[l2 - 1];
        if (l1 == 1)
        {
            if (l2 == 1 && v1s == v2s)
                cn[size_c++] = v1s;
            else if (v1s >= v2s && v1s <= v2e)
            {
                int lb = binarySearch(v1s, vec2, 0, l2);
                if (lb != l2 && vec2[lb] == v1s)
                    cn[size_c++] = v1s;
            }
        }
        else if (l2 == 1 && v2s >= v1s && v2s <= v1e)
        {
            int lb = binarySearch(v2s, vec1, 0, l1);
            if (lb != l1 && vec1[lb] == v2s)
                cn[size_c++] = v2s;
        }
        else if (std::max(hashes[vec1[0]], hashes[vec2[0]]) <= std::min(hashes[vec1[l1 - 1]], hashes[vec2[l2 - 1]]))
        {
            int b1 = 0;
            int f1 = l1;
            int b2 = 0;
            int f2 = l2;
#ifdef TIMECOUNT
            std::chrono::_V2::system_clock::time_point merge_start_time;
            merge_start_time = std::chrono::high_resolution_clock::now();
#endif
            if (vec1[0] < vec2[0])
                b1 = binarySearch(vec2[0], vec1, b1, f1);
            else if (vec2[0] > vec1[0])
                b2 = binarySearch(vec1[0], vec2, b2, f2);
            if (vec1[l1 - 1] > vec2[l2 - 1])
                f1 = binarySearch(vec2[l2 - 1] + 1, vec1, b1, f1);
            else if (vec1[l1 - 1] < vec2[l2 - 1])
                f2 = binarySearch(vec1[l1 - 1] + 1, vec2, b2, f2);
            while (b1 != f1 && b2 != f2)
            {
                if (vec1[b1] == vec2[b2])
                {
                    cn[size_c++] = vec1[b1];
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
        }
    }
    return size_c;
}
#elif CNNBV == 3
int CNRange(int *vec1, int l1, int *vec2, int l2, int *cn, std::vector<int> &hashes)
{
    int size_c = 0;
    if (l1 != 0 && l2 != 0)
    {
        int v1s = vec1[0];
        int v2s = vec2[0];
        int v1e = vec1[l1 - 1];
        int v2e = vec2[l2 - 1];
        if (l1 == 1)
        {
            if (l2 == 1 && v1s == v2s)
                cn[size_c++] = v1s;
            else if (v1s >= v2s && v1s <= v2e)
            {
                int lb = binarySearch(v1s, vec2, 0, l2);
                if (lb != l2 && vec2[lb] == v1s)
                    cn[size_c++] = v1s;
            }
        }
        else if (l2 == 1 && v2s >= v1s && v2s <= v1e)
        {
            int lb = binarySearch(v2s, vec1, 0, l1);
            if (lb != l1 && vec1[lb] == v2s)
                cn[size_c++] = v2s;
        }
        else if (std::max(hashes[vec1[0]], hashes[vec2[0]]) <= std::min(hashes[vec1[l1 - 1]], hashes[vec2[l2 - 1]]))
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
                    cn[size_c++] = vec1[b1];
                    b1++;
                    b2++;
                }
                else if (vec1[b1] < vec2[b2])
                {
                    b1 = binarySearch(vec2[b2], vec1, b1, f1);
                }
                else
                {
                    b2 = binarySearch(vec1[b1], vec2, b2, f2);
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
                    f1 = binarySearch(vec2[f2 - 1] + 1, vec1, b1, f1);
                }
                else
                {
                    f2 = binarySearch(vec1[f1 - 1] + 1, vec2, b2, f2);
                }
            }
            for (int i = revcands.size() - 1; i >= 0; i--)
            {
                cn[size_c++] = revcands[i];
            }
        }
    }
    return size_c;
}
#endif
OrgSubGraphMatchHash::OrgSubGraphMatchHash()
{
    v_num = 0;
    e_num = 0;
    l_num = 0;
}

OrgSubGraphMatchHash::~OrgSubGraphMatchHash()
{
    if (pool_edges != NULL)
        free(pool_edges);
    if (pool_vertices != NULL)
        free(pool_vertices);
    for (int i = 0; i < 2; ++i)
        if (temp_join_res[i] != NULL)
            free(temp_join_res[i]);
}

void OrgSubGraphMatchHash::build(const EdgeVector &_e_v, const std::vector<int> _v_l, const std::vector<int> _v_h)
{
    // std::cout << 'here!'<< std::endl<<std::flush;
    vertex2label = _v_l;
    vertex2hash = _v_h;

    edge_vec.reserve(_e_v.size());
    for (auto &e : _e_v)
        if (e.first != e.second)
            edge_vec.push_back(e);
    // edge_vec = _e_v;

    std::sort(edge_vec.begin(), edge_vec.end(), edge_idpair_cmp);
    edge_vec.erase(std::unique(edge_vec.begin(), edge_vec.end()), edge_vec.end());

    for (auto &e : edge_vec)
    {
        v_num = std::max(v_num, e.first);
        v_num = std::max(v_num, e.second);
    }
    v_num++;
    e_num = (long long)edge_vec.size();
    for (auto l : vertex2label)
        l_num = std::max(l_num, l);
    l_num++;
    // std::cout << 'here!!!'<< std::endl<<std::flush;
    align_malloc((void **)&pool_edges, 32, sizeof(int) * e_num);
    align_malloc((void **)&pool_vertices, 32, sizeof(int) * v_num);

    for (int i = 0; i < 2; ++i)
        align_malloc((void **)&temp_join_res[i], 32, sizeof(int) * v_num);
    // std::cout << 'here!!'<< std::endl<<std::flush;
    graph.resize(v_num);
    int cur_node_idx = 0;
    int prev_u = -1;
    // std::cout << 'here!!!!'<< std::endl<<std::flush;
    auto count = 0;
    for (auto &e : edge_vec)
    {
        if (e.first != prev_u)
        {
            prev_u = e.first;
            graph[e.first].start = cur_node_idx;
        }

        if (vertex2hash[e.second] != INT_MAX)
        {
            // std::cout << "\r"<<e.first<<" "<<e.second<<std::flush;
            graph[e.first].deg++;
        }
        else
        {
            std::cout << "strange" << std::endl
                      << std::flush;
        }

        pool_edges[cur_node_idx++] = e.second;
        // std::cout << "hhh~"<<std::flush;
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

std::vector<std::vector<int>> OrgSubGraphMatchHash::subgraph_matching(const LabelSubgraph &q)
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

    // q.print();

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
        int u = pool_vertices[labels[cur_label].start + i]; // 就是i
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

        int *cand_list = pool_vertices + labels[sg_u_label].start;
        int cand_size = labels[sg_u_label].deg;
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
            temp_join_size[0] = CNRange(cand_list, cand_size,
                                        pool_edges + graph[g_v].start, graph[g_v].deg, temp_join_res[0], vertex2hash);
            // temp_join_size[0] = intersect1(cand_list, cand_size,
            //         pool_edges + graph[g_v].start, graph[g_v].deg, temp_join_res[0]);
            // printf("g_v=%d temp_join_size[0]=%d\n", g_v, temp_join_size[0]);
            // temp_join_size[0] = Copy(pool_edges + graph[g_v].start, graph[g_v].deg, temp_join_res[0]);
            for (size_t j = 1; j < cur_join_idx.size(); ++j)
            {
                g_v = rec[cur_join_idx[j]];
                temp_join_size[j & 1] = CNRange(temp_join_res[(j & 1) ^ 1], temp_join_size[(j & 1) ^ 1],
                                                pool_edges + graph[g_v].start, graph[g_v].deg, temp_join_res[j & 1], vertex2hash);
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
    std::cout << res.size() << std::endl;
    return res;
}

void OrgSubGraphMatchHash::gen_join_order(int u,
                                          std::vector<std::vector<int>> &q_edge_list,
                                          std::vector<int> &visited, std::vector<int> &join_order)
{
    visited[u] = 1;
    join_order.push_back(u); // DFS顺序
    auto &nbrs = q_edge_list[u];
    for (auto v : nbrs)
        if (visited[v] == 0)
            gen_join_order(v, q_edge_list, visited, join_order);
}

void OrgSubGraphMatchHash::reportRatio(std::ostream &out = std::cout)
{
    // cout << "Reduced ratio:" << (float)mins / maxs << endl
    //      << flush;
    // cout << "Reduced ratio(AVG):" << avgs / times << endl
    //      << flush;
    // cout << "Reduced ratio(Actual):" << (float)actcmpr / orgcmpr << endl
    //      << flush;
    // cout << "Binary Search Time:" << bintime.count() << endl
    //      << flush;
    out << "Hash Merge Time:" << mergetime.count() << std::endl
        << std::flush;
    // avgs = 0;
    // times = 0;
    // actcmpr = 0;
    mergetime = std::chrono::duration<double>(0);
    // bintime = chrono::duration<double>(0);
}