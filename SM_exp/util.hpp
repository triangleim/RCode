#ifndef _UTIL_H
#define _UTIL_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>
#include <cerrno>
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <x86intrin.h>
#include <unistd.h>
#include <sys/time.h>
#include <boost/filesystem.hpp>
#include <numeric>
#define CNNBV 3
#define SIMD_STATE 4 // 0:none, 2:scalar2x, 4:simd4x
#define SIMD_MODE 1 // 0:naive 1: filter
// #define TIMECOUNT
typedef int PackBase;
#ifdef SI64
typedef long long PackState;
#else
typedef int PackState;
#endif
const int PACK_WIDTH = sizeof(PackState) * 8;
const int PACK_SHIFT = __builtin_ctzll(PACK_WIDTH);
const int PACK_MASK = PACK_WIDTH - 1;

//const size_t PARA_DEG_M128 = sizeof(__m128i) / sizeof(PackState);
//const size_t PARA_DEG_M256 = sizeof(__m256i) / sizeof(PackState);

const size_t PACK_NODE_POOL_SIZE = 1024000000;

// const int CACHE_LINE_SIZE = sysconf (_SC_LEVEL1_DCACHE_LINESIZE); // in byte.
struct PackNode
{
    PackBase base;
    PackState state;

    PackNode() {};
    PackNode(PackBase _b, PackState _s): base(_b), state(_s) {};
};

struct UVertex
{
    int start, deg;
    UVertex(): start(-1), deg(0) {};
    UVertex(int _s, int _d): start(_s), deg(_d) {};
};

struct DVertex
{
    int out_start, out_deg;
    int in_start, in_deg;

    DVertex(): out_start(-1), out_deg(0), in_start(-1), in_deg(0) {};
};

typedef std::pair<int, int> Edge;
typedef std::vector<std::pair<int,int>> EdgeVector;

struct LabelSubgraph
{
    int v_num, e_num;
    EdgeVector edge_vec;
    std::vector<int> vertex2label;

    LabelSubgraph(int _v, int _e): v_num(_v), e_num(_e) {};
    void print() const;
};

void quit();
std::string extract_filename(const std::string full_filename);
int arg_pos(char *str, int argc, char **argv);
void align_malloc(void **memptr, size_t alignment, size_t size);
EdgeVector load_graph(const std::string path);
void save_graph(const std::string path, const EdgeVector& edge_vec);
std::vector<int> load_vertex_order(const std::string path);
void save_newid(const std::string path, std::vector<int> org2newid);
bool edge_idpair_cmp(const Edge& a, const Edge& b);
inline int binarySearch(int _key, int* _list, int _left, int _right) {
	int _mid;
	while (_left < _right) {
		_mid = (_right + _left) / 2;
		if (_list[_mid] < _key)
            _left = _mid+1;
        else
            _right = _mid;
    }
	return _left;
}
inline int binarySearch2(int _key, int* _list, int _left, int _right) {
	int _mid;
	while (_left < _right) {
		_mid = (_right + _left + 1) / 2;
		if (_list[_mid] > _key)
            _right = _mid-1;
        else
            _left = _mid;
    }
	return _right;
}
inline int Copy(int *set_a, int size_a, int *set_c)
{
    memcpy(set_c, set_a, size_a * sizeof(int));
    return size_a;
}
#endif
