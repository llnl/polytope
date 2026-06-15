// polytope_internal.hh
//
// Put common includes for polytope here that you don't necessarily 
// want exposed in the public interface.
#ifndef POLYTOPE_INTERNAL_HH
#define POLYTOPE_INTERNAL_HH

#include <vector>
#include <map>
#include <set>
#include <iostream>

// An POLY_ASSERT macro, if one isn't already defined.

namespace polytope {

// Forward declare our helper abort method.
void internal_abort();

#ifndef NDEBUG
#define POLY_ASSERT(x) \
  if (!(x)) \
  { \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl; \
    internal_abort(); \
  }
#define POLY_ASSERT2(x, msg) \
  if (!(x)) \
  { \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
    internal_abort(); \
  }
#define POLY_BEGIN_CONTRACT_SCOPE { 
#define POLY_END_CONTRACT_SCOPE }
#else
#define POLY_ASSERT(x)
#define POLY_ASSERT2(x, msg)
#define POLY_BEGIN_CONTRACT_SCOPE if (false) {
#define POLY_END_CONTRACT_SCOPE }
#endif

#define POLY_CONTRACT_VAR(x) if (0 && &x == &x){}

// A requirement contract that is always on to check user input.
#define POLY_VERIFY(x)                                     \
  if (!(x)) {                                              \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl; \
    internal_abort();                                         \
  }
#define POLY_VERIFY2(x, msg)                                            \
  if (!(x)) {                                                           \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
    internal_abort();                                                   \
  }

#define POLY_CHECK(x)                                                   \
  if (!(x)) {                                                           \
    std::cout << "Check " << #x << " failed at " << __FILE__ << ":" << __LINE__ << std::endl; \
    internal_abort();                                                   \
  }

#define POLY_CHECK2(x, msg)                                             \
  if (!(x)) {                                                           \
    std::cout << "Check " << #x << " failed at " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
    internal_abort();                                                   \
  }

namespace internal {

//------------------------------------------------------------------------------
// Comparator to compare std::pair's by their first or second element.
//------------------------------------------------------------------------------
template<typename T1, typename T2>
struct ComparePairByFirstElement {
  bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const {
    return lhs.first < rhs.first;
  }
};

template<typename T1, typename T2>
struct ComparePairBySecondElement {
  bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const {
    return lhs.second < rhs.second;
  }
};

//------------------------------------------------------------------------------
// Edge comparator
//------------------------------------------------------------------------------
template <typename T1, typename T2>
bool
pairCompare(const std::pair<T1,T2> pair1, const std::pair<T1,T2> pair2) {
   return (pair1.first  <  pair2.first      ? true :
           pair1.first  == pair2.first and
           pair1.second <  pair2.second     ? true : false);
}

//------------------------------------------------------------------------------
// Update the map of thingies to unique indices.
//------------------------------------------------------------------------------
template<typename Key,
	 typename Comparator>
int
addKeyToMap(const Key& key, std::map<Key, int, Comparator>& key2id) {
  const typename std::map<Key, int>::const_iterator itr = key2id.find(key);
  int result;
  if (itr == key2id.end()) {
    result = key2id.size();
    key2id[key] = result;
  } else {
    result = itr->second;
  }
  return result;
}

//------------------------------------------------------------------------------
// An implementation of the map specialized to help constructing counters.
// This thing just overloads the index operator to start the count at zero
// for new key values.
//------------------------------------------------------------------------------
template<typename Key,
         typename Comparator = std::less<Key> >
class CounterMap: public std::map<Key, unsigned, Comparator> {
public:
  CounterMap(): std::map<Key, unsigned, Comparator>() {}
  virtual ~CounterMap() {}
  unsigned& operator[](const Key& key) {
    typename std::map<Key, unsigned>::iterator itr = this->find(key);
    if (itr == this->end()) {
      std::map<Key, unsigned, Comparator>::operator[](key) = 0U;
      itr = this->find(key);
    }
    POLY_ASSERT(itr != this->end());
    return itr->second;
  }
};

//------------------------------------------------------------------------------
// Hash two node indices uniquely to represent an edge.
//------------------------------------------------------------------------------
inline
std::pair<int, int>
hashEdge(const int i, const int j) {
  POLY_ASSERT(i != j);
  return i < j ? std::make_pair(i, j) : std::make_pair(j, i);
}

//------------------------------------------------------------------------------
// Return the positive index.
//------------------------------------------------------------------------------
inline
int
positiveID(const int x) {
  return x >= 0 ? x : ~x;
}

//------------------------------------------------------------------------------
// Sort a set of edges around a face so that sequential edges share nodes.
// We allow for one break in the chain (representing on unbounded surface).
// In such a situation we insert the new edge at the beginning of the chain, and
// return "true" indicating that a new edge was created.
//------------------------------------------------------------------------------
inline
bool
computeSortedFaceEdges(std::vector<std::pair<int, int> >& edges,
                       std::vector<int>& result) {
  typedef std::pair<int, int> EdgeHash;

  unsigned nedges = edges.size();
  POLY_ASSERT(nedges >= 2);

  // Invert the mapping, from nodes to edges.
  std::map<int, std::set<unsigned> > nodes2edges;
  internal::CounterMap<int> nodeUseCount;
  unsigned i;
  for (i = 0; i != nedges; ++i) {
    nodes2edges[edges[i].first].insert(i);
    nodes2edges[edges[i].second].insert(i);
    ++nodeUseCount[edges[i].first];
    ++nodeUseCount[edges[i].second];
  }

  // // BLAGO!
  // cerr << "Input edges :";
  // for (unsigned i = 0; i != edges.size(); ++i) cerr << " (" << edges[i].first << " " << edges[i].second << ")";
  // cerr << endl << "nodes2edges: " << endl;
  // for (std::map<int, std::set<EdgeHash> >::const_iterator itr = nodes2edges.begin();
  //      itr != nodes2edges.end();
  //      ++itr) {
  //   cerr << "   " << itr->first << " : ";
  //   for (std::set<EdgeHash>::const_iterator eitr = itr->second.begin();
  //        eitr != itr->second.end();
  //        ++eitr) cerr << " (" << eitr->first << " " << eitr->second << ")";
  //   cerr << endl;
  // }
  // cerr << "nodeUseCount: " << endl;
  // for (internal::CounterMap<int>::const_iterator itr = nodeUseCount.begin();
  //      itr != nodeUseCount.end();
  //      ++itr) {
  //   cerr << "   " << itr->first << " : " << itr->second << endl;
  // }
  // // BLAGO!

  // Look for any edges with one node in the set.  There can be at most
  // two such edges, representing the two ends of the chain.  We introduce a
  // new edge hooking those hanging nodes together, and off we go.
  int lastNode;
  std::vector<int> hangingNodes;
  for (i = 0; i != nedges; ++i) {
    if (nodeUseCount[edges[i].first] == 1 or
        nodeUseCount[edges[i].second] == 1) {
      POLY_ASSERT((nodeUseCount[edges[i].first] == 1 and nodeUseCount[edges[i].second] == 2) or
                  (nodeUseCount[edges[i].first] == 2 and nodeUseCount[edges[i].second] == 1));
      result.push_back(i);
      nodes2edges[edges[i].first].erase(i);
      nodes2edges[edges[i].second].erase(i);
      lastNode = (nodeUseCount[edges[i].first] == 1 ? edges[i].first : edges[i].second);
      hangingNodes.push_back(lastNode);
    }
  }
  POLY_ASSERT(result.size() == 0 or (hangingNodes.size() == 2 and result.size() == 2));

  // If needed create that new edge and put it in the set.
  if (hangingNodes.size() == 2) {
    result.insert(result.begin() + 1, edges.size());
    edges.push_back(internal::hashEdge(hangingNodes[0], hangingNodes[1]));
    ++nedges;
    POLY_ASSERT(result.size() == 3);
  }
  POLY_ASSERT(edges.size() == nedges);

  // Pick a node to start the chain.
  if (hangingNodes.size() == 2) {
    POLY_ASSERT(nodeUseCount[edges[result.back()].first] == 2 or
                nodeUseCount[edges[result.back()].second] == 2);
    lastNode = (nodeUseCount[edges[result.back()].first] == 2 ? 
                edges[result.back()].first :
                edges[result.back()].second);
  } else {
    lastNode = edges[0].first;
  }

  // Walk the remaining edges
  EdgeHash ehash;
  while (result.size() != nedges) {
    POLY_ASSERT(nodes2edges[lastNode].size() > 0);
    result.push_back(*nodes2edges[lastNode].begin());
    ehash = edges[result.back()];
    nodes2edges[ehash.first].erase(result.back());
    nodes2edges[ehash.second].erase(result.back());
    lastNode = (ehash.first == lastNode ? ehash.second : ehash.first);
  }
  
  // // BLAGO!
  // cerr << "Sorted edges : ";
  // for (i = 0; i != nedges; ++i) cerr << " (" << orderedEdges[i].first << " " << orderedEdges[i].second << ")";
  // cerr << endl;
  // // BLAGO!

  // Set the orientation for the ordered edges.
  lastNode = (edges[result[0]].first == edges[result[1]].first ? edges[result[0]].first : edges[result[0]].second);
  for (i = 1; i != nedges; ++i) {
    POLY_ASSERT(edges[result[i]].first == lastNode or edges[result[i]].second == lastNode);
    if (edges[result[i]].first == lastNode) {
      lastNode = edges[result[i]].second;
    } else {
      lastNode = edges[result[i]].first;
      result[i] = ~result[i];
    }
  }

  // That's it.
  POLY_BEGIN_CONTRACT_SCOPE;
  {
    POLY_ASSERT(edges.size() == result.size());
    for (auto i = 0u; i < edges.size(); ++i) {
      auto j = (i + 1u) % edges.size();
      auto ii = result[i];
      auto jj = result[j];
      POLY_CONTRACT_VAR(ii);
      POLY_CONTRACT_VAR(jj);
      POLY_ASSERT((ii >= 0 ? ii : ~ii) < (int)edges.size());
      POLY_ASSERT((jj >= 0 ? jj : ~jj) < (int)edges.size());
      POLY_ASSERT(((ii >= 0 and jj >= 0) and edges[ii].second == edges[jj].first) or
                  ((ii >= 0 and jj <  0) and edges[ii].second == edges[~jj].second) or
                  ((ii <  0 and jj >= 0) and edges[~ii].first == edges[jj].first) or
                  ((ii <  0 and jj <  0) and edges[~ii].first == edges[~jj].second));
    }
  }
  POLY_END_CONTRACT_SCOPE;
  return !(hangingNodes.empty());
}

} //end namespace internal

} //end namespace polytope

// Classes within the library.
#include "ErrorHandler.hh"

#endif
