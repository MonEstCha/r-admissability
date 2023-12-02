#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <utility>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <time.h>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/config/user.hpp>
//#include "matplotlibcpp.h"
#include "Headers.hpp"
#include "FilesOps.hpp"
#include "FlagParser.hpp"
#include "ReadTxt.hpp"

// compile: c++ -I /opt/homebrew/Cellar/boost/1.78.0_1/include SimAnnealing_v2.cpp ReadTxt.cpp FilesOps.cpp FlagParser.cpp -o SimAnnealing_v2 -Ofast -std=c++14
// including matplotlibcpp (not working):  c++ -I /opt/homebrew/Cellar/boost/1.78.0_1/include -I /Applications/Spyder.app/Contents/Resources/lib/python3.9/numpy/core/include -I /System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 SimAnnealing.cpp ReadTxt.cpp FilesOps.cpp FlagParser.cpp -o SimAnnealing -Ofast -std=c++14
// execute: ./SimAnnealing_v2 --in=yeast.csv --rad=1 --heur=wReachLeft

using namespace boost;

const char* usage="\nUSAGE:\n\n\t./appr_tw -i <filename>  [-p] [-f <vertex_id>] [-o] [-e] [-h]\n\n\twhere\n\t-i <filename>:\t\tfile <filename> should be in the *.tsv format\n\t\t\t\t(every line contains two numbers describing an edge),\n\t-o:\t\t\tprint the order found (default: no),\n\t-e: \t\t\tshow progrEss (default: no),\n\t-t:\t\t\tprint times used for computations (default: no),\n\t-c: \t\t\tperform a test,\n\t-d <level>: \t\tprint details of level (1-5) (default: 0),\n\t-w: \t\t\tprint computed width\n\t-r: \t\t\tchoose the next vertex randomly \n\t-h: \t\t\tprint this help\n\n";

/*************** Globals and defintions for wReachLeft heuristic ***************/

vector<int> _wreach_szs;
vector<int> _deg;

struct Vert {
  int id;
  // overload "<" breaking ties by degrees
  bool operator<(const Vert& oth) const {
    if (_wreach_szs[id] != _wreach_szs[oth.id]) { return _wreach_szs[id] > _wreach_szs[oth.id]; }
    if (_deg[id] != _deg[oth.id]) { return _deg[id] > _deg[oth.id]; }
    return id < oth.id;
  }
};

/*************************** Some convenient functions *************************/

/**
 * Function checks if fullString ends with ending
 * @param fullString, String to check for ending
 * @param ending, ending to look for in fullString
 * @return true if check positive, false if negative
 */
bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

void Err() {
  cerr<<"Usage: ./SimAnnealing --in=graph.txtg --rad=radius [--o=output.txt]"<<endl;
  cerr<<"--h for help\n";
  exit(1);
}

template < typename Graph, typename VertexNamePropertyMap >
void printGraph(Graph g){
	using v_it = typename graph_traits<Graph>::vertex_iterator;
	cout << "vertices" << endl;
	pair<v_it,v_it> vs = vertices(g);
	copy(vs.st, vs.second,ostream_iterator<typename graph_traits<Graph>::vertex_descriptor>{cout, "\n"});

	using e_it = typename graph_traits<Graph>::edge_iterator;
	cout << "edges" << endl;
	pair<e_it,e_it> es = boost::edges(g);
	copy(es.first, es.second,ostream_iterator<typename graph_traits<Graph>::edge_descriptor>{cout, "\n"});
}

template < typename Graph, typename VertexNameMap >
int count_adj_vertices(typename graph_traits< Graph >::vertex_descriptor u, const Graph& g,
    VertexNameMap name_map)
{
	int ct = 0;
    typename graph_traits< Graph >::adjacency_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = adjacent_vertices(u, g); vi != vi_end; ++vi){
    	//cout << "in line 121 "<< ct << endl;
    	ct++;
    }
    return ct;
}

template < typename Graph >
void printAdjacentVerts(const Graph& g, typename graph_traits<Graph>::vertex_descriptor v){
	cout << "neighbors of " << v << ": {";
	// print
	typename graph_traits<Graph>::adjacency_iterator vit, vend;
	tie(vit, vend) = adjacent_vertices(v, g);
	copy(vit, vend, ostream_iterator<typename graph_traits<Graph>::vertex_descriptor>{cout, ", "});
	// count
	typename graph_traits<Graph>::adjacency_iterator vi, vi_end, next;
	int ct = 0;
	for (tie(vi, vi_end) = adjacent_vertices(v, g); vi != vi_end; ++vi){
		ct++;
	}
	cout << "}, ct: " << ct << endl;
}

template <typename descVec, typename VertexNameMap, typename Graph>
void printVec(const Graph& g, descVec vec, typename graph_traits< Graph >::vertex_descriptor root, VertexNameMap name_map, string desc){
	cout << desc << " of " << get(name_map, root) << ": {";
	for (typename graph_traits<Graph>::vertex_descriptor vD : vec) {
		int v = get(name_map, vD);
		cout << v << ", ";
	}
	cout << "}" << endl;
}

/**
 * Helper to calculate (symm.) difference between two vector data structures
 * @param v1, v2, sort; already ordered vectors, if sort is set to false
 * @return difference of the vectors
 */
template <typename descVec>
descVec getDifferenceOfVecs(descVec v1, descVec v2){

	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	descVec diff;
	set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(diff));
	return diff;
}

template <typename descVec>
descVec getIntersecOfVecs(descVec v1, descVec v2){

	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());

	descVec intersec;
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(intersec));
	return intersec;
}

template <typename descVec>
descVec getUnionOfVecs(descVec v1, descVec v2){

	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());

	descVec un;
	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(un));
	return un;
}

/*
 * Function calculates the probability of taking a new solution even if it is worse that the current one
 */
float getProb(int t,int wcolOld, int wcolNew, int n ){
	float wcolRatio = (float)(1 - (wcolNew - wcolOld)*0.01);
	float tempRatio = (0.5*n + (float) t) / (float) n;
	return exp(- tempRatio) * wcolRatio;
	/*float wcolRatio = (float)(wcolNew - wcolOld)*6;
	float tempRatio = 1.0 / (float) swapLim;
	return exp(- wcolRatio * tempRatio);*/
}

/*************************** Functions adapted from Nadara et al. (2019) **************************/

/**
 * Function calculates the vertices that contain the root vertex in their potential weakly R-reachable set for a given order
 * @param where_in_order, the position of the vertices in the order
 * @param phase_id, the current position to fill in the ordering
 * @return a vector containing all vertices that have root in their weakly R-reachable set
 */
template < typename Graph, typename VertexNameMap >
vector<typename graph_traits<Graph>::vertex_descriptor> ComputeSingleCluster(const Graph& graph,
                                 vector<int>& where_in_order,
                                 int R,
                                 vector<int>& is_forb,
                                 vector<int>& last_vis,
                                 vector<int>& dis,
                                 typename graph_traits<Graph>::vertex_descriptor u,
								 VertexNameMap name_map,
                                 int phase_id) {
	int root = get(name_map,u);
	//cout << "root: " <<root  << endl;

	vector<typename graph_traits<Graph>::vertex_descriptor> res;
	//if(root == 0) return res;
	if (!is_forb.empty() && is_forb[root]) { return {}; }
	last_vis[root] = phase_id;
	dis[root] = 0;
	vector<typename graph_traits<Graph>::vertex_descriptor> que{u};
	// breadth first search up to depth R
	for (int ii = 0; ii < (int)que.size(); ii++) {
		typename graph_traits<Graph>::vertex_descriptor cur_vD = que[ii];
		int cur_v = get(name_map, cur_vD);
		res.PB(cur_vD);
		if (dis[cur_v] == R) { continue; }
		//for (auto nei : graph[cur_v]) {
		typename graph_traits<Graph>::adjacency_iterator neiD, nei_end;
		for (tie(neiD, nei_end) = adjacent_vertices(cur_vD, graph); neiD != nei_end; ++neiD){
			int nei = get(name_map, *neiD);
			if (last_vis[nei] != phase_id && where_in_order[nei] > where_in_order[root] && (is_forb.empty() || !is_forb[nei])) {
				last_vis[nei] = phase_id;
				que.PB(*neiD);
				dis[nei] = dis[cur_v] + 1;
			}
		}
	}
	return res;
} // Function seems to have side effects!!

/**
 * Function calculates the vertices that contain the root vertex in their potential weakly R-reachable set for a given order
 * @param where_in_order, the position of the vertices in the order
 * @param phase_id, the current position to fill in the ordering
 * @return a vector containing all vertices that have root in their weakly R-reachable set
 */
template < typename Graph, typename VertexNameMap >
vector<typename graph_traits<Graph>::vertex_descriptor> ComputeSingleCluster(const Graph& graph,
                                 vector<int>& where_in_order,
                                 int R,
                                 vector<int>& is_forb,
                                 typename graph_traits<Graph>::vertex_descriptor u,
								 VertexNameMap name_map,
                                 int phase_id) {
	int root = get(name_map,u);
	//cout << "root: " <<root  << endl;
	vector<int> last_vis(num_vertices(graph)+1);
	vector<int> dis(num_vertices(graph)+1);
	vector<typename graph_traits<Graph>::vertex_descriptor> res;
	//if(root == 0) return res;
	if (!is_forb.empty() && is_forb[root]) { return {}; }
	last_vis[root] = phase_id;
	dis[root] = 0;
	vector<typename graph_traits<Graph>::vertex_descriptor> que{u};
	// breadth first search up to depth R
	for (int ii = 0; ii < (int)que.size(); ii++) {
		typename graph_traits<Graph>::vertex_descriptor cur_vD = que[ii];
		int cur_v = get(name_map, cur_vD);
		res.PB(cur_vD);
		if (dis[cur_v] == R) { continue; }
		//for (auto nei : graph[cur_v]) {
		typename graph_traits<Graph>::adjacency_iterator neiD, nei_end;
		for (tie(neiD, nei_end) = adjacent_vertices(cur_vD, graph); neiD != nei_end; ++neiD){
			int nei = get(name_map, *neiD);
			if (last_vis[nei] != phase_id && where_in_order[nei] > where_in_order[root] && (is_forb.empty() || !is_forb[nei])) {
				last_vis[nei] = phase_id;
				que.PB(*neiD);
				dis[nei] = dis[cur_v] + 1;
			}
		}
	}
	return res;
} // Function seems to have side effects!!

/**
 * Function calculates all weakly reachable sets for the given graph and order
 * @return a vector with the weakly reachable sets
 */
template < typename Graph, typename VertexNameMap, typename descVec >
vector<descVec> ComputeAllWReach(const Graph& graph,
									 VertexNameMap name_map,
                                     vector<int>& where_in_order,
                                     int R,
                                     vector<int> is_forb, descVec dVdummy) {
	int n = num_vertices(graph);
	vector<int> last_vis(n + 1, -1);
	vector<int> dis(n + 1);
	vector<descVec> res(n + 1);
	int ct = 1;
	typename graph_traits< Graph >::vertex_iterator root, end;
	for(tie(root, end) = vertices(graph); root != end; ++root){
		descVec cluster = ComputeSingleCluster(graph, where_in_order, R, is_forb, last_vis, dis, *root, name_map, ct);

		for (typename graph_traits<Graph>::vertex_descriptor vD : cluster) {
			int v = get(name_map, vD);
			res[v].PB(*root);
		}
		ct++;
	}
	return res;
}

/**
 * Function calculates the sizes of the weakly reachable sets for the given graph and order
 * @return a vector with the sizes of the weakly reachable sets
 */
template < typename Graph, typename VertexNameMap >
vector<int> ComputeWreachSzs(const Graph& graph, vector<int>& where_in_order, int R, VertexNameMap name_map) {
	int n = num_vertices(graph);
	vector<int> wreach_sz(n+1);
	vector<int> last_vis(n + 1, -1);
	vector<int> dis(n + 1);
	vector<int> is_forb;
	int ct = 1;
	typename graph_traits< Graph >::vertex_iterator root, end;
	for(tie(root, end) = vertices(graph); root != end; ++root){
		int rootName = get(name_map, *root);
		vector<typename graph_traits<Graph>::vertex_descriptor> cluster = ComputeSingleCluster(graph, where_in_order, R, is_forb, last_vis, dis, *root, name_map, ct);
		for (typename graph_traits<Graph>::vertex_descriptor v : cluster) {
		  // increase the wcol for all vertices in the cluster of root as they have
		  // root in their weakly r-reachable set
			wreach_sz[get(name_map, v)]++;
		}
		ct++;
	}
	return wreach_sz;
}

// Returns graph where u, v are connected iff dis(u, v) <= R
template < typename Graph, typename VertexNameMap, typename descVec >
Graph PowerGraph(const Graph& graph, int R, std::unordered_set<int>& forb, VertexNameMap name_map, descVec& vD_pg) {
	int n = num_vertices(graph) - 1;
	Graph pow_graph;
	typename property_map < Graph, vertex_name_t >::type name_map_pg = get(vertex_name, pow_graph);
	typename property_traits< typename property_map < Graph, vertex_name_t >::type >::value_type name;
	map<typename graph_traits<Graph>::vertex_descriptor,typename graph_traits<Graph>::vertex_descriptor> descMap;
	// for each vertex in graph add one to pow_graph and remember the mapping
	typename graph_traits< Graph >::vertex_iterator vi, vi_end;
	// add a dummy at pos. 0, since there is no shrinked id 0 and thus neither a corr. vertex descriptor
	vD_pg.push_back(-1);
	int ct = 1;
	for(tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi){

		typename graph_traits<Graph>::vertex_descriptor vd = add_vertex(pow_graph);
		descMap[*vi] = vd;
		name = ct;
		put(name_map_pg, vd, name);
		vD_pg.push_back(vd);
		ct++;
	}

	vector<int> last_vis(n + 1);
	vector<int> dis(n + 1);
	vector<int> where_in_order(n + 1);
	vector<int> is_forb;
	if (!forb.empty()) {
		is_forb.resize(n + 1);
		for (auto v : forb) {
		  is_forb[v] = 1;
		}
	}
	int ct2 = 1;
	typename graph_traits< Graph >::vertex_iterator rootD, end;
	for(tie(rootD, end) = vertices(graph); rootD != end; ++rootD){
		int root = get(name_map, *rootD);
		where_in_order[root] = -1; // hack, to make it think root is before everybody  in order

		vector<typename graph_traits<Graph>::vertex_descriptor> cluster = ComputeSingleCluster(graph, where_in_order, R, is_forb, last_vis, dis, *rootD, name_map, ct2);

		for (typename graph_traits<Graph>::vertex_descriptor v : cluster) {
			// get mapped vertices in pow_graph
			if(descMap[*rootD] != descMap[v] && !edge(descMap[v],descMap[*rootD],pow_graph).second)
				add_edge(descMap[*rootD], descMap[v], pow_graph);
		}
		where_in_order[root] = 0;
		ct2++;
	}

	return pow_graph;
}

/**
 * Function computes the degeneracy ordering for the passed graph
 * @param graph, the graph to calculate the degeneracy for
 * @param R, the radius
 * @return a pair containing the degeneracy and the corresponding ordering
 */
template < typename Graph, typename VertexNameMap, typename descVec>
pair<vector<int>, vector<int>> Degeneracy(const Graph& graph, int R, VertexNameMap name_map, descVec& removed_orderD) {
	int n = num_vertices(graph);
	std::unordered_set<int> dummy_forb;
	Graph pow_graph;
	descVec vD_pg;
	// pow_graph: only vertices of distance up to R are connected; format is an adjacency matrix
	pow_graph = PowerGraph(graph, R, dummy_forb, name_map, vD_pg);
	typename property_map < Graph, vertex_name_t >::type name_map_pg = get(vertex_name, pow_graph);

	int degeneracy = 0;
	// vector with n+1=graph.size() elements
	vector<int> degree(n + 1);
	// bucket[i] contains the vertex v if it has degree i
	vector<set<int>> buckets(n + 1);

	vector<bool> already_removed(n + 1);
	// nonempty_buckets contains the degrees recurring in the graph
	set<int> nonempty_buckets;
	// init data structures
	typename graph_traits< Graph >::vertex_iterator vi, vi_end;
	for(tie(vi, vi_end) = vertices(pow_graph); vi != vi_end; ++vi){
		int deg = out_degree(*vi, pow_graph);
		int v = get(name_map_pg, *vi);
		buckets[deg].insert(v);
		nonempty_buckets.insert(deg);
		degree[v] = deg;
	}

	vector<int> removed_order;
	vector<int> where_in_order(n + 1);
	int count = 1;
	while (!nonempty_buckets.empty()) {
		// get currently lowest degree
		int wh_bucket = *(nonempty_buckets.begin());

		degeneracy = max(degeneracy, wh_bucket);
		// get current vertex, that is the first one in list of lowest degree
		int v_to_remove = *(buckets[wh_bucket].begin());

		typename graph_traits<Graph>::vertex_descriptor v_to_removeD = vD_pg[v_to_remove];
		removed_order.PB(v_to_remove);
		where_in_order[v_to_remove] = n - count + 1;
		removed_orderD.PB(v_to_removeD);
		already_removed[v_to_remove] = true;
		buckets[wh_bucket].erase(v_to_remove);
		// if all vertices with degree wh_bucket have been looked at remove the degree wh_bucket from nonempty_buckets
		if (buckets[wh_bucket].empty()) {
			nonempty_buckets.erase(wh_bucket);
		}
		// update the current vertex's neighbor
		typename graph_traits<Graph>::adjacency_iterator neiD, nei_end;
		for (tie(neiD, nei_end) = adjacent_vertices(v_to_removeD, pow_graph); neiD != nei_end; ++neiD){
			int nei = get(name_map_pg, *neiD);
			if (already_removed[nei]) {
				continue;
			}
			buckets[degree[nei]].erase(nei);
			if (buckets[degree[nei]].empty()) {
				nonempty_buckets.erase(degree[nei]);
			}
			degree[nei]--;
			assert(degree[nei] >= 0);
			// in case bucket for new degree of nei was empty
			if (buckets[degree[nei]].empty()) {
				nonempty_buckets.insert(degree[nei]);
			}
			// since degree has decreased, the vertex has to be added to the next lower bucket
			buckets[degree[nei]].insert(nei);
		}
		count++;
	}
	cout << "degeneracy: " << degeneracy << endl;
	buckets.clear();
	degree.clear();
	already_removed.clear();
	nonempty_buckets.clear();
	reverse(removed_order.begin(), removed_order.end());
	reverse(removed_orderD.begin(), removed_orderD.end());

	return {removed_order, where_in_order};
}

template < typename Graph, typename VertexNameMap, typename descVec>
pair<vector<int>, vector<int>> ByWReachLeft(const Graph& graph, int R, VertexNameMap name_map, descVec& orderD, descVec vD){

	int n = num_vertices(graph);
	vector<int> where_in_order(n + 1, n + 1); // hacky hack to set where_in_order to n+1 for all not decided vertices
	vector<int> put(n + 1);
	vector<int> order;
	_wreach_szs.resize(n + 1);
	_deg.resize(n + 1);
	vector<int> last_vis(n + 1);
	vector<int> dis(n + 1);
	// doesn't contain anything, but expected as param from ComputeSingleCluster
	vector<int> is_forb_dummy;
	// contains all vertices, changes over time
	set<Vert> verts;

	typename graph_traits< Graph >::vertex_iterator vi, vi_end;
	for(tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi){
		int j = get(name_map, *vi);
		// save degree
	    _deg[j] = out_degree(*vi, graph);
	    // push vertex to data structure
	    verts.insert({j});
	}

	for (int i = 1; i <= n; i++) {
		// get next vertex to place in order: is first one of vertex set
		// which is ordered already by size of weakly reachable sets
	    int where_best = verts.begin()->id;
	    // erase it from set
	    verts.erase(verts.begin());
	    // remember where and that it was put in the order
	    where_in_order[where_best] = i;
	    put[where_best] = 1;
	    orderD.PB(vD[where_best]);
	    order.PB(where_best);
	    // calculate from which vertices it is weakly reachable (BFS)
	    descVec cluster = ComputeSingleCluster(graph, where_in_order, R, is_forb_dummy, last_vis, dis, vD[where_best], name_map, i);

	    for (typename graph_traits<Graph>::vertex_descriptor v: cluster) {
			// ignore the current vertex to be placed if found in the cluster
	    	int x = get(name_map, v);
			if (x == where_best) { continue; }
			// delete cluster vertex x in verts
			auto it = verts.find({x});
			verts.erase(it);
			// x' wreach_set now contains one more vertex
			_wreach_szs[x]++;
			// needs to be newly inserted to set
			// < was overloaded, so it is used for the insertion
			verts.insert({x});
		}
	}
	return {order, where_in_order};
}
/**
 * Function identifies the weak coloring number of the passed graph and the given order
 * and all the vertices with a size of their weakly reachable sets equal to it
 */
int ComputeWcol(vector<int> wreach_szs) {
	  int res = 0;
	  for (auto x : wreach_szs) {
	    res = max(x, res);
	  }
	  return res;
}


int main(int argc, char** argv) {

	/******************************* DEFINTIONS ********************************/

	typedef adjacency_list< listS, // Store out-edges of each vertex in a std::list
	vecS, // Store vertex set in a std::vector
	undirectedS, // The graph is directed
	property< vertex_name_t, int > // Add a vertex property
	> UndirectedGraph;
	typedef vector<graph_traits<UndirectedGraph>::vertex_descriptor> descVec;

	UndirectedGraph g; // use default constructor to create empty graph
	int n = 0;
	property_map < UndirectedGraph, vertex_name_t >::type name_map = get(vertex_name, g);
	typename property_traits< property_map < UndirectedGraph, vertex_name_t >::type >::value_type name;

	/******************************* GRAPH READ-IN ********************************/

	if (argc == 2 && string(argv[1]) == "--h") {
	  cerr<<"Usage: ./SimAnneal --in=graph.txtg --rad=radius --heur=heuristic [--o=output.txt]"<<endl;
	  cerr<<"rad - for rad=1 it computes typical degeneracy,"<<endl;
	  cerr<<"      for rad>1 it computes degeneracy on an auxiliary graph where"<<endl;
	  cerr<<"      two vertices are connected if their distance in original graph was <=R"<<endl;
	  cerr<<"      what serves as some kind of naive heuristic for wcol_R"<<endl;
	  cerr<<"o - if you want to print order in not default output file\n";
	  return 1;
	}

	GraphReader reader;
	/* Files */
	string graph_file, output_file, file_format, anneal_trace;
	/* directories */
	string graph_dir;
	string graph_name, rad_str = "1", heuristic;

	int R;
	try {
		FlagParser flag_parser;
		flag_parser.ParseFlags(argc, argv);
		graph_file = flag_parser.GetFlag("in", true);
		heuristic = flag_parser.GetFlag("heur", true);

		// set format
		if(hasEnding(graph_file, format_txtg)){
			file_format = format_txtg;
		}
		else if(hasEnding(graph_file, format_csv)){
			file_format = format_csv;
		}

		assert(graph_file.find(file_format) == graph_file.size() - file_format.size());

		int last_slash = -1;
		for (int i = 0; i < (int)graph_file.size(); i++) {
			if (graph_file[i] == '/') {
				last_slash = i;
			}
		}
		graph_dir = graph_file.substr(0, last_slash + 1);
		graph_name = graph_file.substr(last_slash + 1, (int)graph_file.size() - file_format.size() - last_slash - 1);

		rad_str = flag_parser.GetFlag("rad", true);
		try {
			R = stoi(rad_str);
		} catch (...) {
			cerr<<"Error: Radius must be a positive integer\n";
		}

		output_file = graph_dir + "orders/" + graph_name + ".simAnneal" + rad_str + ".txt";

		string cand_output_file = flag_parser.GetFlag("o", false);
		if (!cand_output_file.empty()) {
			output_file = cand_output_file;
		}
		flag_parser.Close();
	} catch (string err) {
		cerr<<"Error: "<<err<<endl;
		Err();
	}


	pair<int,vector<pair<string, string>>> res = reader.ReadGraphEdges(graph_file, file_format);
	n = res.st - 1;
	cout << "n: " << n << endl;
	descVec vD;
	// add a dummy at pos. 0, since there is no shrinked id 0 and thus neither a corr. vertex descriptor
	vD.push_back(-1);
	typename graph_traits < UndirectedGraph >::vertex_descriptor u;
	for(int v=1; v <= n; v++){
		// u = v-1, if the vertices of the input graph are denoted by integers
		u = add_vertex(g);
		name = v;
		put(name_map, u, name);
		vD.push_back(u);
		//cout << get(name_map, u) << endl;
	}

	// add edges
	for(auto e: res.nd){
		if(vD[reader.shrink_indices[e.st]] != vD[reader.shrink_indices[e.nd]])
			add_edge(vD[reader.shrink_indices[e.st]], vD[reader.shrink_indices[e.nd]], g);
	}
	cout << "edges added " << endl;

	/*pair<UndirectedGraph::edge_iterator, UndirectedGraph::edge_iterator> es = edges(g);
	copy(es.first, es.second, std::ostream_iterator<UndirectedGraph::edge_descriptor>{cout, "\n"});*/

	/******************************* SETUP OF INITIAL ORDER ********************************/

	int wcol = 0, wcol_test = 0;
	vector<int> wreach_szs(n+1,0);
	vector<int> order;
	descVec maxWcolVertD;
	vector<int> last_vis(n + 1);
	vector<int> dis(n + 1);
	// doesn't contain anything, but expected as param from ComputeSingleCluster
	vector<int> is_forb_dummy;
	std::unordered_set<int> forb;
	vector<int> where_in_order(n + 1);
	descVec orderD;
	vector<descVec>  wreachSets(n + 1);
	descVec dVdummy(1);

	if(heuristic != "none"){
		// orderD is passed by reference and modified in called functions
		if(heuristic == "deg"){
			tie(order, where_in_order) = Degeneracy(g, R, name_map, orderD);
			wreach_szs = ComputeWreachSzs(g, where_in_order, R, name_map);
		}
		else if(heuristic == "wReachLeft"){
			tie(order, where_in_order) = ByWReachLeft(g, R, name_map, orderD, vD);
			wreach_szs = _wreach_szs;

		}
		wcol = ComputeWcol(wreach_szs);
		wreachSets = ComputeAllWReach(g, name_map, where_in_order, R, is_forb_dummy, dVdummy);

	}
	else{
		// calculate the initial order, which simply is that of the vertex descriptor list
		// along the way we also calculate the size of the weakly R-reachable sets
		int count = 1;

		graph_traits< UndirectedGraph >::vertex_iterator i, end;
		for(tie(i, end) = vertices(g); i != end; ++i){
			int iName = get(name_map, *i);
			where_in_order[iName] = count;
			orderD.PB(*i);
			order.PB(iName);
			descVec cluster = ComputeSingleCluster(g, where_in_order, R, is_forb_dummy, last_vis, dis, *i, name_map, count);

			// for each vertex in the cluster add i to its weakly reachable set and increment the set size counter
			for(int v=0; v<cluster.size(); v++){
				int memberName = get(name_map, cluster[v]);
				wreachSets[memberName].PB(*i);
				wreach_szs[memberName]++;
				// compute wcol and remember each vertex v with wreachSets[v].size() == wcol
				if(wreach_szs[memberName] > wcol){
					maxWcolVertD.clear();
					maxWcolVertD.PB(cluster[v]);
					wcol = wreach_szs[memberName];
				}
				else if(wreach_szs[memberName] == wcol){
					maxWcolVertD.PB(cluster[v]);
				}
			}
			count++;

		}
	}

	cout << "initial order done" << endl;

	cout << "wcol at start: " << wcol << endl;
	cout << "check: wcol = "<< ComputeWcol(ComputeWreachSzs(g, where_in_order, R, name_map)) << endl;

	/******************************* SIMULATED ANNEALING ********************************/
	/* Try with reheating */

	/* initialize random seed: */
	srand (time(NULL));

	anneal_trace = graph_dir + "anneal_traces/" + graph_name + "_simAnneal" + rad_str + "_trace.txt";
	ofstream trace;
	InitOfstream(trace, anneal_trace);

	trace << "#t,il,wcol,swaps,rdVal,prob" << endl;
	trace << 1 << ",," << wcol << ",0,,0"<< endl;
	float swapLim_init = (float) n / (float) 10;
	int schedule_Start = 0.4*n;
	int schedule_End = 0.9*n;

	int ilLast = 0, tLast = 0;
	float slope = 0.002;
	int compCt = 0, df=0;
	int maxComputations = n * 400;

	while(compCt < maxComputations && schedule_Start < schedule_End){ // computation loop
		for(int t = schedule_Start; t < schedule_End; t++){ // schedule loop
			// temperature: how often to swap in relation to iterations
			int swapLim = max(1, (int) (swapLim_init * exp(-slope*t)));

			// trials to find a suitable neighbor
			for(int il=0; il< t; il++){ // trial loop
				// by mult. with (float) n/ (float) t, rdVal increases for lower temperatures
				float rdVal = min(0.99, (rand() % 100 + 35.0) * 0.01); //* (float) n/ (float) t;
				descVec orderD_copy = orderD;
				vector<descVec> wreachSets_copy = wreachSets;
				vector<int> where_in_order_copy = where_in_order, wreach_szs_copy = wreach_szs, order_copy = order;
				bool wcolInc = false;
				float leftAvg = 0.0;
				float rightAvg = 0.0;
				int left = 0, right = 0;
				for(int swapI=1; swapI < swapLim; swapI++){ // swap loop

					// pick a random neighbor, i.e. randomly swap two vertices in the given order
					// left and right denote the pos. of the vertices before the swap
					left = rand() % n;
					right = rand() % n;
					//cout << "left: " << left << ", right: " << right << endl;
					if(left == right){
						continue;
					}
					if(left > right){
						int temp = left;
						left = right;
						right = temp;
						//cout << "left: " << left << ", right: " << right << endl;
					}
					descVec clusterOldLeft = ComputeSingleCluster(g, where_in_order_copy, R, is_forb_dummy, orderD_copy[left], name_map, order_copy[left]);
					descVec clusterOldRight = ComputeSingleCluster(g, where_in_order_copy, R, is_forb_dummy, orderD_copy[right], name_map, order_copy[right]);
					// where_in_order has to be adapted alongside swapping
					swap(orderD_copy[left], orderD_copy[right]);
					swap(order_copy[left], order_copy[right]);
					where_in_order_copy[order_copy[right]] = right+1;
					where_in_order_copy[order_copy[left]] = left+1;
					// the clusters of vertices at positions left and right resp. after the swap
					descVec clusterNewLeft = ComputeSingleCluster(g, where_in_order_copy, R, is_forb_dummy, orderD_copy[right], name_map, order_copy[right]);
					descVec clusterNewRight = ComputeSingleCluster(g, where_in_order_copy, R, is_forb_dummy, orderD_copy[left], name_map, order_copy[left]);

					descVec diffLeft = getDifferenceOfVecs(clusterOldLeft, clusterNewLeft);
					// the vertex moved to the left loses the vertices that appear in its cluster now
					descVec diffRight = getDifferenceOfVecs(clusterNewRight, clusterOldRight);

					wreach_szs_copy[order_copy[right]] = wreach_szs_copy[order_copy[right]] + diffLeft.size();
					// update the weakly reachable sets for the vertices in diffLeft
					// they lose the vertex moved to the right
					for(typename graph_traits < UndirectedGraph >::vertex_descriptor neileft: diffLeft){
						int memberName = get(name_map, neileft);
						if(wreachSets_copy[memberName].size()  - 1 != getDifferenceOfVecs(wreachSets_copy[memberName], {orderD_copy[right]}).size()){
							printVec(g, wreachSets_copy[memberName], neileft, name_map, "wreach set ");
							printVec(g, getDifferenceOfVecs(wreachSets_copy[memberName], {orderD_copy[right]}), orderD_copy[right], name_map, "diff set after removing");
							printVec(g, clusterNewLeft, orderD_copy[right], name_map, "wrong cluster");
							cout << "wiom: " << where_in_order_copy[memberName] << ", " << right+1 << endl;
							df++;

						}
						wreachSets_copy[memberName] = getDifferenceOfVecs(wreachSets_copy[memberName], {orderD_copy[right]});
						wreach_szs_copy[memberName]--;
						// update the weakly reachable set for the vertex that was moved to the right
						wreachSets_copy[order_copy[right]].PB(neileft);

					}
					// update the weakly reachable sets for the vertices in diffRight
					// they gain the vertex moved to the left
					for(typename graph_traits < UndirectedGraph >::vertex_descriptor neiright: diffRight){
						int memberName = get(name_map, neiright);
						if(memberName != order_copy[right]){
							if(wreachSets_copy[memberName].size()  + 1 != getUnionOfVecs(wreachSets_copy[memberName], {orderD_copy[left]}).size()){
								cout << "t, il, swap: " << t << ", " << il << ", " << swapI << endl;
								printVec(g, wreachSets_copy[memberName], neiright, name_map, "wreach set ");
								printVec(g, getUnionOfVecs(wreachSets_copy[memberName], {orderD_copy[left]}), orderD_copy[left], name_map, "union set after adding");
								printVec(g, clusterNewRight, orderD_copy[left], name_map, "wrong cluster");
								cout << "wiom: " << where_in_order_copy[memberName] << ", " << left+1 << endl;
								df++;

							}
							wreachSets_copy[memberName] = getUnionOfVecs(wreachSets_copy[memberName], {orderD_copy[left]});
							wreach_szs_copy[memberName]++;
							wreach_szs_copy[order_copy[left]]--;
						}
					}
					// update the weakly reachable set for the vertex that was moved to the left
					wreachSets_copy[order_copy[left]] = getDifferenceOfVecs(wreachSets_copy[order_copy[left]], diffRight);
				}

				int wcolNew = ComputeWcol(wreach_szs_copy);

				wcolInc = wcolNew >= wcol;
				float prob = getProb(t,wcol, wcolNew, n);

				if(!wcolInc || prob >= rdVal){
					// found better solution or accepting worse one
					// update wcol and data structures
					wcol = wcolNew;
					where_in_order = where_in_order_copy, wreach_szs = wreach_szs_copy, wreachSets = wreachSets_copy;
					orderD = orderD_copy, order = order_copy;
					if(wcolInc && prob >= rdVal){
						trace << t << "," << il << "," << wcol << "," << swapLim << "," << rdVal << "," << prob << endl;
					}
					else{
						trace << t << "," << il << "," << wcol << "," << swapLim << endl;
					}
					break; // trial loop
				}
				ilLast = il;
			} // end trial loop
			// since we've been stuck in a local minimum we need to "heat up" again
			// we leave the loop at the current low temperature and do another round of annealing
			if(ilLast + 1 >= t){
				cout << "heat up -- compCt: " << compCt << ", wcol: " << wcol << endl;
				trace << "###heat up" << endl;
				schedule_Start += (t - schedule_Start)*0.15;
				break; // schedule loop
			}
			tLast = t;
		} // end schedule loop
		compCt += tLast;
	} // end computation loop
	trace.close();

	cout << "wcol final: "<< wcol << endl;
	vector<int> ws = ComputeWreachSzs(g, where_in_order, R, name_map);
	//assert(wcol == ComputeWcol(ComputeAllWReach(g, name_map, where_in_order, R, is_forb_dummy, dVdummy)));
	cout << "check: wcol = "<< ComputeWcol(ws) << endl;
	for(int j=0; j < order.size(); j++){
		assert(order[j] = orderD[j] + 1);
	}
	ofstream out;
	InitOfstream(out, output_file);
	cout << "Order final" <<endl;
	for (int v : order) {
		//cout<< v <<endl;
		out << reader.GetOriginalFromMapped(v) << " ";
	}
	out << endl;
	out.close();
}

