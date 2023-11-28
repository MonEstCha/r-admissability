/*
 * Admissability_BoostVersion.cpp
 *
 *  Created on: 02.11.2023
 *      Author: monique
 */

#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include "Headers.hpp"
#include "FilesOps.hpp"
#include "FlagParser.hpp"
#include "ReadTxt.hpp"

using namespace boost;


const string format_dat = ".dat";

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
  cerr<<"Usage: ./Degeneracy --in=graph.txtg --rad=radius [--o=output.txt]"<<endl;
  cerr<<"--h for help\n";
  exit(1);
}

template < typename Graph, typename VertexNameMap >
int execDFS(typename graph_traits< Graph >::vertex_descriptor u, const Graph& g,
	    VertexNameMap name_map, int putToOrder[], int visited[], int R){
	int rootN = get(name_map, u);
	visited[rootN] = 1;
	int admissVert = 0;
	R--;
	// does root have a neighbor outside of the order?
	typename graph_traits < Graph >::adjacency_iterator vi, vi_end;
	for (boost::tie(vi, vi_end) = adjacent_vertices(u, g); vi != vi_end; ++vi){
		int child = get(name_map, *vi);
		// heuristical behavior: just take the first non-visited one
		if(visited[child] == 0){
			if(putToOrder[child] == 0){
				visited[child] = 1;
				return 1;
			}
			// else traverse further on vertices in order
			if(R > 0){
				admissVert = execDFS(*vi, g, name_map, putToOrder, visited, R);
				if (admissVert == 1){
					return 1;
				}
			}
		}
	}

	visited[rootN] = 0;
	return 0;
}

template < typename Graph, typename VertexNameMap >
int approxAdmOnWholeGraph(typename graph_traits< Graph >::vertex_descriptor u, const Graph& g,
	    VertexNameMap name_map, int putToOrder[], int n, int R){
	int approx_rAdmiss = 0;
	int visited[n+1];
	for(int i=0; i<n+1; i++){
		visited[i] = 0;
	}
	// looking at radius 1 means just looking at neighbors
	R--;
	visited[get(name_map, u)] = 1;
	typename graph_traits < Graph >::adjacency_iterator vi, vi_end;
	for (boost::tie(vi, vi_end) = adjacent_vertices(u, g); vi != vi_end; ++vi){
		int i = get(name_map, *vi);
		if(putToOrder[i] == 0){
			// get 1-reachable vertices, i.e. neighbors that are not in order yet and therefore greater than u w.r.t. order
			approx_rAdmiss++;
		}
		else if(R > 0){
			approx_rAdmiss += execDFS(*vi, g, name_map, putToOrder, visited, R);
		}
		// make sure vertex is not visited again in case it was released
		visited[i] = 1;
	}

	return approx_rAdmiss;
}

template < typename Graph, typename VertexNameMap >
void output_adjacent_vertices(std::ostream& out,
    typename graph_traits< Graph >::vertex_descriptor u, const Graph& g,
    VertexNameMap name_map)
{
    typename graph_traits< Graph >::adjacency_iterator vi, vi_end;
    out << get(name_map, u) << " -> { ";
    for (boost::tie(vi, vi_end) = adjacent_vertices(u, g); vi != vi_end; ++vi)
        out << get(name_map, *vi) << " ";
    out << "}" << std::endl;
}

template < typename Graph, typename VertexNameMap >
int count_adj_vertices(typename graph_traits< Graph >::vertex_descriptor u, const Graph& g,
    VertexNameMap name_map, int putToOrder[])
{
	int ct = 0;
    typename graph_traits< Graph >::adjacency_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = adjacent_vertices(u, g); vi != vi_end; ++vi){
    	if(putToOrder[get(name_map, *vi)] == 0)
    		ct++;
    }
    return ct;
}

template < typename NameMap > class name_equals_t
{
public:
    name_equals_t(const std::string& n, NameMap map)
    : m_name(n), m_name_map(map)
    {
    }
    template < typename Vertex > bool operator()(Vertex u) const
    {
        return get(m_name_map, u) == m_name;
    }

private:
    std::string m_name;
    NameMap m_name_map;
};

// object generator function
template < typename NameMap >
inline name_equals_t< NameMap > name_equals(
    const std::string& str, NameMap name)
{
    return name_equals_t< NameMap >(str, name);
}

int main(int argc, char** argv)
{

	typedef adjacency_list< listS, // Store out-edges of each vertex in a std::list
	vecS, // Store vertex set in a std::vector
	undirectedS, // The graph is directed
	property< vertex_name_t, int > // Add a vertex property
	> UndirectedGraph;

	UndirectedGraph g; // use default constructor to create empty graph
	int n = 0;
	property_map < UndirectedGraph, vertex_name_t >::type name_map = get(vertex_name, g);

	if (argc == 2 && string(argv[1]) == "--h") {
	  cerr<<"Usage: ./Degeneracy --in=graph.txtg --rad=radius [--o=output.txt]"<<endl;
	  cerr<<"rad - for rad=1 it computes typical degeneracy,"<<endl;
	  cerr<<"      for rad>1 it computes degeneracy on an auxiliary graph where"<<endl;
	  cerr<<"      two vertices are connected if their distance in original graph was <=R"<<endl;
	  cerr<<"      what serves as some kind of naive heuristic for wcol_R"<<endl;
	  cerr<<"o - if you want to print order in not default output file\n";
	  return 1;
	}

	GraphReader reader;
	string graph_file, output_file, file_format;
	int R;
	try {
		FlagParser flag_parser;
		flag_parser.ParseFlags(argc, argv);
		graph_file = flag_parser.GetFlag("in", true);
		  //debug(graph_file);

		// set format
		if(hasEnding(graph_file, format_txtg)){
			file_format = format_txtg;
		}
		else if(hasEnding(graph_file, format_csv)){
			file_format = format_csv;
		}
		else if(hasEnding(graph_file, format_dat)){
			file_format = format_dat;
		}

		assert(graph_file.find(file_format) == graph_file.size() - file_format.size());

		int last_slash = -1;
		for (int i = 0; i < (int)graph_file.size(); i++) {
			if (graph_file[i] == '/') {
				last_slash = i;
			}
		}
		string graph_dir = graph_file.substr(0, last_slash + 1);
		string graph_name = graph_file.substr(last_slash + 1, (int)graph_file.size() - file_format.size() - last_slash - 1);

		string rad_str = flag_parser.GetFlag("rad", true);
		try {
			R = stoi(rad_str);
		} catch (...) {
			cerr<<"Error: Radius must be a positive integer\n";
		}

		output_file = graph_dir + "orders/" + graph_name + ".admBoost" + rad_str + ".txt";
		//debug(output_file);
		string cand_output_file = flag_parser.GetFlag("o", false);
		if (!cand_output_file.empty()) {
			output_file = cand_output_file;
			//debug(output_file);
		}
		flag_parser.Close();
	} catch (string err) {
		cerr<<"Error: "<<err<<endl;
		Err();
	}

	typename property_traits< property_map < UndirectedGraph, vertex_name_t >::type >::value_type name;

	pair<int,vector<pair<string, string>>> res = reader.ReadGraphEdges(graph_file, file_format);
	n = res.st - 1;
	vector<typename graph_traits < UndirectedGraph >::vertex_descriptor> vD;
	vD.push_back(-1);
	typename graph_traits < UndirectedGraph >::vertex_descriptor u;
	for(int v=1; v <= n; v++){
		u = add_vertex(g);
		name = v;
		put(name_map, u, name);
		vD.push_back(u);
	}
	// add edges

	for(auto e: res.second){
		add_edge(vD[reader.shrink_indices[e.first]], vD[reader.shrink_indices[e.second]], g);
	}
	cout << "edges added " << endl;

	int putToOrder[n+1];
	for(int v=0; v<n+1; v++){
		putToOrder[v] = 0;
	}

	vector<int> order;
	int minAdmiss = n-1, maxMinAdmiss = 0, minAdmissV = n-1, adm=0;
	for(int k=0; k<n; k++){ // order step
		minAdmiss = n-1;
		graph_traits< UndirectedGraph >::vertex_iterator i, end;
		for(tie(i, end) = vertices(g); i != end; ++i){
			if(get(name_map, *i) == 0){
				continue;
			}
			int v = get(name_map, *i);
			if(putToOrder[v] == 0){
				adm = approxAdmOnWholeGraph(*i, g, get(vertex_name, g), putToOrder, n, R);
				if(minAdmiss > adm){
					minAdmiss = adm;
					minAdmissV = v;
				}
			}
		}

		if(minAdmiss > maxMinAdmiss){
			maxMinAdmiss = minAdmiss;
		}
		order.push_back(minAdmissV);
		putToOrder[minAdmissV] = 1;

	}

	cout << "maxMinAdmiss: " << maxMinAdmiss <<endl;

	ofstream out;
	InitOfstream(out, output_file);
	reverse(order.begin(), order.end());
	cout << "Order final" <<endl;
	for (int v : order) {
		out << reader.GetOriginalFromMapped(v) << " ";
	}
	out << endl;
	out.close();
	return 0;
}









