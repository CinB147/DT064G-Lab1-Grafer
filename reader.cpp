#include "reader.h"
#include <iostream>
#include <fstream>
#include <vector>

enum token {
    COMMENT, META, EDGE, END_OF_FILE
};

token get_line_type(std::istream& is) {
    switch (is.peek()) {
    case std::istream::traits_type::eof():      return END_OF_FILE;
    case '#':                                   return COMMENT;
    case 'M':                                   return META;
    };
    return EDGE;
}

edge read_edge(std::istream& is) {
    edge e;
    is >> e.n1 >> e.n2 >> e.weight;
    std::getline(is, e.description);
    return e;
}

void read_meta(std::istream& is, meta_t& meta) {
    char discard;
    node_id_t vertex_id;
    std::string name;
    is >> discard >> vertex_id; // is = discard = vertex_id which means is = vertex_id
    std::getline(is, name);
    meta[vertex_id] = name; // stores name in meta map using vertex_id as the key
}

adjacency_list_t parse_file(std::string filename) {
    std::ifstream in(filename);
    token l;
    meta_t meta;
    edge_list_t edge_list;
    while ((l = get_line_type(in)) != END_OF_FILE) {
        edge e;
        switch (l) {
        case token::EDGE:
            e = read_edge(in);
            edge_list.push_back(e);
            break;
        case token::META:
            read_meta(in, meta);
            break;
        default:
            std::string comment;
            std::getline(in, comment);
        }
    }
    return std::make_pair(meta, edge_list);

}

void display_adjacency_list(adjacency_list_t adj_list) {
    std::cout << "Meta:" << std::endl;
    for (auto& meta_pair : adj_list.first) {
        std::cout << meta_pair.first << " " << meta_pair.second << std::endl;
    }
    std::cout << std::endl << "Edges:" << std::endl;
    for (auto const& edge : adj_list.second) {
        std::cout << edge.n1 << " " << edge.n2 << " " << edge.weight << " " << edge.description << std::endl;
    }
}


