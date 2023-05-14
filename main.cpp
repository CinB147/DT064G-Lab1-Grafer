#include <iostream>
#include <locale>
#include <stack>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <limits>
#include "reader.h"




std::unordered_map<int, bool> visited;
std::unordered_map<int, int> parent;
std::unordered_map<int, double> dist;
std::stack<int> s;
typedef int node_id_t;
typedef double weight_t;
typedef std::pair<node_id_t, weight_t> neighbor_t;
typedef std::vector<std::vector<neighbor_t>> adjacency_list;

const int INF = INT_MAX;

bool is_valid_graph(const adjacency_list_t& adj_list);
void dfs(const adjacency_list_t& adj_list, node_id_t node, std::unordered_set<node_id_t>& visited_nodes);

bool is_validbfs(const std::vector<node_id_t>& nodes, const std::vector<edge>& edges);

void validbfs(node_id_t node, std::vector<bool>& visited, std::vector<std::vector<node_id_t>>& adj_list);
void dfs(int start, int end, const adjacency_list_t& adj_list);
void print_shortest_path(int start, int end);
void dijkstra(const adjacency_list_t& adj_list, node_id_t start_node, node_id_t end_node);




int main() {

    std::locale::global(std::locale("en_US.utf8")); //supports the character encoding sv-lang
    adjacency_list_t adj_list = parse_file("adjlist.txt");

    display_adjacency_list(adj_list);

    //Extract edge list
    edge_list_t edges = adj_list.second;

    //Check if it's a valid graph
    std::vector<node_id_t> nodes;
    for (const auto& meta_pair : adj_list.first) {
        nodes.push_back(meta_pair.first);
    }

    std::cout << "-----------------dfs------------------ \n";


    if (is_valid_graph(adj_list)) {
        std::cout << "The graph is valid" << std::endl;
    }
    else {
        std::cout << "The graph is invalid" << std::endl;

    }
    std::cout << "-----------------bfs------------------ \n";


    // Check if the graph is valid
    bool is_graph2 = is_validbfs(nodes, edges);

    // Print the result
    if (is_graph2) {
        std::cout << "The graph is valid." << std::endl;
    }
    else {
        std::cout << "The graph is not valid." << std::endl;
    }

    std::cout << "------------------search----------------- \n";

    int start;
    int end;

    std::cout << "Enter start node : ";
    std::cin >> start;
    std::cout << "Enter end node : ";
    std::cin >> end;
    std::cout << "------------------dfs-search----------------- \n";


    dfs(start, end, adj_list);

    if (visited[end]) {
        print_shortest_path(start, end);
    }
    else {
        std::cout << "There is no path from node " << start << " to node " << end << "\n";
    }


    std::cout << "------------------dijkstra-search----------------- \n";

    dijkstra(adj_list, start, end);

    return 0;
}


bool is_valid_graph(const adjacency_list_t& adj_list) {
    // Create an unordered set to keep track of visited nodes, unordered_set is faster than using a vector for keeping track 
    // of visited nodes is that unordered_set is implemented using a hash table(binary heap), which provides constant-time complexity for the most common operations
    // hash table is a container that can store a pair value where, it can take 'data' key as input and return int value of data key
    std::unordered_set<node_id_t> visited_nodes;

    // Start DFS from the first node in the adjacency list
    dfs(adj_list, adj_list.first.begin()->first, visited_nodes);

    // if numbers of visited node is the same as adj_list size ,returns true
    return visited_nodes.size() == adj_list.first.size();
}

void dfs(const adjacency_list_t& adj_list, node_id_t node, std::unordered_set<node_id_t>& visited_nodes) {
    // set current node as visited
    visited_nodes.insert(node);

    // Iterate over all neighbors of the current node, loop over all edges in the edge_list_t type of adj_list. 
    for (auto& edge : adj_list.second) {
        if (edge.n1 == node && visited_nodes.count(edge.n2) == 0) {  // 'edge.n1 == node' ensures that the edge is edge is connected to one of the nodes. 
                                                                     // ex edge is connects to nodes A and B, which is a link between A and B, then this edge is incident to both nodes A and B.
            // Recursively call DFS on unvisited neighbors
            dfs(adj_list, edge.n2, visited_nodes);
        }

        //The count() function of the std::unordered_set container returns the number of elements with a given key in the container.
        //In this case, the key is the node ID represented by edge.n1.
        //if  visited_nodes.count(edge.n1) == 0 it's mean that the node has not been yet visisted and it will returns dfs on that node 
        //otherwise it will move to the next edge check its adjacent nodes. it's means that if there is a connection between edge and nodes.
        else if (edge.n2 == node && visited_nodes.count(edge.n1) == 0) {
            dfs(adj_list, edge.n1, visited_nodes);
        }
    }
}



//--------BFS---------------//
bool is_validbfs(const std::vector<node_id_t>& nodes, const std::vector<edge>& edges) {
    // Create the adjacency list
    std::vector<std::vector<node_id_t>> adj_list(nodes.size());
    for (const auto& edge : edges) {
        adj_list[edge.n1].push_back(edge.n2);
        adj_list[edge.n2].push_back(edge.n1);
    }

    // Create a vector to keep track of visited nodes
    std::vector<bool> visited(nodes.size(), false);

    // Start BFS from the first node
    validbfs(nodes[0], visited, adj_list);

    // Check if all nodes are reachable from the first node
    for (bool is_visited : visited) {
        if (!is_visited) {
            return false;
        }
    }

    return true;
}

void validbfs(node_id_t node, std::vector<bool>& visited, std::vector<std::vector<node_id_t>>& adj_list) {
    std::queue<node_id_t> q;
    q.push(node);
    visited[node] = true;

    while (!q.empty()) {
        //check the point at the first node in queue
        node_id_t curr_node = q.front();
        //removes current node that been visited on queue
        q.pop();

        //the function checks if the neighbor has not been visited before by checking the corresponding value in the visited vector.
        //If the neighbor has not been visited, the function marks it as visited and adds it to the queue.

        for (node_id_t neighbour : adj_list[curr_node]) {

            //adding unvisited neighbors to the queue 
            if (!visited[neighbour]) {
                visited[neighbour] = true;
                q.push(neighbour);
            }
        }
    }
}




//----------------//


void print_shortest_path(int start, int end) {
    std::cout << "Shortest path from node " << start << " to node " << end << " is: ";
    std::vector<int> path;
    int current = end;
    while (current != start) {
        path.push_back(current);
        current = parent[current];
    }
    path.push_back(start);
    for (auto iter = path.rbegin(); iter != path.rend(); ++iter) {
        std::cout << *iter;
        if (*iter != end) {
            std::cout << " -> ";
        }
    }
    std::cout << "\n";
    std::cout << "Total distance is: " << dist[end] << "\n";
}

void dfs(int start, int end, const adjacency_list_t& adj_list) {
    visited[start] = true;
    s.push(start);
    while (!s.empty()) {
        int current = s.top();
        s.pop();
        for (auto& edge : adj_list.second) {
            if (edge.n1 == current && !visited[edge.n2]) {
                visited[edge.n2] = true;
                dist[edge.n2] = dist[current] + edge.weight;
                parent[edge.n2] = current;
                s.push(edge.n2);
                if (edge.n2 == end) {
                    return;
                }
            }
        }
    }
}

//-----------------------//

void dijkstra(const adjacency_list_t& adj_list, node_id_t start_node, node_id_t end_node) {
    std::priority_queue<std::pair<double, node_id_t>, std::vector<std::pair<double, node_id_t>>, std::greater<>> pq;
    std::unordered_map<node_id_t, double> dist;
    std::unordered_map<node_id_t, node_id_t> parent;

    // Initialize distances to infinity and parent pointers to null
    for (auto const& meta_pair : adj_list.first) {
        dist[meta_pair.first] = INF;
        parent[meta_pair.first] = -1;
    }

    // Set distance to starting node to 0 and add it to the priority queue
    dist[start_node] = 0;
    pq.push({ 0, start_node });

    while (!pq.empty()) {
        node_id_t curr_node = pq.top().second;
        pq.pop();

        // If we've reached the end node, stop searching
        if (curr_node == end_node) {
            break;
        }
        node_id_t neighbor;
        bool updated = false;
        for (const auto& edge : adj_list.second) {
            if (edge.n1 == curr_node || edge.n2 == curr_node) { //will check at both way ex 19 to 3 and 3 to 19
                if (edge.n1 == curr_node) { neighbor = edge.n2; }
                else { neighbor = edge.n1; }
      
                double weight = edge.weight;
                if (dist[curr_node] + weight < dist[neighbor]) {
                    dist[neighbor] = dist[curr_node] + weight;
                    parent[neighbor] = curr_node;
                    updated = true;
                }
            }
        }

        if (updated) {
            // The shortest path to at least one neighbor node was updated,
            // so we need to re-calculate the priorities of all nodes in the priority queue.
            for (const auto& edge : adj_list.second) {
                node_id_t neighbor = edge.n2;
                if (dist[neighbor] < INF) {
                    pq.push({ dist[neighbor], neighbor });
                }
            }
        }

    }

    // If we didn't find a path to the end node, it's not reachable
    if (dist[end_node] == INF) {
        std::cout << "No path found from node " << start_node << " to node " << end_node << std::endl;
        return;
    }

    // Construct the path from end node to start node using parent pointers
    std::vector<node_id_t> path;
    node_id_t curr = end_node;
    while (curr != -1) {
        path.push_back(curr);
        curr = parent[curr];
    }
    std::reverse(path.begin(), path.end());

    // Display the shortest path and distance
    std::cout << "Shortest path from node " << start_node << " to node " << end_node << ": ";
    for (auto it = path.begin(); it != path.end(); ++it) {
        std::cout << *it;
        if (it != path.end() - 1) { //prevent "->" ends of node
            std::cout << " -> ";
        }
    }
    std::cout << std::endl;
    std::cout << "Shortest distance: " << dist[end_node] << std::endl;
}
