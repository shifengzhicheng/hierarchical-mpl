#include <tcl.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <sstream>
using namespace std;


typedef map<int, vector<int> > Graph;

static Graph graph;  // 存储读取的图
static set<int> visited;  // 记录访问过的节点
static ofstream logFile;  // 日志文件
int ReadGraphCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]) {
    if (argc != 2) {
        Tcl_SetResult(interp, const_cast<char *>("Usage: read_graph filename"), TCL_STATIC);
        return TCL_ERROR;
    }

    const char *filename = argv[1];
    ifstream infile(filename);
    if (!infile) {
        string err = "Cannot open file: ";
        err += filename;
        Tcl_SetResult(interp, const_cast<char *>(err.c_str()), TCL_VOLATILE);
        return TCL_ERROR;
    }

    graph.clear();
    int nodeCount;
    infile >> nodeCount;

    string line;
    getline(infile, line);  // 读取剩余的换行符

    while (getline(infile, line)) {
        if (line.empty()) continue;
        istringstream iss(line);
        int node;
        iss >> node;
        int neighbor;
        while (iss >> neighbor) {
            graph[node].push_back(neighbor);
        }
    }

    infile.close();
    Tcl_SetResult(interp, const_cast<char *>("Graph loaded successfully"), TCL_STATIC);
    return TCL_OK;
}

void dfs(int node) {
    visited.insert(node);
    logFile << "Visited node: " << node << endl;
    for (int neighbor : graph[node]) {
        if (visited.find(neighbor) == visited.end()) {
            dfs(neighbor);
        }
    }
}

int DFSCmd(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[]) {
    if (argc != 3) {
        Tcl_SetResult(interp, const_cast<char *>("Usage: dfs start_node log_filename"), TCL_STATIC);
        return TCL_ERROR;
    }

    int startNode = atoi(argv[1]);
    const char *logFilename = argv[2];

    if (graph.find(startNode) == graph.end()) {
        Tcl_SetResult(interp, const_cast<char *>("Start node not found in graph"), TCL_STATIC);
        return TCL_ERROR;
    }

    visited.clear();
    logFile.open(logFilename);
    if (!logFile) {
        string err = "Cannot open log file: ";
        err += logFilename;
        Tcl_SetResult(interp, const_cast<char *>(err.c_str()), TCL_VOLATILE);
        return TCL_ERROR;
    }

    dfs(startNode);
    logFile.close();

    Tcl_SetResult(interp, const_cast<char *>("DFS traversal completed"), TCL_STATIC);
    return TCL_OK;
}

extern "C" int Graphextension_Init(Tcl_Interp *interp) {
    if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL) {
        return TCL_ERROR;
    }

    // 注册命令
    Tcl_CreateCommand(interp, "read_graph", ReadGraphCmd, NULL, NULL);
    Tcl_CreateCommand(interp, "dfs", DFSCmd, NULL, NULL);

    // 提供包信息
    Tcl_PkgProvide(interp, "GraphExtension", "1.0");

    return TCL_OK;
}
