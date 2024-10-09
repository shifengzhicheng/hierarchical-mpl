# 加载扩展库
load ./build/lib/libgraphextension.so Graphextension


# 读取图文件
set graph_file "graph.txt"
set read_result [read_graph $graph_file]
puts $read_result

# 执行DFS遍历，起始节点为1，日志文件为"dfs.log"
set dfs_result [dfs 1 "dfs.log"]
puts $dfs_result
