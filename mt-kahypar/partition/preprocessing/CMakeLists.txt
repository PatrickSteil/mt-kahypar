set(PreprocessingSources
        community_detection/parallel_louvain.cpp
        community_detection/local_moving_modularity.cpp
        filtering/filter_graph.cpp
        filtering/dimacs_io.cpp
        filtering/graph_contraction.cpp
        filtering/parallel_connectivity.cpp
        filtering/cut_signatures.cpp
        filtering/tiny_cut_detection.cpp
        filtering/dinic_max_flow.cpp
        filtering/natural_cut_detection.cpp
        filtering/fragment_assembly.cpp)

target_sources(MtKaHyPar-Sources INTERFACE ${PreprocessingSources})
target_sources(MtKaHyPar-ToolsSources INTERFACE ${PreprocessingSources})
