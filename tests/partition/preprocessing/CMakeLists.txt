target_sources(mtkahypar_tests PRIVATE
        louvain_test.cc
        filtering/filter_graph_test.cc
        filtering/dimacs_io_test.cc
        filtering/union_find_test.cc
        filtering/graph_contraction_test.cc
        filtering/parallel_connectivity_test.cc
        filtering/cut_signatures_test.cc
        filtering/tiny_cut_detection_test.cc
        filtering/dinic_max_flow_test.cc
        filtering/natural_cut_detection_test.cc
        )
