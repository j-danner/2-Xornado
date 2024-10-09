add_library(
    graph
    #misc
    misc.hpp
    #solve
    solve.hpp
    solve.cpp
    #xlit
    LA/lineral.hpp
    LA/lineral.cpp
    #xsys
    LA/lineqs.hpp
    LA/lineqs.cpp 
    #graph
    graph/graph.hpp
    graph/graph_al.hpp
    graph/graph_al.cpp
    graph/graph_lhgr.hpp
    graph/graph_lhgr.cpp
    #impl_graph
    impl_graph.hpp
    impl_graph.cpp
    #vert_label
    vl/vl.hpp
    vl/vl_trie.hpp
    vl/vl_trie.cpp
    vl/vl_hmap.hpp
    vl/vl_hmap.cpp
    #hashmaps
    parallel-hashmap/phmap.h
    robin_hood-3.11.5/robin_hood.h
)

#find openmp package
find_package(OpenMP REQUIRED)

target_link_libraries(graph PUBLIC ${EXTRA_LIBS} PRIVATE OpenMP::OpenMP_CXX m4ri)

