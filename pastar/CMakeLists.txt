

add_subdirectory( pastar_functions )

set (PASTAR_SRC backtrace.cpp Coord.cpp CoordHash.cpp Cost.cpp HeuristicHPair.cpp lz4sup.cpp msa_options msa_pastar_main.cpp Node.cpp PairAlign.cpp PAStar.cpp read_fasta.cpp Sequences.cpp TimeCounter.cpp WeightedSP.cpp)

add_executable(pastar ${PASTAR_SRC})
target_link_libraries(pastar pastar_functions)