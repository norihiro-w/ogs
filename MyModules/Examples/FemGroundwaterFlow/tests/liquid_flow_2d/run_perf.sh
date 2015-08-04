perf record --call-graph dwarf  ~/ogs/ogs6/ogs6-thmc/BuildClangRelease/bin/ogs6 q_quad
perf report --call-graph
