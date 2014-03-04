kpath_centrality: betweenness.o fibheap.o kpath.o readgml.o main_kpath.o
	g++ betweenness.o fibheap.o kpath.o readgml.o main_kpath.o -o kpath_centrality
fibheap.o: fibheap.h fibheap.cpp
kpath.o: kpath.h kpath.cpp
readgml.o: readgml.h readgml.cpp
betweenness.o: betweenness.h betweenness.cpp
main_kpath.o: readgml.h kpath.h betweenness.h main_kpath.cpp
