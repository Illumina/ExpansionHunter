# Path Families

## Files

* [include/graphcore/PathFamily.hh](/include/graphcore/PathFamily.hh)
* [src/graphcore/PathFamily.hh](/src/graphcore/PathFamily.cpp)
* [tests/PathFamilyTest.cpp](/tests/PathFamilyTest.cpp)

## Definitions

If enumerating all [paths](paths.md) in a set becomes impractical, 
we can define a set of paths implicitly by a set of edges.
We call sets of paths defined in this way *path families*.
We construct a path family from
a set of edges *F = \{ f<sub>1</sub>, f<sub>2</sub>, ... f<sub>l</sub> \}*. 
A path *P* is in the path family *F* iff. at least one node in *P* 
is connected via an edge in *F* and for any node *p<sub>n</sub>* in *P*
both of these conditions hold:

* No incoming edge into *p<sub>n</sub>* exists in *F* or edge *(p<sub>n-1</sub>, p<sub>n</sub>)* is in *F*.
* No outgoing edge from *p<sub>n</sub>* exists in *F* or edge *(p<sub>n</sub>, p<sub>n+1</sub>)* is in *F*.

Informally, a path family F consists of paths that enter and exit nodes 
only through edges in F when such edges exist. 
A path does not have to cover the entire path family 
(either all nodes or all edges), but may start in the middle or 
enter/exit through any edge into a node that has no incoming/outgoing
edges in *F*.

In graph-tools, path families can be created as follows:

```c++
Graph G = makeDeletionGraph("TTT", "AT", "CCCCC");
PathFamily family(&G);
// family of all possible alleles
family.addEdge(0, 1);
family.addEdge(1, 2);
family.addEdge(0, 2);
```

The primary operation provided by path families is a query
operation which tells us whether a given path is contained in 
a path family:

```c++
Path path(&graph, 0, { 0, 1, 2 }, 0);
assert( family.containsPath(path) );
```

We can also test if a path family includes another family:

```c++
PathFamily family2(&G);
family2.addEdge(0, 1);
family2.addEdge(1, 2);

assert( family.includes(family2) );
assert( !family2.includes(family) );
```
# Path Family Operations

## Files

* [include/graphcore/PathFamilyOperations.hh](/include/graphcore/PathFamilyOperations.hh)
* [src/graphcore/PathFamilyOperations.cpp](/src/graphcore/PathFamilyOperations.cpp)
* [tests/PathFamilyOperationsTest.cpp](/tests/PathFamilyOperationsTest.cpp)

## Operations

### Generate path segments in a family
```c++
list<Path> getPathSegmentsForFamily(PathFamily const& family);
```

This function generates a set of maximum-length path segments for a family,
where each path that is the result of merging two path segments
is also contained in the path family.

### Enumerate path combinations in a family

```c++
bool enumeratePathCombinationsInFamily(
    PathFamily const& family, list<Path> const& segments, std::list<Path>* paths, size_t maxPaths);
```

Given a set of path segments from a family, this function
performs restricted enumeration of all paths that can be
obtained by merging path segments.

### Generate a set of maximal paths

```c++
bool getMaximalPathsForFamily(
    PathFamily const& family, list<Path>* paths, size_t maxPaths);
```

This function combines segment generation and path enumeration into a
single function call:

* It generates all path segments for a path family.
* It then enumerates all (or a maximum number of ) merged paths from
  these segments.
