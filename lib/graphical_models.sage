def create_function_graphical_model(G, same_color_list):
    """
    Creates the f_list for a Gaussian graphical model, where the functions
    are given by the adjugate matrix.

    Input:
        G (graph class of sage): graph
        same_color_list (list): list of prescribed color classes, each consisting
            either of vertices or of edges. A singleton class may be
            included explicitly, e.g. (k,) or ((u,v),), but this is optional:
            any vertex or edge not appearing in same_color_list is treated as having
            its own color automatically.
            If you use a coloring produced by the other functions in this file,
            combine the vertex color classes and edge color classes into a single list
            before passing it to this function.

That would make the convention completely explicit.
    Output:
        list: list of functions in QQ[k_i_j]
    """
    validate_same_color_list(G, same_color_list)
    vertices = list(G.vertices())
    n = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}
    edges = list(G.edges(labels=False))

    adjacency_mat = adjacency_matrix_from_labels(vertices, edges)
    K_mats_indices = []

    for same_color in same_color_list:
        mat = K_matrix_same_color_from_labels(vertices, same_color)

        first = same_color[0]
        if first in idx:
            i = idx[first]
            indices = (i, i)
        elif isinstance(first, (list, tuple)) and len(first) == 2:
            u, v = first
            i, j = idx[u], idx[v]
            if i > j:
                i, j = j, i
            indices = (i, j)
        else:
            raise ValueError("Invalid color class format.")

        K_mats_indices.append((mat, indices))
        adjacency_mat = adjacency_mat - mat

    K_mats_indices += K_mats_indices_given_adj_mat(adjacency_mat)

    K_vars = [var(f'k_{i}_{j}') for (_, (i, j)) in K_mats_indices]

    R = PolynomialRing(QQ, names=K_vars)
    K_gens = R.gens()

    K_mat = zero_matrix(R, n)
    for i in range(len(K_mats_indices)):
        K_mat += K_gens[i] * K_mats_indices[i][0]

    Sigma_list = []
    for i in range(n):
        for j in range(i + 1):
            Sigma_list.append((-1)^(i + j) * first_minor(K_mat, i, j))

    return Sigma_list
        
def validate_same_color_list(G, same_color_list):
    """
    Validate that same_color_list is a list/tuple of nonempty, pairwise disjoint
    color classes, where each class consists entirely of vertices or entirely of
    edges of G.

    Vertices and edges may be omitted from same_color_list; omitted elements are
    treated later as singleton color classes. It is also allowed to include
    singleton color classes explicitly.
    """
    vertices = set(G.vertices())
    edges = set((u, v) if u <= v else (v, u) for u, v in G.edges(labels=False))

    if not isinstance(same_color_list, (list, tuple)):
        raise TypeError("same_color_list must be a list or tuple of color classes.")

    seen_vertices = set()
    seen_edges = set()

    for same_color in same_color_list:
        if not isinstance(same_color, (list, tuple)):
            raise TypeError(f"Each color class must be a list or tuple, got {same_color}.")

        if len(same_color) == 0:
            raise ValueError("Color classes must not be empty.")

        first = same_color[0]

        # Vertex color class
        if first in vertices:
            for t in same_color:
                if t not in vertices:
                    raise ValueError(f"Mixed or invalid color class: {same_color}")
                if t in seen_vertices:
                    raise ValueError(f"Vertex {t} appears in more than one color class.")
                seen_vertices.add(t)

        # Edge color class
        elif isinstance(first, (list, tuple)) and len(first) == 2:
            for t in same_color:
                if not (isinstance(t, (list, tuple)) and len(t) == 2):
                    raise ValueError(f"Mixed or invalid color class: {same_color}")
                u, v = t
                e = (u, v) if u <= v else (v, u)
                if e not in edges:
                    raise ValueError(f"Edge {t} is not an edge of G.")
                if e in seen_edges:
                    raise ValueError(f"Edge {e} appears in more than one color class.")
                seen_edges.add(e)

        else:
            raise ValueError(f"Invalid color class: {same_color}")           

def adjacency_matrix_from_labels(vertices, edges):
    """
    Creates the adjacency matrix of a graph on the given vertex set.
    """
    vertices = list(vertices)
    n = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}

    A = zero_matrix(QQ, n)
    for i in range(n):
        A[i, i] = 1

    for u, v in edges:
        i, j = idx[u], idx[v]
        A[i, j] = 1
        A[j, i] = 1

    return A


def K_matrix_same_color_from_labels(vertices, same_color):
    """
    Creates the matrix corresponding to one color class of either vertices or edges.
    """
    vertices = list(vertices)
    n = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}

    mat = zero_matrix(QQ, n)

    if not same_color:
        raise ValueError("Color class must not be empty.")

    first = same_color[0]

    if first in idx:
        # vertex color class
        for t in same_color:
            if t not in idx:
                raise ValueError(f"Mixed or invalid color class: {same_color}")
            i = idx[t]
            mat[i, i] = 1

    elif isinstance(first, (list, tuple)) and len(first) == 2:
        # edge color class
        for t in same_color:
            if not (isinstance(t, (list, tuple)) and len(t) == 2):
                raise ValueError(f"Mixed or invalid color class: {same_color}")
            u, v = t
            i, j = idx[u], idx[v]
            mat[i, j] = 1
            mat[j, i] = 1

    else:
        raise ValueError(f"Invalid element in color class: {first}")

    return mat

def K_mats_indices_given_adj_mat(A):
    """
    Given an adjacency-style matrix A, create the remaining basis matrices
    and their index pairs.
    """
    n = A.nrows()
    K_mats_indices = []

    for i in range(n):
        for j in range(i, n):
            if A[i, j] == 1:
                mat = zero_matrix(QQ, n)
                mat[i, j] = 1
                mat[j, i] = 1
                K_mats_indices.append((mat, (i, j)))

    return K_mats_indices
    
def first_minor(A, j, k):
    """
    Computes the jk minor of a matrix A

    Input:
        A (matrix): matrix
        j (int): row index of the minor
        k (int): column index of the minor

    Output:
        scalar: the (j,k)-minor of A
    """

    submatrix = A.delete_rows([j]).delete_columns([k])     

    return submatrix.det()
    
    
def non_triangle_regular_blockgraph_colorings_up_to_automorphism(G):
    """
    Given a graph G, computes all possible colorings with distinct vertex and edge colors.
    Colorings are computed up to isomorphisms.
    If the graph G is a block graph, then triangle-regular colorings are omitted, as we know that these give toric ideals.
    Note: This only works for smaller graphs, as the number of colorings becomes too large.

    Input:
        graph: a graph G, see class https://doc.sagemath.org/html/en/reference/graphs/index.html

    Output:
        list: list of tuples, where each tuple two tuples, one with tuples of vertices having the same color and the second one with tuples of edges with the same color. 
                for example, the list can contain the element (((0, 1), (2,)), (((0, 1), (0, 2), (1, 2)),)), then the vertices 0, 1 have the same color and all edges have the same color.
    """
    if not G.is_block_graph():
        return colorings_up_to_automorphism(G)
    output = []
    for t in colorings_up_to_automorphism(G):
        if not vertex_regular(G,t[0],t[1]) or not edge_regular(G,t[0],t[1]) or not edge_triangle_regular(G,t[0],t[1]):
            output.append(t)
    return output
    

def colorings_up_to_automorphism(G):
    """
    Given a graph G, computes all possible colorings with distinct vertex and edge colors.
    Colorings are computed up to isomorphism.
    Note: This only works for smaller graphs, as the number of colorings becomes too large.

    Input:
        graph: a graph G, see class https://doc.sagemath.org/html/en/reference/graphs/index.html

    Output:
        list: list of tuples, where each tuple two tuples, one with tuples of vertices having the same color and the second one with tuples of edges with the same color. 
                for example, the list can contain the element (((0, 1), (2,)), (((0, 1), (0, 2), (1, 2)),)), then the vertices 0, 1 have the same color and all edges have the same color. 
    """
    vertex_partition = vertexcolorings_up_to_automorphism(G)

    known = set()        # set of normalized partitions we've seen
    reps = []            # one representative per orbit

    for pe in SetPartitions(G.edges(labels=False)):
        pe_blocks = [set(samecolor) for samecolor in pe]

        pe_norm = normalize_partition(pe_blocks)

        for pv in vertex_partition:
            is_equivalent = False
            
            for g in G.automorphism_group():
                im = (partition_vertex_perm(pv,g) , partition_edge_perm(pe_norm, g))
                if im in known:
                    is_equivalent = True
                    break
            if not is_equivalent:
                reps.append((pv,pe_norm))
                for g in G.automorphism_group():
                    known.add((partition_vertex_perm(pv,g) , partition_edge_perm(pe_norm, g)))
    
    return reps


def normalize_partition(p):
    """
    Normalize a partition of edges/vertices so it becomes a canonical tuple-of-tuples.
    If edges are given as triples (u,v,label),
    the label is dropped, leaving just (u,v) with u<=v.
    """
    blocks = []
    for block in p:
        clean_block = []
        for e in block:
            if isinstance(e, tuple) and len(e) == 3:
                u, v, _ = e   # drop the label
                if u > v:     # make orientation canonical
                    u, v = v, u
                clean_block.append((u, v))
            else:
                clean_block.append(e)
        blocks.append(tuple(sorted(clean_block)))
    return tuple(sorted(blocks))

def partition_vertex_perm(pv, perm):
    """
    Apply permutation `perm` (a permutation object from G.automorphism_group())
    to each element of the partition and return the normalized tuple form.
    """
    mapped_blocks = [tuple(sorted(perm(v) for v in block)) for block in pv]
    return tuple(sorted(mapped_blocks))

def partition_edge_perm(pe, g):
    """
    Apply permutation `perm` (a permutation object from G.automorphism_group())
    to each element of the partition and return the normalized tuple form.
    """
    mapped_blocks = [tuple(sorted( (g(e[0]) , g(e[1])) if g(e[0]) < g(e[1]) else (g(e[1]) , g(e[0])) for e in block)) for block in pe]
    return tuple(sorted(mapped_blocks))

def vertexcolorings_up_to_automorphism(G):
    known = set()        # set of normalized partitions we've seen
    reps = []            # one representative per orbit

    for p in SetPartitions(G.vertices()):
        p_blocks = [set(samecolor) for samecolor in p]
        p_norm = normalize_partition(p_blocks)

        # Check whether any automorphic image is already known
        is_equivalent = False
        for g in G.automorphism_group():
            im = partition_vertex_perm(p_blocks, g)
            if im in known:
                is_equivalent = True
                break
        if not is_equivalent:
            reps.append(p_norm)
            for g in G.automorphism_group():
                known.add(partition_vertex_perm(p_blocks, g))

    return reps


    

    
def vertex_regular(G, vertex_partition, edge_partition):
    """
    Check whether a coloring is vertex-regular.

    For each vertex, record how many incident edges it has of each edge color.
    Then vertices in the same vertex color class must have the same such profile.
    """
    n_edge_colors = len(edge_partition)

    # For each vertex, store counts of incident edge colors
    nb_colors_vertices = {
        v: [0 for _ in range(n_edge_colors)]
        for v in G.vertices()
    }

    for j, edge_color_class in enumerate(edge_partition):
        for e in edge_color_class:
            u, v = e[0], e[1]
            nb_colors_vertices[u][j] += 1
            nb_colors_vertices[v][j] += 1

    for vertex_color_class in vertex_partition:
        if len(vertex_color_class) <= 1:
            continue

        reference_profile = nb_colors_vertices[vertex_color_class[0]]
        for v in vertex_color_class[1:]:
            if nb_colors_vertices[v] != reference_profile:
                return False

    return True


def edge_regular(G, vertex_partition, edge_partition):
    """
    Check whether a coloring is edge-regular.

    For an edge e = (u,v), its profile is the pair of endpoint vertex colors,
    counted as a vector over vertex color classes.

    Since each endpoint belongs to exactly one vertex color class, this profile
    records how many endpoints of e lie in each vertex color class.

    Then edges in the same edge color class must have the same such profile.
    """
    n_vertex_colors = len(vertex_partition)

    # Map each vertex to its vertex color index
    vertex_color_of = {}
    for j, vertex_color_class in enumerate(vertex_partition):
        for v in vertex_color_class:
            vertex_color_of[v] = j

    # For each edge, store how many endpoints lie in each vertex color class
    edge_profiles = {}

    for e in G.edges(labels=False):
        u, v = e
        profile = [0 for _ in range(n_vertex_colors)]
        profile[vertex_color_of[u]] += 1
        profile[vertex_color_of[v]] += 1

        # Store under canonical orientation
        edge_key = (u, v) if u <= v else (v, u)
        edge_profiles[edge_key] = profile

    for edge_color_class in edge_partition:
        if len(edge_color_class) <= 1:
            continue

        e0 = edge_color_class[0]
        e0_key = (e0[0], e0[1]) if e0[0] <= e0[1] else (e0[1], e0[0])
        reference_profile = edge_profiles[e0_key]

        for e in edge_color_class[1:]:
            e_key = (e[0], e[1]) if e[0] <= e[1] else (e[1], e[0])
            if edge_profiles[e_key] != reference_profile:
                return False

    return True


def edge_triangle_regular(G, vertex_partition, edge_partition):
    """
    Check whether a coloring is edge-triangle-regular.

    For each edge e = (u,v), look at all vertices k such that {u,v,k} forms a triangle.
    For every such triangle, record the pair of colors of the two other edges:
        (u,k) and (v,k),
    sorted so that the order does not matter.

    Then edges in the same edge color class must have the same multiset of such pairs.
    """
    n_edge_colors = len(edge_partition)

    # Map each edge to its edge color index
    edge_color_of = {}
    for j, edge_color_class in enumerate(edge_partition):
        for e in edge_color_class:
            u, v = e[0], e[1]
            edge_key = (u, v) if u <= v else (v, u)
            edge_color_of[edge_key] = j

    def canonical_edge(u, v):
        return (u, v) if u <= v else (v, u)

    def triangle_profile(u, v):
        """
        Return the sorted list of color-pairs arising from triangles containing edge (u,v).
        """
        profile = []
        Nu = set(G.neighbors(u))
        Nv = set(G.neighbors(v))

        # Common neighbors give triangles
        for k in Nu.intersection(Nv):
            if k == u or k == v:
                continue

            color1 = edge_color_of[canonical_edge(u, k)]
            color2 = edge_color_of[canonical_edge(v, k)]

            profile.append(tuple(sorted((color1, color2))))

        profile.sort()
        return profile

    for edge_color_class in edge_partition:
        if len(edge_color_class) <= 1:
            continue

        e0 = edge_color_class[0]
        reference_profile = triangle_profile(e0[0], e0[1])

        for e in edge_color_class[1:]:
            if triangle_profile(e[0], e[1]) != reference_profile:
                return False

    return True


def colorings_at_most_two_colors_for_vertices_and_two_for_edges(list_of_colorings):
    return [col for col in list_of_colorings if len(col[0]) <= 2 and len(col[1]) <= 2]



def one_edge_color_one_vertex_color(G):
    vertices = tuple(sorted(G.vertices()))
    edges = tuple(sorted((u, v) if u <= v else (v, u) for u, v in G.edges(labels=False)))
    return ((vertices,), (edges,))