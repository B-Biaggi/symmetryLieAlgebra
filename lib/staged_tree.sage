from itertools import permutations, product
from copy import deepcopy

"""
We encode staged trees as follows: The tree needs two information, a dictionary containing the list of internal nodes pointing to their data, which is a dictionary containing their stage. The nodes are in a tuple format describing the path tho the root, so for example (1,) is the root or (1,2,1) is the note with depth 2 which is the first child or the internal node (1,2), which is itself the second child of the root (1,2). Secondly, a tree needs the dictionary stage_children, which contains for all stages how many outgoing edges branch at the given internal node. 

For example, the following is the Caterpillar tree (Table 2): 
caterpillar = {(1,): {'stage':1},
    (1,1): {'stage':1},
    (1,1,2): {'stage':1},
    (1,1,2,3): {'stage':1},} 
stage_children = {1:3} 
"""

def tree_polynomials_hom(tree, stage_children, root=(1,)):
    """
    Computes the homogenous polynomial to a staged three. 

    Input:
        tree (dictonary): describes all the internal nodes of the tree and their stage
        stage_children (dictonary): number of branches of a node of a given stage
    Ouptut:
        list: list of homogenous polynomials, see Section 6
    """
    R, gens = polynomial_ring_hom(stage_children)
    paths = all_paths(tree, stage_children, root)
    # total degree = maximum path length - 1 (number of edges)
    total_degree = max(len(p)-1 for p in paths)
    polys = [path_polynomial_hom(p, tree, stage_children, R, gens, total_degree) for p in paths]
    return polys


def trees_up_to_symmetry(max_depth, stage_children, stages, quotient_stage_labels=False):
    """
    Generate all staged rooted trees up to child-permutation symmetry
    with depth <= max_depth.

    Input:
        max_depth (int): Maximum depth of an internal node, where the root (1,) has depth 1.
        stage_children (dict): Dictionary {stage: integere denoting how many children}.
        stages (list): Allowed stage labels for internal nodes.
        quotient_stage_labels (bool, optional): If True, also quotient by permutations of stage labels
            via canonical_form(tree, stage_children, stages).
            If False, only quotient by child symmetry.

    Output:
        list of dict: Trees in dictionary representation.
    """

    def cf(tree):
        if quotient_stage_labels:
            return canonical_form(tree, stage_children, stages)
        else:
            return canonical_form(tree, stage_children)

    # Start with all possible one-node trees
    all_seen = {}
    frontier = []

    for s in stages:
        tree = {(1,): {'stage': s}}
        key = cf(tree)
        if key not in all_seen:
            all_seen[key] = tree
            frontier.append(tree)

    # Breadth-first generation by repeated leaf expansion
    while frontier:
        next_frontier = []

        for tree in frontier:
            # expand_one_leaf only adds internal nodes at leaves whose parent
            # has depth < max_depth, so resulting internal nodes have depth <= max_depth
            expanded_trees = expand_one_leaf(tree, stage_children, stages, max_depth)

            for new_tree in expanded_trees:
                key = cf(new_tree)
                if key not in all_seen:
                    all_seen[key] = new_tree
                    next_frontier.append(new_tree)

        frontier = next_frontier

    return list(all_seen.values())

def expand_tree_k_times(tree, stage_children, stages, k, quotient_stage_labels=False):
    """
    Expand tree exactly k times, keeping depth <= original depth,
    eliminating symmetric duplicates.
    """
    original_depth = tree_depth(tree)

    def cf(tr):
        if quotient_stage_labels:
            return canonical_form(tr, stage_children, stages)
        else:
            return canonical_form(tr, stage_children)

    current_set = {cf(tree): tree}

    for _ in range(k):
        next_set = {}
        for tr in current_set.values():
            expansions = expand_one_leaf(tr, stage_children, stages, original_depth)
            for new_tr in expansions:
                key = cf(new_tr)
                if key not in next_set:
                    next_set[key] = new_tr
        current_set = next_set

    return list(current_set.values())




def keep_only_trees_with_all_stages_at_least(trees, stages, at_least_many_times):
    """
    Filter trees so that every stage in `stages` appears at least `at_least_many_times` times.

    Input: 
        trees (list of dict): Trees in dictionary representation.
        stages (list of int): stages that must be present.
        at_least_many_times (int): Minimum number of occurrences required for EACH stage.

    Output:
    list of dict: Trees satisfying the condition.
    """
    required = set(stages)
    k = at_least_many_times

    filtered = []
    for tree in trees:
        # Count occurrences of each stage
        counts = {}
        for node in tree:
            t = tree[node]['stage']
            counts[t] = counts.get(t, 0) + 1

        # Check: every required stage appears at least k times
        if all(counts.get(t, 0) >= k for t in required):
            filtered.append(tree)

    return filtered


# -----------------------------
# Helper functions to create the polynomials:
# -----------------------------


def all_paths(tree, stage_children, root=(1,)):
    """
    Compute all root->leaf paths including implicit leaf edges.
    """
    paths = []
    stack = [(root, [root])]

    while stack:
        node, path = stack.pop()

        t = tree[node]["stage"]
        num_edges = stage_children[t]

        # Generate all children (internal or leaf)
        children = [node + (i,) for i in range(1, num_edges + 1)]

        # Check which children are internal
        internal_children = [c for c in children if c in tree]

        # If no children are internal, all children are leaves
        if len(internal_children) == 0:
            # append a path for each leaf edge
            for leaf_child in children:
                paths.append(path + [leaf_child])
        else:
            # For each child, extend the path
            for c in children:
                if c in tree:
                    # internal child: continue recursion
                    stack.append((c, path + [c]))
                else:
                    # leaf child: end path here
                    paths.append(path + [c])

    return paths


def polynomial_ring_hom(stage_children):
    """
    Creates the polynomial ring over QQ with as many varialbes as stages.
    """
    names = []
    for t, m in stage_children.items():
        for j in range(1, m):  # only m-1 real variables
            names.append(f"x_{t}_{j}")
    names.append("z")  # homogenizing variable
    R = PolynomialRing(QQ, names=names)
    gens = R.gens_dict()
    return R, gens



def edge_variable(t, k, stage_children, gens):
    m = stage_children[t]
    if k < m:
        return gens[f"x_{t}_{k}"]
    else:
        return gens["z"] - sum(gens[f"x_{t}_{j}"] for j in range(1, m))


def path_polynomial_hom(path, tree, stage_children, R, gens, total_degree):
    """
    Converts path to homogeneous polynomial
    """
    poly = R(1)
    for i in range(len(path) - 1):
        parent = path[i]
        child = path[i+1]
        t = tree[parent]["stage"]
        edge = child[-1]
        poly *= edge_variable(t, edge, stage_children, gens)
    # homogenize by multiplying z^(total_degree - deg)
    deg = poly.degree()
    poly *= gens["z"]**(total_degree - deg)
    return poly



# -----------------------------
# Helper functions to create trees;
# -----------------------------

def canonical_form(tree, stage_children, stages=None):
    """
    Canonical form of a staged rooted tree under:
      - global child permutations for each stage
      - optionally, stage relabelings among stages with equal arity
    """

    # rigid encoding of a dictionary-tree
    def rigid_code(tr):
        return tuple(sorted((node, data['stage']) for node, data in tr.items()))

    # apply one global child permutation per stage
    def apply_child_permutations(tr, child_perm_map):
        new_tree = {}

        for node, data in tr.items():
            new_node = [node[0]]  # root starts with 1
            current_path = (node[0],)

            for edge in node[1:]:
                parent_stage = tr[current_path]['stage']
                perm = child_perm_map[parent_stage]
                new_edge = perm[edge]
                new_node.append(new_edge)
                current_path = current_path + (edge,)

            new_tree[tuple(new_node)] = {'stage': data['stage']}

        return new_tree

    # relabel stages
    def apply_stage_relabel(tr, stage_map):
        return {node: {'stage': stage_map[data['stage']]} for node, data in tr.items()}

    # allowed stage relabelings: only inside equal-arity groups
    if stages is None:
        stage_maps = [{t: t for t in stage_children}]
    else:
        groups = {}
        for t in stages:
            groups.setdefault(stage_children[t], []).append(t)

        group_keys = sorted(groups.keys())
        perm_lists = [list(permutations(groups[a])) for a in group_keys]

        stage_maps = []
        for choices in product(*perm_lists):
            mp = {}
            for a, perm in zip(group_keys, choices):
                grp = groups[a]
                for old, new in zip(grp, perm):
                    mp[old] = new
            stage_maps.append(mp)

    # allowed child permutations: one global permutation for each stage
    child_perm_lists = []
    stage_order = sorted(stage_children.keys())
    for t in stage_order:
        m = stage_children[t]
        perms_t = []
        for p in permutations(range(1, m + 1)):
            perm_dict = {i + 1: p[i] for i in range(m)}
            perms_t.append(perm_dict)
        child_perm_lists.append(perms_t)

    best = None

    for stage_map in stage_maps:
        relabeled_tree = apply_stage_relabel(tree, stage_map)

        for choices in product(*child_perm_lists):
            child_perm_map = {t: choices[i] for i, t in enumerate(stage_order)}
            transformed = apply_child_permutations(relabeled_tree, child_perm_map)
            code = rigid_code(transformed)

            if best is None or code < best:
                best = code

    return best


def tree_depth(tree):
    return max(len(node) for node in tree)


def expand_one_leaf(tree, stage_children, stages, original_depth):
    """
    Expand one implicit leaf into an internal node,
    assigning any allowed stage, only if depth does not exceed original_depth.
    """
    leaves = []

    for node in tree:
        node_stage = tree[node]['stage']
        arity = stage_children[node_stage]

        if len(node) >= original_depth:
            continue

        for i in range(1, arity + 1):
            child = node + (i,)
            if child not in tree:
                leaves.append(child)

    result = []
    for leaf in leaves:
        for s in stages:
            new_tree = deepcopy(tree)
            new_tree[leaf] = {'stage': s}
            result.append(new_tree)

    return result


