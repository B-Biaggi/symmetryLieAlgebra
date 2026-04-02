# COMPUTING THE CONTINUOUS SYMMETRIES OF A PARAMETRIZED VARIETY

SageMath code for Algorithms presented in the paper https://arxiv.org/abs/????.????
The algorithms computes the symmetry Lie algebra of a parametrized variety. Using the symmetry Lie algebra,
there is second algorithm which decides wheter the vanishish ideal of the parametrized variety is biniomial after a linear
change of coordinates. 
The paper contains computational experiments on staged tree models and colored Gaussian graphical models. In both cases, we include
the code to create the parametrization and a Jupyter notebook showing our computations.

There are four different files containing code:
- `symmetry_Lie_algebra.sage` contains Algorithm 1 computing the symmetry Lie algebra.
- `binomiality.sage` contains Algorithm 2 testing binomilyt of the vanishing ideal.
- `staged_trees.sage` contains the code to create the parametrization of a colored Gaussian graphical model. 
- `graphical_models.sage` contains the code to create the parametrization of a colored Gaussian graphical model.


The following Jupyter notebooks contain the examples of Section 6:
- `binary_two-staged_tree_models.ipynb` contains the computation done in Example 6.1. Decides for all 127 two state binary tree models with depth at most 3 if they are binomial or not.
- `caterpillar_tree_expansion.ipynb` contains the computation done in Example 6.2. Decides for all expansions of the caterpillar tree at one or two leaves if it is biniomial or not.
- `K4_colorings.ipynb` contains the computation done in Example 6.3. Decides for all 215 colorings of the complete graph K<sub>4</sub> if they define a binomial model.
- `C4_colorings.ipynb` contains the computation done in Example 6.4. Decides for all 74 colorings of the cycle C<sub>4</sub> if they define a binomial model.'


## Documentation

### symmetry\_Lie\_algebra.sage
```
symmetry_lie_algebra_random(flist , epsilon = 0)
```

This is Algorithm 1 computing the symmetry Lie algebra of a parametrized variety. The algorithm follows the ideas described there and in Remark 3.3: 

There are two options for the sampling of points. Either, one can enter an epsilon as described in Proposition 3.1 or sample from a smaller set of points [0,1,-1,2,-2,1/2,-1/2, 1/3, 2/3, -1/3, -2/3, 3,-3]. 

After the setup, we first try to guess the rank of the Jacobian at random points. We take the maximum rank at 30 random points.

The algorithm works as follows. We start by setting L to be the full space of n x n matrices. In every step, we randomly generate a point p, evaluate the image of the point using the rational maps f\_i and compute the Jacobian J at this point. We denote the image of p under f by q. If the Jacobian has smaller rank than expected, we start over again. Else we solve for which matrices A in L, it holds that Aq is an element of the image of J.

We repeat the while loop either n^2 times or end if all elements A in L satisfy that Aq is an element of the image of J. In the second case, we repeat 5 times before we end the loop.

This way, we do not necessary have to repeat the loop n^2 times. Compared to Algorithm 1, we solve the linear equation system for L in every step instead of only in the end. 

Input:
- flist (list): A list of rational functions over QQ[x\_1,..., x\_m]
- epsilon (rational, optional): Enter epsilon in (0,1) to have probability at leat 1-epsilon (see Prop 3.1). If left empty reduce random choices to [0,1,-1,2,-2,1/2,-1/2, 1/3, 2/3, -1/3, -2/3, 3,-3]

Output: 
- list of Matrices: A basis of the space of the symmetry Lie algebra


### binomiality.sage
```
test_binomiality_random(flist , B,  epsilon = 0)
```
This is Algorithm 2 computing the symmetry Lie algebra of a parametrized variety. It proceeds as described. Needs as an input the output from `symmetry_lie_algebra_random`.  

Note that the semisimple part is computed using Hensel lift.

Input: 
- flist (list): A list of rational functions over QQ[x_1,...,x_m]
- B (list of Matrices): Input the basis of the symmetry Lie algebra given by symmetry\_lie\_algebra\_random
- epsilon (rational, optional): Enter epsilon in (0,1) to have probability as in the paper
        
Output:
- boolean

### graphical\_models.sage
This contains different funtions to create Gaussian graphical models and the paramatrization.

We describe colored graphs as follows: Each graph has a number n, the number of vertices. Additionally, we need the list of edges and a list of list, where each list contains vertices of edges of the same color. For example: 

n=3, \
edges = [[1,2],[2,3]], \
 same\_color\_list= [[1,3],[[1,2],[2,3]]

 describes a graph with three vertices 1,2,3, edges \{1,2\} and \{2,3\} and where the vertices 1,3 have the same color and the edges \{1,2\}, \{2,3\} have the same color. In total, this example has 3 different colors. 

```
create_function_graphical_model(n,edges,same_color_list)
```
Creates the f_list for a Gaussian graphical model, where the functions are given by the adjugate matrix.

Some other usefull function to experiment are:

```
non_triangle_regular_blockgraph_colorings_up_to_automorphism(G)
```
Given a graph G, computes all possible colorings with distinct vertex and edge colors. Colorings are computed up to isomorphisms. If the graph G is a block graph, then triangle-regular colorings are omitted, as we know that these give binomial ideals (see https://arxiv.org/abs/2507.06437).

```
colorings_up_to_automorphism(G)
```
Given a graph G, computes all possible colorings with distinct vertex and edge colors. Colorings are computed up to isomorphism.





### staged\_tree.sage

We encode staged trees as follows: The tree needs two information, a dictionary containing the list of internal nodes pointing to their data, which is a dictionary containing their stage. The nodes are in a tuple format describing the path tho the root, so for example (1,) is the root or (1,2,1) is the note with depth 2 which is the first child or the internal node (1,2), which is itself the second child of the root (1,2). Secondly, a tree needs the dictionary stage\_children, which contains for all stages how many outgoing edges branch at the given internal node. 

For example, the following is the Caterpillar tree (Table 2): \
caterpillar = {
    (1,): {'stage':1},
    (1,1): {'stage':1},
    (1,1,2): {'stage':1},
    (1,1,2,3): {'stage':1},
} \
stage_children = {1:3} 
```
tree_polynomials_hom(tree, stage_children, root=(1,))
```
Computes the homogenous polynomial of a staged three. 

```
trees_up_to_symmetry(max_depth, stage_children, stages, quotient_stage_labels=False)
```
Generate all staged rooted trees up to child-permutation symmetry
with depth <= max_depth.
If input 'quotient_stage_labels=True', it also gives trees up to symmetry under permuting stages with the same number of children. 

```
expand_tree_k_times(tree, stage_children, stages, k, quotient_stage_labels=False):
```
Expand tree exactly k times, keeping depth <= original depth,
eliminating symmetric duplicates.

```
keep_only_trees_with_all_stages_at_least(trees, stages, at_least_many_times)
```
Filter trees so that every stage in `stages` appears at least `at_least_many_times` times.



## Examples

### Example 6.1 (Binary two-staged tree models)
The Jupyter notebook `binary_two-staged_tree_models.ipynb` contains the computation of Example 6.1. 
Using our setup, there exist two stages [1,2] and both stages have two childre, so stages_children = {1:2, 2:2}. 
We compute all binary trees with depth at most 3, so we use 'trees\_up\_to\_symmetrie(4, stages\_children , stages)'. 
For each of the trees, we then compute wheter the model is binomial or not.

### Example 6.2 (Expansions of the caterpillar tree)
The Jupyter notebook `caterpillar_tree_expansion.ipynb` contains the computation of Example 6.2. 

In our setting, the Caterpillar tree is described by

caterpillar = {
    (1,): {'type':1},
    (1,1): {'type':1},
    (1,1,2): {'type':1},
    (1,1,2,3): {'type':1},
}

stage\_children = {1:3}

stages = [1]

Using the function `expand_tree_k_times`, we create all expansion at one, two or three leaves and decide whether the model is binomial or not.

### Example 6.3 (Colorings of the complete graph K<sub>4</sub>)

We create all colorings of the complete graph with 4 vertices of which we know that the vanishing ideal is not binomial (without a linear change of coordinates.)

There are 215 colorings. For all of those, we compute the symmetry Lie algebra and decide whether the vanishing ideal is binomial. 

### Example 6.4 (Colorings of the cycle C<sub>4</sub>)

We create all colorings of the cycle with 4 vertices.

There are 61 colorings. For all of those, we compute the symmetry Lie algebra and decide whether the vanishing ideal is binomial. 

## Tests

In the folder `test`, there are some tests regarding the code:

`test_uncolored_graphs_4_vertices.ipynb` \
Compute the symmetry Lie algebra of all connected uncolored graphs with 4 vertices and compares the result to the ones of Kahle-Vill in https://arxiv.org/pdf/2408.14323 p. 13.

`test_uncolored_graphs_5_vertices.ipynb` \
Compute the symmetry Lie algebra of all connected uncolored graphs with 5 vertices and compares the result to the ones of Kahle-Vill in https://arxiv.org/pdf/2408.14323 p. 14.

`test_trees_up_to_symmetry.ipynb`\
Some tests regarding the function `trees_up_to_symmetry` to check whether the ones listed are the ones expected.

`test_expand_tree_k_times`\
Some tests regarding the function `expand_tree_k_times` to check whether the ones listed are the ones expected.

`test_colorings`\
Some tests regarding the function `colorings_up_to_automorphism` and `non_triangle_regular_blockgraph_colorings_up_to_automorphism` to check whether they create the expected number of colorings.