import random
import math
import ast
import pickle
import matplotlib.pyplot as plt
import time
from itertools import combinations, permutations
import itertools
NE, NW, SE, SW = 'NE', 'NW', 'SE', 'SW'
from sympy import Matrix
import numpy as np

#===============================================================================
#   Sum Tree
#===============================================================================


class SumTree:
    "An array with quick prefix/suffix sum, using a binary tree"

    def __init__(self, n):
        self.length = n
        self.arrays = [(2**i)*[0] for i in range(math.ceil(math.log(n,2))+1)][::-1]
        
    def add(self, i, val):
        for x in range(len(self.arrays)):
            self.arrays[x][int(i // 2**x)] += val  

    def sum_suffix(self, i):
        return sum([self.arrays[x][i//2**x + 1]
                    for x in range(len(self.arrays)-1) if i//2**x % 2 == 0]) 

    def sum_prefix(self, i):
        return sum([self.arrays[x][i//2**x - 1]
                    for x in range(len(self.arrays)-1) if i//2**x % 2 == 1]) 

#===============================================================================
#   Vertex-Edge Implementation
#===============================================================================



def vertex(permutation, tree):
    '''
        permutation : one-line encoding of the 'large' permutation
        tree : (
                 (EDGE_LABEL_1, CHILD_1),
                 ...
                 (EDGE_LABEL_k, CHILD_k) )
    '''
    n = len(permutation)
    A = [1] * n
    
    for edge_label, child in tree:
        child_result = edge(permutation, child, edge_label)
        A = [a * c for a, c in zip(A, child_result)]
    
    return A

def edge(permutation, tree, edge_label):
    '''
        permutation : one-line encoding of the 'large' permutation
        tree : (
                 (EDGE_LABEL_1, CHILD_1),
                 ...
                 (EDGE_LABEL_k, CHILD_k) )
        edge_label : edge label 
    '''
    n = len(permutation)
    A = vertex(permutation, tree)
    C = [0] * n
    B = SumTree(n)
    S, W = edge_label

    if W == 'W':
        for i in range(n):
            idx = permutation[i]
            if S == 'S':
                C[i] = B.sum_prefix(idx)
            else:
                C[i] = B.sum_suffix(idx)
            B.add(idx,A[i]) 
    else:
        for i in reversed(range(n)):
            idx = permutation[i]
            if S == 'S':
                C[i] = B.sum_prefix(idx)
            else:
                C[i] = B.sum_suffix(idx)
            B.add(idx,A[i]) 

    return C

#===============================================================================
#   Printing corner trees
#===============================================================================

def pretty_print_tree(tree, depth=0):
    children = tree
    indent = "    " * depth # 4 spaces per level

    for direction, subtree in children:
        print(f"{indent}{direction}")
        pretty_print_tree(subtree, depth + 1)

#===============================================================================
#   Corner trees: patterns expansion
#===============================================================================


def shape(c):
    sc = sorted(c)
    return [sc.index(x) for x in c]

def count_brute_force(pattern, permutation):
    count = 0
    for sigma in list(combinations(permutation, len(pattern))):
        if list(pattern) == shape(sigma):
            count += 1
    return count

def expand_in_patterns(cornertree, n_vertex):
    
    combos = []
    for i in range(1, n_vertex+1):

        combo = []
        for sigma in permutations(range(i)): 

            coef = sum(vertex(sigma,cornertree))

            for j in range(1,i):
     
                for t, tau in enumerate(permutations(range(j))):
                    
                    
                    coef -= count_brute_force(tau, sigma) * combos[j-1][t]
            combo.append(coef)
        combos.append(combo)
    return combos

def test_counting_corner(corner_tree, n_vertex):
    n = 20
    permutation = random.sample(range(n), n)
    x = sum(vertex(permutation,corner_tree))
    a = expand_in_patterns(corner_tree, n_vertex)
    y = 0
    for z in range(len(a)):
        for t, tau in enumerate(permutations(range(z+1))):
            y += a[z][t]*count_brute_force(tau,permutation)
    assert x == y
    

#===================================================================================================================
#   Algorithm for pure west corner trees occurrences, and type_A, type_A_not_B occurrences of Tree_{5/3} 
#===================================================================================================================

from collections import defaultdict

def hash_tree(tree):
    return hash(tree)

def countW(permutation, tree, row=None, chunk_size=None, not_B=False):
    '''
        Depending on the presence of 'row' and the value of 'not_B' the function
        counts:

            - (row=None) The array containing at each position i the number of
              occurences of tree (pure west corner tree) in the permutation where the root is
              mapped to (i,permutation(i)).

            - (row!=None, not_B=False)
              Here, tree is assumed to be an element of CornerTree_{5/3};
              it has two disinguished nodes called MN and ME.
              For us tree is rooted in ME!

              Let t be the tree double poset in Tree_{5/3} where a node "MAX" is attached to the north-east
              of 'tree'. The function returns an array containing at each position i
              the number of occurences of the double poset where MAX is mapped to
              (i, permutation(i)) IF permutation(i) is in the slice [row, row+chunk_size)
              (otherwise the entry is 0) and MN is mapped (in vertical direction) into [0,row).

            - (row!=None, not_B=True)
              Same as the (row!=None, not_B=False) case but with the additional constraint
              that ME is mapped (in horizontal direction) into [j * chunk_size, (j+1) * chunk_size),
              where j is such that v is mapped (in horizontal direction) into [j * chunk_size, (j+1) * chunk_size).

        permutation : the 'large' permutation
    '''

    A_v_dict = {}
    B_e_dict = {}
    touched_dict = defaultdict(bool)

    def vertex_W(permutation, tree, i):
        if not tree in A_v_dict:
            A_v_dict[tree] = len(permutation) * [1]
        if not touched_dict[(tree, i)]:
            touched_dict[(tree, i)] = True
            for label, child in tree:
                result = edge_W(permutation, child, label, i)
                A_v_dict[tree][i] = A_v_dict[tree][i]*result

        return A_v_dict[tree][i]

    def edge_W(permutation, child, edge_label, i):
        edge_identifier = (child, edge_label)

        if not edge_identifier in B_e_dict:
            B_e_dict[edge_identifier] = SumTree(len(permutation))

        idx = permutation[i]
        if not touched_dict[(edge_identifier, i)]:
            touched_dict[(edge_identifier, i)] = True
            B_e_dict[edge_identifier].add(idx, vertex_W(permutation, child, i))
        if edge_label[0] == "S":
            return B_e_dict[edge_identifier].sum_prefix(idx)
        else:
            return B_e_dict[edge_identifier].sum_suffix(idx)
           
    ret = []
    if not_B:
        countCT = 0
        for i in range(len(permutation)):
            y = permutation[i]
            if i % chunk_size == 0:
                countCTback = countCT
            if y < row:
                countCT += vertex_W(permutation, tree, i)
            if row <= y < row + chunk_size:
                ret.append( countCT-countCTback )
            else:
                ret.append(0)
    elif row is not None:
        countCTW = 0
        for i in range(len(permutation)):
            if permutation[i] < row:
                countCTW += vertex_W(permutation, tree, i)
            if row <= permutation[i] < row + chunk_size:
                ret.append( countCTW )
            else:
                ret.append(0)
    else:
        for i in range(len(permutation)):
            ret.append( vertex_W(permutation, tree, i) )
    return ret
    

def test_counting_corner_pure_west(corner_tree, n_vertex):
    n = 12
    permutation = random.sample(range(n), n)
    x = sum(countW(permutation, corner_tree, row=None, chunk_size=None, not_B=False))
    a = expand_in_patterns(corner_tree, n_vertex)
    y = 0
    for z in range(len(a)):
        for t, tau in enumerate(permutations(range(z+1))):
            y += a[z][t]*count_brute_force(tau,permutation)
    assert x == y

#===============================================================================
#   Product Tree
#===============================================================================

def dyadic_range(l, r):
    x = 1
    while l != r:
        if l % (2*x) != 0:
            yield (l,l+x)
            l += x
        if r % (2*x) != 0:
            yield (r-x,r)
            r -= x
        x *= 2

class ProductTree:
    "A 2-dim array with quick box sum"
    
    def __init__(self, n):
        self.length = n
        self.logn = int(math.ceil(math.log(max(1,n))/math.log(2)))
        self.table = {}

    def add(self, x, y, value):
        for i in range(self.logn + 1):
            for j in range(self.logn + 1):
                key = (x-x%2**i, x-x%2**i+2**i,
                       y-y%2**j, y-y%2**j+2**j)
                self.table[key] = self.table.get(key,0) + value



    def sum_box(self, x1, x2, y1, y2):
        count = 0
        for x_range in dyadic_range(x1, x2):
            for y_range in dyadic_range(y1, y2):
                count += self.table.get(x_range + y_range, 0)
        return count    

def invperm(p):
    d = {y:x for x,y in enumerate(p)} 
    return [d[y] for y in range(len(p))]


#===================================================================================================================
#   Algorithm for type_not_A_not_B occurrences of Tree_{5/3} 
#===================================================================================================================


def count_Box(permutation,left_trees,middle_tree,right_trees):
     
    "left trees are the corner trees dangling from ME,"
    " middle tree is the corner trees with root between ME and MN"
    " right trees are the corner trees dangling from MN"
    
    n = len(permutation)
    
    chunk = int(n**(1/3.))
    permT = invperm(permutation)
    count = 0
   
    
        
    middle_box = ProductTree(n)
    right_boxes = [ProductTree(n) for _ in right_trees]
    left_boxes = [ProductTree(n) for _ in left_trees]

    for i in range(len(left_trees)):
        outcome = countW(permutation, left_trees[i], row=None, chunk_size=None, not_B=False)
        for x,y in enumerate(permutation):
            left_boxes[i].add(x,y,outcome[x])
    for i in range(len(right_trees)):
        outcome = countW(permutation, right_trees[i], row=None, chunk_size=None, not_B=False)
        for x,y in enumerate(permutation):
            right_boxes[i].add(x,y,outcome[x])
    outcome = countW(permutation, middle_tree, row=None, chunk_size=None, not_B=False)
    for x,y in enumerate(permutation):
            middle_box.add(x,y,outcome[x])
   
        
    
    for x4,y4 in enumerate(permutation):
        col,row = x4-x4%chunk,y4-y4%chunk
        for x1,y1 in enumerate(permutation[col:x4], col):
            for y3,x3 in enumerate(permT[row:y4], row):
                if x3 < x1 and y3 > y1:
                  
                    left_product = right_product = 1
                    for i in range(len(left_boxes)):
                        left_product*= left_boxes[i].sum_box(0,x1,0,y1)

                    for i in range(len(right_boxes)):
                        right_product *= right_boxes[i].sum_box(0,x3,0,y3)
                    
                    count =count+ left_product*right_product*middle_box.sum_box(x3+1,x1,y1+1,y3) 
                    
    
    return count

#===================================================================================================================
#   Algorithm for occurrences of Tree_{5/3} 
#===================================================================================================================

def count_gen(permutation,tree,inverse_tree,left_trees,middle_tree,right_trees):
    "tree is the corner tree rooted in ME, once we remove MAX from the element of Tree_{5/3} of which we count the occurrences"
    "inverse_tree is the corner tree obtained after swapping and exchanging labels of ME and MN"
    "left trees are the corner trees dangling from ME,"
    " middle tree is the corner trees with root between ME and MN"
    " right trees are the corner trees dangling from MN"
    n = len(permutation)
    chunk = int(n**(1/3.))
    permT = invperm(permutation)
    count = 0
    for row in range(0, n, chunk):
        count += sum(countW(permutation, tree, row=row, chunk_size=chunk, not_B=False))
    for col in range(0, n, chunk):
        count += sum(countW(permT, inverse_tree, row=col, chunk_size=chunk, not_B=True))
    count += count_Box(permutation,left_trees,middle_tree,right_trees)
     
    return count


def test_counting_tree_five_third():
    tree = ((NW, ((NW, ((SW,()),)),(SW, ()))),(SW, ((SW, ()),(SW, ()))))
    n_vertex_tree = 8
    inverse_tree = ((SW,()),(NW,((SW, ()),(NW, ((SW,((SW, ()),(SW, ()))),)))))
    pretty_print_tree(tree)
    pretty_print_tree(inverse_tree)
    left_trees = [((SW,()),(SW,()))]
    middle_tree =((SW,()),)
    right_trees = [()]
    n = 20
    permutation = random.sample(range(n), n)
    y = count_gen(permutation,tree,inverse_tree,left_trees,middle_tree,right_trees)
    a = expand_in_patterns(tree, n_vertex_tree)
    x = 0
    for z in range(len(a)):
        for t, tau in enumerate(permutations(range(z+1))):
            x += a[z][t]*count_brute_force(tau+(len(tau),),permutation)
    print(x,y)
    assert x == y

#===================================================================================================================
#   Enumerating isomorphism classes of corner trees with n vertices without repetition
#===================================================================================================================


def adj_pairs(a):
    return zip(a, a[1:])

def partitions(n):
    if n <= 0:
        yield []
        return
    for x in itertools.product([False, True], repeat=n-1):
        bars = [0] + [i for i, val in enumerate(x, 1) if val] + [n]
        yield [y - z for z, y in adj_pairs(bars)]

def canonical_corner_trees(n):
    if n == 1:
        return [()]
    
    result = []
    directions = ['NE', 'NW', 'SE', 'SW']
    
    for p in partitions(n-1):
        for tuple_cts in itertools.product(*[canonical_corner_trees(q) for q in p]):
            for a in itertools.product(directions, repeat=len(p)):
                tree = [((a[j], tuple_cts[j]),) if tuple_cts[j] else ((a[j], ()),) for j in range(len(p))]
                
                if all(x <= y for x, y in adj_pairs([str(i) for i in tree])):
                    canonical = tuple(itertools.chain.from_iterable(tree))
                    result.append(canonical)
    
    return result


#===================================================================================================================
#   Counting 112 vectors at level five in \tilde{O}(n^{5/3}) time and \tilde{O}(n) space
#===================================================================================================================


def unique_elements(lst):
    unique_list = []
    
    for item in lst:
        if item not in unique_list:
            unique_list.append(item)
    
    return unique_list


def matrix_CT(n):
    trees_n = canonical_corner_trees(n)
    list_of_vectors = []
    for y in trees_n:
        list_of_vectors.append(expand_in_patterns(y, n)[n-1])
    list_of_vectors = unique_elements(list_of_vectors)
    return(list_of_vectors)


def addmax(vector_pattern,n):

    one_line = []

    for i, tau in enumerate(permutations(range(n))):
        one_line.append([tau+(len(tau),),vector_pattern[i]])
    
    vect = math.factorial(n+1)*[0]

    for x in one_line:
        for j, sigma in enumerate(permutations(range(n+1))):
            if sigma == x[0]:
                vect[j] = x[1]
    return(vect)

def rev1(p):
    return p[::-1]

def rev2(p):
    n = len(p) - 1
    return [n - x for x in p]

def rev_vect(vector_pattern,n):

    one_line = []

    for i, tau in enumerate(permutations(range(n))):
        one_line.append([tuple(rev1(tau)),vector_pattern[i]])
    
    vect = math.factorial(n)*[0]

    for x in one_line:
        for j, sigma in enumerate(permutations(range(n))):
            if sigma == x[0]:
                vect[j] = x[1]
    return(vect)

def comp_vect(vector_pattern,n):

    one_line = []

    for i, tau in enumerate(permutations(range(n))):
        one_line.append([tuple(rev2(tau)),vector_pattern[i]])
    
    vect = math.factorial(n)*[0]
    
    for x in one_line:
        for j, sigma in enumerate(permutations(range(n))):
            if sigma == x[0]:
                vect[j] = x[1]
    return(vect)

def inv_vect(vector_pattern,n):

    one_line = []

    for i, tau in enumerate(permutations(range(n))):
        one_line.append([tuple(invperm(tau)),vector_pattern[i]])
    
    vect = math.factorial(n)*[0]
  
    for x in one_line:
        for j, sigma in enumerate(permutations(range(n))):
            if sigma == x[0]:
                vect[j] = x[1]
    return(vect)

def identity_vect(vector_pattern,n):
    return vector_pattern


def compose(f, g):
    return lambda vector_pattern, n: f(g(vector_pattern, n), n)


# D4 action: all operations
# operations = {
#     'e': identity_vect,
#     'I': inv_vect,
#     'R': rev_vect,
#     'C': comp_vect,
#     'IR': compose(inv_vect, comp_vect),
#     'IC': compose(inv_vect, comp_vect),
#     'RC': compose(rev_vect, comp_vect),
#     'IRC': compose(inv_vect, compose(rev_vect, comp_vect))
# }

#D4 actions: identity plus the ones that kick us out of Tree_{5/3}
operations = {
    'e': identity_vect,
    'R': rev_vect,
    'C': comp_vect,
    'RC': compose(rev_vect, comp_vect)
}


def new_direction_level_5():
    
    A = Matrix([[item for item in row] for row in matrix_CT(5)])
    A = A.T
    ct1 = ((SW,()), (NW, ((NW, ()),)))
    ct2 = ((NW, ((NW, ()), (SW, ()))),)
    ct3 = ((NW, ((NW, ((SW, ()),)),)),)
    ct4 = ((NW, ()), (SW, ((SW, ()),)))
    ct5 = ((NW, ((SW, ((SW, ()),)),)),)
    ct6 = ((NW, ((SW, ()),)), (SW, ()))
    ct7 = ((NW, ((SW, ()), (SW, ()))),)
    ct8 = ((NW, ()), (SW, ()), (SW, ()))
    corner_tree_five_third = [ct1, ct2, ct3, ct4, ct5, ct6, ct7, ct8]
    number_new_dir = 0
    for name, op in operations.items():
        for i in corner_tree_five_third:
            z = op(addmax(expand_in_patterns(i, 4)[3], 4),5)
            column_vector = Matrix([x for x in z])
            rank = A.rref()[1]
            new_matrix_rank =A.row_join(column_vector).rref()[1]
            if new_matrix_rank > rank:
                print(name,i)
                number_new_dir = number_new_dir+1
                A = A.row_join(column_vector)
    return number_new_dir

#===================================================================================================================
#   Level two, three, four and five profiles of a permutation
#===================================================================================================================

def profile_level_two_brute_force(permutation):
    occurrences  = {}
    for sigma in permutations(range(2)):
        occurrences[sigma] = count_brute_force(sigma,permutation)
    return(occurrences)


def profile_level_three_brute_force(permutation):
    occurrences  = {}
    for sigma in permutations(range(3)):
        occurrences[sigma] = count_brute_force(sigma,permutation)
    return(occurrences)

def profile_level_four_brute_force(permutation):
    occurrences  = {}
    for sigma in permutations(range(4)):
        occurrences[sigma] = count_brute_force(sigma,permutation)
    return(occurrences)

def profile_level_five_brute_force(permutation):
    occurrences  = {}
    for sigma in permutations(range(5)):
        occurrences[sigma] = count_brute_force(sigma,permutation)
    return(occurrences)

with open('profile_2.pkl', 'rb') as f:
    profile_2_loaded = pickle.load(f)
    all_keys_upto_2 = set()
    for lc in profile_2_loaded:
        all_keys_upto_2.update( lc.keys() )

with open('profile_3.pkl', 'rb') as f:
    profile_3_loaded = pickle.load(f)
    all_keys_upto_3 = set()
    for lc in profile_3_loaded:
        all_keys_upto_3.update( lc.keys() )

with open('profile_4.pkl', 'rb') as f:
    profile_4_loaded = pickle.load(f)
    all_keys_upto_4 = set()
    for lc in profile_4_loaded:
        all_keys_upto_4.update( lc.keys() )

with open('profile_5.pkl', 'rb') as f:
    profile_5_loaded = pickle.load(f)
    all_keys_upto_5 = set()
    for lc in profile_5_loaded:
        all_keys_upto_5.update( lc.keys() )

def profile_level_two_with_ct(perm):
    occurrences_fast = {}
    for key in all_keys_upto_2:
        occurrences_fast[key] = sum(vertex(perm,key[1]))
    occurrences  = {}
    for t,sigma in enumerate(permutations(range(2))):
            vector = profile_2_loaded[t]
            occurrences[sigma] = 0
            for key in vector.keys():
                occurrences[sigma] = occurrences[sigma] + occurrences_fast[key]*vector[key]
    return occurrences    
def profile_level_three_with_ct(perm):
    occurrences_fast = {}
    for key in all_keys_upto_3:
        occurrences_fast[key] = sum(vertex(perm,key[1]))
    occurrences  = {}
    for t,sigma in enumerate(permutations(range(3))):
            vector = profile_3_loaded[t+2]
            occurrences[sigma] = 0
            for key in vector.keys():
                occurrences[sigma] = occurrences[sigma] + occurrences_fast[key]*vector[key]
    return occurrences   
def profile_level_four_with_ct_and_3214(perm):
    occurrences_fast = {}
    for key in all_keys_upto_4:
        if key[0]=='corner tree':
            occurrences_fast[key] = sum(vertex(perm,key[1]))
        else:
            data = ast.literal_eval(key[1])
            occurrences_fast[key] = count_gen(perm,data[1],data[2],data[3],data[4],data[5])
    occurrences  = {}
    for t,sigma in enumerate(permutations(range(4))):
            vector = profile_4_loaded[t+8]
            occurrences[sigma] = 0
            for key in vector.keys():
                occurrences[sigma] = occurrences[sigma] + occurrences_fast[key]*vector[key]
    return occurrences

    
def profile_level_five_with_dp(perm, only_5=True):
    ''' This is slow: it costs O(n^{5}). Only to check that our new 12 directions work'''
    occurrences_fast = {}
    #these are elements in Tree_{5/3}
    for key in all_keys_upto_5:
        if key[0]=='tree algo':
            data = ast.literal_eval(key[1])
            if data[0] == 'e':
                occurrences_fast[key] = count_gen(perm,data[1],data[2],data[3],data[4],data[5])
            elif data[0] == 'R':
                occurrences_fast[key] = count_gen(rev1(perm),data[1],data[2],data[3],data[4],data[5])
            elif data[0] == 'C':
                occurrences_fast[key] = count_gen(rev2(perm),data[1],data[2],data[3],data[4],data[5])
            elif data[0] == 'RC':
                occurrences_fast[key] = count_gen(rev1(rev2(perm)),data[1],data[2],data[3],data[4],data[5])
        elif key[0]== 'pure pattern':
                occurrences_fast[key] = count_brute_force(key[1],perm)
        else:
        #corner trees
            occurrences_fast[key] = sum(vertex(perm,key[1]))

    occurrences  = {}
    if only_5:
        for t,sigma in enumerate(permutations(range(5))):
            vector = profile_5_loaded[t+33]
            occurrences[sigma] = 0
            for key in vector.keys():
                occurrences[sigma] = occurrences[sigma] + occurrences_fast[key]*vector[key]
    else:
        for t,sigma in enumerate(itertools.chain( *[ permutations(range(k)) for k in range(1,5+1)] )):
            vector = profile_5_loaded[t]
            occurrences[sigma] = 0
            for key in vector.keys():
                occurrences[sigma] = occurrences[sigma] + occurrences_fast[key]*vector[key]

    return occurrences


def test_profiles_2():
    n = 30
    permutation = random.sample(range(n), n)
    profile_two = profile_level_two_with_ct(permutation)
    assert profile_two==profile_level_two_brute_force(permutation)


def test_profiles_3():
    n = 30
    permutation = random.sample(range(n), n)
    profile_three = profile_level_three_with_ct(permutation)
    assert profile_three==profile_level_three_brute_force(permutation)

def test_profiles_4():
    n = 30
    permutation = random.sample(range(n), n)
    profile_four = profile_level_four_with_ct_and_3214(permutation)
    assert profile_four==profile_level_four_brute_force(permutation)

    
def test_profiles_5():
    n = 30
    permutation = random.sample(range(n), n)
    pro_dp = profile_level_five_with_dp(permutation)
    pro = profile_level_five_brute_force(permutation)
    assert pro_dp == pro

    pro_dp_all = profile_level_five_with_dp(permutation, only_5=False)
    assert pro_dp_all[(0,1)] + pro_dp_all[(1,0)] == n * (n-1) // 2
    print(pro_dp_all)



if __name__ == '__main__':
    test_counting_corner((), 1)
    test_counting_corner(((NE,()),), 2)
    test_counting_corner(((NE,()),(NW,())), 3)
    test_counting_corner(((NE,((SE,()),)),(NW,())), 4)
    test_counting_corner_pure_west((), 1)
    test_counting_corner_pure_west(((NW,((SW,()),)),(NW,())), 4)
    test_profiles_2()
    test_profiles_3()
    test_profiles_4()
    test_profiles_5()
    

    
