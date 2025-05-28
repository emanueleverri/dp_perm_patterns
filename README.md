# Efficient counting of permutation patterns via double posets: `dp_count_patterns.py`

## Overview
This document explains the purpose of `dp_count_patterns.py`.  

- `profile_2.pkl` contains the basis expansion of patterns of order 2 using corner trees up to two vertices

- `profile_3.pkl` contains the basis expansion of patterns of order 2 and order 3 using corner trees up to three vertices

- `profile_4.pkl` contains the basis expansion of patterns of order 2, order 3, and order 4 using corner trees up to four vertices and the pattern 3214

- `profile_5.pkl` contains the basis expansion of patterns of order 2, order 3, order 4, and order 5 using corner trees up to five vertices, the pattern 3214, the twelve double posets countable in \tilde{O}(n^{5/3}) and the patterns: 12435, 12453, 13245, 13254, 13425, 14235, 14325, 14352.

- `new_directions()` yields the three double posets in Tree_{5/3} that can be counted in \tilde{O}(n^{5/3}) time together with the double posets obtained by letting D_{4} acting on them. This yields the 12 double posets that can be counted in \tilde{O}(n^{5/3}). The function actually outputs the underlying pure west corner trees.

- `profile_level_two_with_ct` computes the 2 profile of a permutation of length n in \tilde{O}(n) time

- `profile_level_three_with_ct` computes the 3 profile of a permutation of length n in \tilde{O}(n) time

- `profile_level_four_with_ct` computes the 4 profile of a permutation of length n in \tilde{O}(n^{5/3}) time

- `profile_level_five_with_dp(perm, only_5=True)` computes the 5 profile of a permutation of length n in O(n^{5}) time. It is not efficient. it is just a sanity check that we can span all directions at level five

