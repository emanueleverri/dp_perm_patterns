# Efficient Counting of Permutation Patterns via Double Posets

`dp_count_patterns.py` provides tools for efficiently counting permutation patterns using **double posets**. Precomputed profiles (`profile_2.pkl` to `profile_5.pkl`) store basis expansions for patterns of various orders.  

See the underlying paper for details: [*Efficient Counting of Permutation Patterns via Double Posets*](https://arxiv.org/abs/2408.08293).

---

## Profiles

- **`profile_2.pkl`**  
  Basis expansions of patterns of **order 2**, using corner trees with up to **2 vertices**.

- **`profile_3.pkl`**  
  Basis expansions of patterns of **order 2 and 3**, using corner trees with up to **3 vertices**.

- **`profile_4.pkl`**  
  Basis expansions of patterns of **order 2, 3, and 4**, using corner trees with up to **4 vertices**, plus the pattern `3214`.

- **`profile_5.pkl`**  
  Basis expansions of patterns of **order 2, 3, 4, and 5**, using corner trees with up to **5 vertices**, plus:  
  - Pattern `3214`  
  - Twelve double posets countable in $\tilde{O}(n^{5/3})$ – these correspond to the three double posets in $\text{Tree}_{5/3}$ that can be counted in $\tilde{O}(n^{5/3})$, together with all double posets obtained by the action of $D_4$.  
  - Additional eight level 5 patterns: `12435, 12453, 13245, 13254, 13425, 14235, 14325, 14352`  

> **Note:** These level 5 patterns are included **only for sanity checks**. The computation is not efficient (simply the brute-force algorithm) and is used to verify that all directions at level 5 can be spanned.  
> Any **eight linearly independent directions** at level 5 would suffice to cover all missing directions; these patterns are just one convenient choice.

---

## Functions

- **`new_directions()`**  
  Returns the three double posets in $\text{Tree}_{5/3}$ that can be counted in $\tilde{O}(n^{5/3})$, along with all double posets obtained by the action of $D_4$. Outputs the **underlying pure west corner trees**.

- **`profile_level_two_with_ct(perm)`**  
  Computes the **2-profile** of a permutation of length $n$ in $\tilde{O}(n)$ time.

- **`profile_level_three_with_ct(perm)`**  
  Computes the **3-profile** of a permutation of length $n$ in $\tilde{O}(n)$ time.

- **`profile_level_four_with_ct(perm)`**  
  Computes the **4-profile** of a permutation of length $n$ in $\tilde{O}(n^{5/3})$ time.

- **`profile_level_five_with_dp(perm, only_5=True)`**  
  Computes the **5-profile** of a permutation of length $n$ in $O(n^5)$ time.  
  > ⚠️ **Warning:** This function is **intentionally inefficient** and is meant only as a **sanity check** for spanning all directions at level 5.
