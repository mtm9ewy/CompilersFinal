import sys
from itertools import permutations

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 create_optimizations_list <path to file containing initial optimizations>")
        print("If there are no initial optimizations, pass an empty file")
        exit()

    # All 16 optimizations
    with open("passes.txt", "r") as file:
        all_optimizations = [line.strip() for line in file.readlines()]

    # Optimizations that we always run
    base_optimizations = [
        "-mem2reg", 
        "-loop-simplify", 
        "-lcssa", 
        "-loop-rotate", 
        "-licm"
    ]

    # Optimizations that we found to be the "best" that we want to use
    # for future searches
    with open(sys.argv[1], "r") as file:
        initial_optimizations = [line.strip().split()[len(base_optimizations):] for line in file.readlines()]

    
    # The optimizations that are left to chose from
    if len(initial_optimizations) != 0:
        remaining_optimizations = [] 
        for i, optimization_list in enumerate(initial_optimizations):
            remaining_optimizations.append(list(set(all_optimizations) - set(optimization_list)))
    else:
        remaining_optimizations = [all_optimizations]


    n_remaining = len(remaining_optimizations[0])
    with open(f"optimization_permutations_{n_remaining}.txt", "w") as file:
        for perm in permutations(range(0, n_remaining), 4):
            if len(initial_optimizations) == 0:
                file.write(" ".join(base_optimizations + [remaining_optimizations[0][j] for j in perm]))
                file.write("\n")
                continue

            for i, optimizations in enumerate(initial_optimizations):
                file.write(" ".join(base_optimizations + optimizations + [remaining_optimizations[i][j] for j in perm]))
                file.write("\n")

