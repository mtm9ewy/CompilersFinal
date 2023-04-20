import ast
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        exit()

    base_optimizations = [
        "-mem2reg", 
        "-loop-simplify", 
        "-lcssa", 
        "-loop-rotate", 
        "-licm"
    ]

    perm = sys.argv[1]
    perm = list(ast.literal_eval(perm))

    passes = []
    with open("passes.txt", "r") as file:
        for line in file.readlines():
            passes.append(line.strip())

    optimizations = [passes[i] for i in perm]
    optimizations = base_optimizations + optimizations
    optimizations = " ".join(optimizations)

    print(optimizations)
        


