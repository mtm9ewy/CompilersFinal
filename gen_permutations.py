from itertools import permutations
total_options = int(input("Type the number of options you are choosing from: "))
set_size = int(input("Type the number of values per permutation you want: "))

p = list(permutations(range(total_options), set_size))
file = open("permutations.txt", 'w')
for perm in p:
    file.write(str(perm) + "\n")
file.close()


