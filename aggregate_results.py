import csv
import sys

RESULTS_DIR="temp_results"


def process_benchmark(name, data):
    with open(f"{RESULTS_DIR}/{i}_{name}.csv", "r", newline="") as file:
        reader = csv.reader(file)
        header = next(reader)
        timings = next(reader)

        for j in range(0, len(header)):
            data[f"{name}_{header[j]}"] = timings[j]

    with open(f"{RESULTS_DIR}/{i}_{name}_size.txt", "r") as file:
        _ = file.readline()
        size, _ = file.readline().split()
        data[f"{name}_size"] = size

    with open(f"{RESULTS_DIR}/{i}_{name}_compile.txt", "r") as file:
        for line in file.readlines():
            if "real" in line:
                _, timing = line.split()
                break
        data[f"{name}_compile"] = timing


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 aggregate_results.py <Permutation # of 1st result> <Permutation # of last result>")
        exit()

    start = int(sys.argv[1])
    stop = int(sys.argv[2])

    results = []
    for i in range(start, stop + 1):
        data = {"permutation_number": i}

        for benchmark in ["gcc", "tsvc", "npb", "lcals"]:
            process_benchmark(benchmark, data)

        results.append(data)


    with open("results.csv", "w", newline="") as file:
        writer = csv.DictWriter(file, results[0].keys())

        writer.writeheader()
        for result in results:
            writer.writerow(result)

