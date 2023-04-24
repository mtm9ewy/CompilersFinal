import csv
import sys

RESULTS_DIR="results"


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
    if len(sys.argv) != 2:
        print("Usage: python3 aggregate_results.py <Number of result files>")

    num_files = int(sys.argv[1])

    results = []
    for i in range(1, num_files + 1):
        data = {"permutation_number": i}

        for benchmark in ["gcc", "tsvc", "npb", "lcals"]:
            process_benchmark(benchmark, data)

        results.append(data)


    with open("results.csv", "w", newline="") as file:
        writer = csv.DictWriter(file, results[0].keys())

        writer.writeheader()
        for result in results:
            writer.writerow(result)

