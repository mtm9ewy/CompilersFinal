import csv
from datetime import datetime, timedelta


def append_to_results(header_name, results, entry):
    if "size" in header_name:
        results[1] = int(entry)
    elif "compile" in header[i]:
        result = datetime.strptime(entry, "%Mm%S.%fs")
        results[2] = timedelta(
            minutes=result.minute, 
            seconds=result.second, 
            microseconds=result.microsecond
        ).total_seconds()
    else:
        results[0] += float(entry)


def process_sorted_results(results, scores):
    for i, result in enumerate(results):
        permutation_number = result[3]
        scores[permutation_number - 1][0] += i


if __name__ == "__main__":
    results_file = "results/results_16.csv"
    passes_file = "optimization_permutations_16.txt"
    best_file = "results/best_optimizations_16.txt"

    gcc_results = []
    tsvc_results = []
    npb_results = []
    lcals_results = []

    with open(results_file, "r", newline="") as file:
        reader = csv.reader(file)
        header = next(reader)

        for row in reader:
            permutation_number = int(row[0])

            # Each result is sum of runtime, executable size, compile time, permutation #
            gcc_result = [0, 0, 0, permutation_number]
            tsvc_result = [0, 0, 0, permutation_number]
            npb_result = [0, 0, 0, permutation_number]
            lcals_result = [0, 0, 0, permutation_number]

            for i, entry in enumerate(row):
                header_name = header[i]
                if "gcc" in header_name:
                    append_to_results(header_name, gcc_result, entry)
                elif "tsvc" in header_name:
                    append_to_results(header_name, tsvc_result, entry)
                elif "npb" in header_name:
                    append_to_results(header_name, npb_result, entry)
                elif "lcals" in header_name:
                    append_to_results(header_name, lcals_result, entry)

            gcc_results.append(gcc_result)
            tsvc_results.append(tsvc_result)
            npb_results.append(npb_result)
            lcals_results.append(lcals_result)

        # Sort the results for each benchmark, the sort is by
        # 1) Benchmark runtime, 2) Benchmark executable size, 3) Benchmark compile time
        gcc_results.sort()
        tsvc_results.sort()
        npb_results.sort()
        lcals_results.sort()

        # Score[i] is the sum of the rankings for permutation i + 1 and permutation i + 1
        scores = [[0, i + 1] for i in range(0, len(gcc_results))]
        process_sorted_results(gcc_results, scores)
        process_sorted_results(tsvc_results, scores)
        process_sorted_results(npb_results, scores)
        process_sorted_results(lcals_results, scores)

        # Sort by the sum of the rankings across the benchmarks, a lower score is better
        scores.sort()

        # Get the passes corresponding to each permutation number
        with open(passes_file, "r") as file:
            passes_list = [line.strip() for line in file.readlines()]

        with open(best_file, "w") as file:
            for i in range(0, 5):
                score, permutation_num = scores[i]
                print(f"Score: {score}\tPermutation #: {permutation_num}\tOptimizations: {passes_list[permutation_num - 1]}")
                file.write(passes_list[permutation_num - 1])
                file.write("\n")

