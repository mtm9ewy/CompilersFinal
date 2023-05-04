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
        if permutation_number not in scores:
            scores[permutation_number] = 0
        scores[permutation_number] += i


def get_optimizations(permutation_number):
    print(permutation_number)
    n_remaining, permutation_number = permutation_number.split("-")
    # Get the passes corresponding to each permutation number
    passes_file = f"optimization_permutations_{n_remaining}.txt"
    with open(passes_file, "r") as file:
        passes_list = [line.strip() for line in file.readlines()]

    return passes_list[int(permutation_number) - 1]


if __name__ == "__main__":
    best_file = f"results/best_optimizations_final.txt"

    gcc_results = []
    tsvc_results = []
    npb_results = []
    lcals_results = []

    for n in [16, 12, 8, 4]:
        results_file = f"results/results_{n}.csv"
        with open(results_file, "r", newline="") as file:
            reader = csv.reader(file)
            header = next(reader)

            for row in reader:
                permutation_number = f"{n}-{row[0]}"

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

    # Score[i] is the sum of the rankings for permutation i
    scores = {}
    process_sorted_results(gcc_results, scores)
    process_sorted_results(tsvc_results, scores)
    process_sorted_results(npb_results, scores)
    process_sorted_results(lcals_results, scores)

    # Sort by the sum of the rankings across the benchmarks, a lower score is better
    scores = [(score, permutation_number) for permutation_number, score in scores.items()]
    scores.sort()

    with open(best_file, "w") as file:
        # TODO: Modify this for loop to get the top n results
        for i in range(0, 1):
            score, permutation_num = scores[i]
            optimizations = get_optimizations(permutation_num)
            print(f"Score: {score}\tPermutation #: {permutation_num}\tOptimizations: {optimizations}")
            file.write(optimizations)
            file.write("\n")

