import csv
import sys

def parse_gcc_loops_results(filepath):
    results = []
    with open(filepath, "r", newline="") as file:
        reader = csv.reader(file)
        # Consume the header
        next(reader)
        print(reader)
        for line in reader:
            # Loop Name, Time, Units
            results.append((line[0], line[1]))

    return results


def parse_npb_results(filepath):
    results = []
    with open(filepath, "r") as file:
        for line in file.readlines():
            if "Benchmark Completed" in line:
                test_name = line.split()[0]

            if "Time in seconds = " in line:
                results.append((test_name, line.split("=")[1].strip()))

    return results


def parse_lcals_results(directory_path):
    results = []

    for loop_type in ["Forall_Lambda", "Raw"]:
        filepath = f"{directory_path}/{loop_type}-meantime.txt"
        with open(filepath, "r", newline="") as file:
            reader = csv.reader(file)
            # Read the first header
            next(reader)

            # Read the second header
            loop_lengths = next(reader)
            loop_lengths = loop_lengths[1:]

            for line in reader:
                for i, loop_length in enumerate(loop_lengths):
                    results.append((f"{loop_type}_{line[0].strip()}_{loop_length.strip()}", line[i + 1]))

    return results


if __name__ == "__main__":
    if len(sys.argv) != 4:
        exit()

    perm, gcc_loops_results_path, lcals_results_dir_path, npb_results_path = sys.argv
    
    results = parse_gcc_loops_results(gcc_loops_results_path)
    results = parse_lcals_results(lcals_results_dir_path)
    results = parse_npb_results(npb_results_path)
