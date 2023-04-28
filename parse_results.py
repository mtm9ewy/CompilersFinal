import csv
import sys

def parse_gcc_loops_results(filepath):
    results = []
    with open(filepath, "r", newline="") as file:
        reader = csv.reader(file)
        # Consume the header
        next(reader)
        for line in reader:
            # Loop Name, Time, Units
            results.append((line[0].strip(), line[1].strip()))

    return results


def parse_npb_results(filepath):
    results = []
    with open(filepath, "r") as file:
        for line in file.readlines():
            if "Benchmark Completed" in line:
                test_name = line.split()[0].strip()

            if "Time in seconds = " in line:
                results.append((test_name, line.split("=")[1].strip()))

    return results


def parse_lcals_results(directory_path):
    results = []

    loop_type = "Raw"
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
                results.append((f"{loop_type}_{line[0].strip()}_{loop_length.strip()}", line[i + 1].strip()))

    return results


def parse_tsvc_results(filepath):
    results = []
    with open(filepath, "r") as file:
        for i, line in enumerate(file.readlines()):
            # Consume headers
            if i < 2:
                continue
            loop_name, time, _ = line.split()
            results.append((loop_name.strip(), time.strip()))

    return results


def write_results(filename, results):
    with open(f"temp_results/{filename}", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([result[0] for result in results])
        writer.writerow([result[1] for result in results])


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(sys.argv)
        exit()


    _, perm, gcc_loops_results_path, lcals_results_dir_path, npb_results_path, tsvc_results_path = sys.argv
    
    results = parse_gcc_loops_results(gcc_loops_results_path)
    write_results(f"{perm}_gcc.csv", results)

    results = parse_lcals_results(lcals_results_dir_path)
    write_results(f"{perm}_lcals.csv", results)

    results = parse_npb_results(npb_results_path)
    write_results(f"{perm}_npb.csv", results)

    results = parse_tsvc_results(tsvc_results_path)
    write_results(f"{perm}_tsvc.csv", results)
