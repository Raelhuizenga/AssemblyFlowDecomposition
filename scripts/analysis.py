import sys
import glob
import json
import matplotlib.pyplot as plt

def main():
    # t, g_length = read_results_from_file_length("output/solutions/genomelength/*.json")
    # plot_genomelength_vs_time(t, g_length)

    t, theoretical, haps = read_results_from_file()
    plot_haps_vs_time(t, theoretical, haps)


def read_results_from_file():
    solution_files = glob.glob("output/solutions/numhaplotypes/*.json")  # Change path if needed

    times = []
    theoretical_haps = []
    num_haps = []
    num_haps_theoretical = 2

    # Read each file and extract the time
    for file in solution_files:
        with open(file, "r") as f:
            data = json.load(f)
            times.append(data["time"])
            num_haps.append(len(data["weights"]))
            theoretical_haps.append(num_haps_theoretical)
            num_haps_theoretical += 1
    print(times)
    return times, theoretical_haps, num_haps



def read_results_from_file_length(file_path):
    solution_files = glob.glob(file_path)  # Change path if needed

    times = []
    genome_lengths = []
    genome_length = 10

    # Read each file and extract the time
    for file in solution_files:
        with open(file, "r") as f:
            data = json.load(f)
            times.append(data["time"])
            genome_lengths.append(genome_length)
            genome_length += 50
    return times, genome_lengths

def plot_genomelength_vs_time(times, genome_lengths):

    plt.plot(genome_lengths, times)
    plt.xlabel("Genome length")
    plt.ylabel("Time (s)")
    plt.title("Time vs genome length")
    plt.show()
    plt.savefig('output/plot_genome_length.png')


def plot_haps_vs_time(times, theoretical_haps, num_haps):
    plt.plot(theoretical_haps, times)
    plt.xlabel("Number of haplotypes")
    plt.ylabel("Time (s)")
    plt.title("Time vs theortical number of haplotypes")
    plt.show()
    plt.savefig('output/plot_theorethical_haps.png')

    plt.scatter(num_haps, times, label="Real")
    plt.xlabel("Number of haplotypes")
    plt.ylabel("Time (s)")
    plt.title("Time vs found number of haplotypes")
    plt.show()
    plt.savefig('output/plot_found_haps.png')

if __name__ == '__main__':
    sys.exit(main())