def calculate_n50(fasta_path):
    contig_lengths = []
    with open(fasta_path) as f:
        length = 0
        for line in f:
            if line.startswith('>'):
                if length > 0:
                    contig_lengths.append(length)
                    length = 0
            else:
                length += len(line.strip())
        if length > 0:
            contig_lengths.append(length)

    contig_lengths.sort(reverse=True)
    total_length = sum(contig_lengths)
    half_length = total_length / 2

    running_sum = 0
    for length in contig_lengths:
        running_sum += length
        if running_sum >= half_length:
            return length


n50 = calculate_n50('../data/data1/contig.fasta')
print("N50 for data1:", n50)
n50 = calculate_n50('../data/data2/contig.fasta')
print("N50 for data2:", n50)
n50 = calculate_n50('../data/data3/contig.fasta')
print("N50 for data3:", n50)
# n50 = calculate_n50('data4/contig.fasta')
# print("N50 for data4:", n50)
