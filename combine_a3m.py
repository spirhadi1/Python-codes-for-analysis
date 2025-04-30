import os
import argparse
from Bio import SeqIO

def combine_unpaired_a3m(chain_a_path, chain_b_path, output_path):
    records_a = list(SeqIO.parse(chain_a_path, "fasta"))
    records_b = list(SeqIO.parse(chain_b_path, "fasta"))

    L1 = len(str(records_a[0].seq))
    L2 = len(str(records_b[0].seq))

    header_line = f"#{L1},{L2}\t1,1"
    combined_lines = [header_line]

    # Add combined reference line
    combined_lines.append(">101\t102")
    combined_lines.append(str(records_a[0].seq) + str(records_b[0].seq))

    # Chain A block
    for record in records_a:
        combined_lines.append(f">{record.id}")
        combined_lines.append(str(record.seq) + '-' * L2)

    # Chain B block
    for i, record in enumerate(records_b):
        header = ">102" if i == 0 else f">{record.id.replace(chr(9), ' ').split()[0]}"
        combined_lines.append(header)
        combined_lines.append('-' * L1 + str(record.seq))

    with open(output_path, "w") as f:
        f.write("\n".join(combined_lines) + "\n")

    print(f"✅ Combined and written to: {output_path}")


def parse_args():
    parser = argparse.ArgumentParser(description="Combine two A3M files into ColabFold unpaired multimer format.")
    parser.add_argument('--chain_a', type=str, required=True, help='Path to the first A3M file (e.g. Rad50)')
    parser.add_argument('--chain_b', type=str, required=True, help='Path to the second A3M file')
    parser.add_argument('--output', type=str, required=True, help='Path to the combined output A3M file')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    combine_unpaired_a3m(args.chain_a, args.chain_b, args.output)

