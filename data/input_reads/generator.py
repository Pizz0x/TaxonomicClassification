import random


def main():
    random.seed(10086)
    ACGT = "ACGT"
    with open("r2_1.fa", "w") as f1, open("r2_2.fa", "w") as f2:
        for i in range(5000):
            # Generate two reads for paired-end
            seq1 = ''.join(ACGT[random.randint(0, 3)] for _ in range(300))
            seq2 = ''.join(ACGT[random.randint(0, 3)] for _ in range(300))
            
            # Write to r2_1.fa
            f1.write(f">instance{i+1}/1\n")
            f1.write(f"{seq1}\n")
            
            # Write to r2_2.fa
            f2.write(f">instance{i+1}/2\n")
            f2.write(f"{seq2}\n")


if __name__ == '__main__':
    main()