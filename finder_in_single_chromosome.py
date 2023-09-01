import twobitreader, re, time, math, sys, pandas, Bio
import subprocess as sp
from tqdm import tqdm
from Bio import Seq


chr_name = "chrX"

def remove_keys_except(diction, keys_to_keep):
    keys_to_remove = [key for key in diction if key not in keys_to_keep]
    for key in keys_to_remove:
        diction.pop(key, None)

f1 = open("/home/shan/cas9cpf1/atac_open_finder/whole_open_region","r")
f2 = open("/home/shan/cas9cpf1/atac_open_finder/open_region/"+chr_name+"_open_region","r")
sequence = f1.readline()
chr_sequence = f2.readline()
f1.close()
f2.close()

path = "/home/shan/cas9cpf1/ref/genome.2bit"
WGS = twobitreader.TwoBitFile(path)

# [chr1-X]
keep_chr = ['chr' + str(i) for i in range(1, 23)]
keep_chr.append('chrX')
remove_keys_except(WGS, keep_chr)

length = len(WGS[chr_name])
length_open = len(chr_sequence)

whole_sequence = WGS[chr_name][0:length]
whole_sequence = whole_sequence.upper()
whole_sequence_rev = Seq.reverse_complement(whole_sequence)

print("Makine whole sequence is done!")

candidate = []

def find_unique(args):
    index_start, index_end = args
    local_candidate = []

    i = index_start

    while i < index_end - 27:

        test_raw = chr_sequence[i:i+27]
        test = test_raw[4:24]

        if test_raw[26:27] == "N":
            i += 53

        else:
            first_index = whole_sequence.find(test)
            second_index = whole_sequence.find(test,first_index+1)

            test_prefix = test_raw[0:3]

            if test_prefix == "TTT":
                if test_raw[3] != "T":
                    if second_index == -1:
                        rev_first_index = whole_sequence_rev.find(test)
                        if rev_first_index == -1:
                            test_suffix = test_raw[25:27]
                            if test_suffix == "GG":
                                print(test_raw,first_index,i)
                                local_candidate.append([test_raw,first_index])

                    skip_search = test_raw[3:27]
                    skip_index = skip_search.find("TTT")
                    if skip_index != -1:
                        i += (skip_index + 3)

                    else:
                        i += 25

                else:
                    i += 1

            else:
                skip_index = test_raw.find("TTT")
                if skip_index != -1:
                    i += skip_index

                else:
                    i += 1

    return local_candidate


if __name__ == '__main__':
    from multiprocessing import Pool

    print("running")

    start_time = time.time()

    num_cores = 10
    chunk_size = (length_open + num_cores - 1) // num_cores

    ranges = [(i, min(i + chunk_size, length_open)) for i in range(0, length_open, chunk_size)]
    print(ranges)


    with Pool(processes=num_cores) as pool:
        results = list(pool.imap(find_unique, ranges))


    print("Running is done! Start to assemble the candidate sequence!")

    candidate = [result for sublist in results for result in sublist]

    print(candidate)

    fw1 = open("/home/shan/cas9cpf1/atac_open_finder/candidate_seq/"+chr_name+"_candidate_sequence","w")
    fw1.write("[")
    for sequence in candidate:
        write_string = '["'
        write_string += sequence[0]
        write_string += '",'
        write_string += str(sequence[1])
        write_string += "]"
        write_string += ","
        fw1.write(write_string)
    fw1.write("]")
    fw1.close()


    end_time = time.time()
    execution_time = end_time - start_time

    print(f"Run time: {execution_time:.5f}s")

