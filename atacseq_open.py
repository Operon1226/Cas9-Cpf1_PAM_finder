import twobitreader, re, time, math, sys, pandas
import subprocess as sp
from tqdm import tqdm
import Bio

def atac_open_list():
    
    # bed_file = input("Please input Peak called ATAC-seq Data\n")

    bed_file = "/home/intern/GEL_ATAC_cellLines/HCT116-1.processed.files/Peak_calling/HCT116-1_peaks.narrowPeak"

    df = pandas.read_csv(bed_file,sep='\t',header=None)

    order_of_chr = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15',
                    'chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
    
    df[0] = pandas.Categorical(df[0],categories=order_of_chr,ordered=True)
    
    sorted_df = df.sort_values(by=[0,1],ascending=[True,True])

    # [chr1-X]
    keep_chr = ['chr' + str(i) for i in range(1, 23)]
    keep_chr.append('chrX')

    index = 0
    row_index = 0
    row_index_list = []

    while True:
        while True:
            if keep_chr[index] == sorted_df.iloc[row_index,0]:
                row_index += 1
            else:   
                break
        row_index_list.append(row_index-1)
        if index == 22:
            break
        index += 1
        
    open_list = dict()
    row_index = 0
    
    for i in range(len(row_index_list)):
        if i == 0:
            open_list[keep_chr[i]] = [0,row_index_list[i]]
            row_index = row_index_list[i] + 1
        
        elif row_index_list[i] != row_index_list[i-1]:
            open_list[keep_chr[i]] = [row_index,row_index_list[i]]
            row_index = row_index_list[i] + 1
    
    print(open_list)
    
    return (sorted_df, open_list, keep_chr)

def remove_keys_except(diction, keys_to_keep):
    keys_to_remove = [key for key in diction if key not in keys_to_keep]
    for key in keys_to_remove:
        diction.pop(key, None)

def search_tool(chr_name,sorted_df,open_list,keep_chr):
    
    # path = input("Please input path of genome.2bit\n")
    
    path = "/home/shan/cas9cpf1/ref/genome.2bit"
    
    WGS = twobitreader.TwoBitFile(path)

    remove_keys_except(WGS, keep_chr)
    
    # sequence of open region
    sequence_open = ""
    for i in range(open_list[chr_name][0],open_list[chr_name][1]+1):
        sequence_open += WGS[chr_name][sorted_df.iloc[i,1]:sorted_df.iloc[i,2]]
        sequence_open += "NNNNNNNNNNNNNNNNNNNNNNNNNNN"
    
    sequence_open = sequence_open.upper()
    
    # # make one fasta file for all open region
    # f = open("total.open.fasta","a")
    # f.write(">"+chr_name+"\n")
    # for i in range(math.ceil(len(sequence_open)/50)):
    #     f.write(sequence_open[i*50:(i+1)*50]+"\n")
    # f.close()
    
    # # make one fasta file for each chromosome
    # f = open(chr_name+".open.fasta","w")
    # f.write(">"+chr_name+"\n")
    # for i in range(math.ceil(len(sequence_open)/50)):
    #     f.write(sequence_open[i*50:(i+1)*50]+"\n")
    # f.close()

    f = open("./open_region/"+chr_name+"_open_region","w")
    f.write(sequence_open)
    f.close()

    # # make only string file
    # f = open("whole_open_region","a")
    # f.write(sequence_open)
    # f.write("EEEEEEEEEE")
    # f.close()
    
    print(f"Length of {chr_name}: {len(sequence_open)}")
    return len(sequence_open)
    
if __name__ == "__main__":
    
    start_time = time.time()

    (sorted_df,open_list,keep_chr) = atac_open_list()
    
    # b = 0
    # for chromosome in keep_chr:
    #     a = search_tool(chromosome,sorted_df,open_list,keep_chr)
    #     b += a
    
    for chromosome in keep_chr:
        if chromosome != 'chr22':
            a = search_tool(chromosome,sorted_df,open_list,keep_chr)
    
        
    # print(f"Total length of open region: {b}")
    end_time = time.time()
    execution_time = end_time - start_time

    print(f"Run time: {execution_time:.5f}s")
