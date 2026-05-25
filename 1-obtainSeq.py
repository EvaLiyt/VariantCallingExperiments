
import pysam

## file path
input_file = ('/home/yil218/data/GCF_000001405.40_GRCh38.p14_genomic.fna')
output_file = ("/home/yil218/data/10K.fna")

## read the sequence and the name
with pysam.FastaFile(input_file) as fasta:
    name = fasta.references[0]
    complete_seq = fasta.fetch(name)

seq_length = 10000

## choose the sequence region
start_position = 44800001
end_position = start_position+ seq_length -1
print(f"start_index: {start_position}")
print(f"end_index: {end_position}")

## get the filtered sequence
filtered_seq = complete_seq[start_position:end_position+1].upper()

## write the new file
with open(output_file, "w") as outfile:
    outfile.write(">chr1\n")
    outfile.write(filtered_seq+"\n")
