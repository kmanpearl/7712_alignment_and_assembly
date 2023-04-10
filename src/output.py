import pandas as pd


def get_longest_contig(contigs):
    length_longest_contig = 0
    for contig in contigs:
        if len(contig.sequence) > length_longest_contig:
            length_longest_contig = len(contig.sequence)
            longest_contig_id = contig.contig_id
            longest_contig_sequence = contig.sequence
            longest_contig = contig
    with open("output/ALLELES.fasta", "w") as file:
        file.write(f">contig_{longest_contig_id}\n{longest_contig_sequence}")
    return longest_contig


def format_aln_file(longest_contig):
    reads = longest_contig.aligned_reads
    sseqid = []
    qseqid = []
    sstart = []
    send = []
    qstart = []
    qend = []
    for read in reads:
        sseqid.append(read.read_id)
        qseqid.append(f"contig{str(read.contig_id)}")
        sstart.append(read.read_start)
        send.append(read.read_stop)
        qstart.append(read.contig_start)
        qend.append(read.contig_stop)
    aln = pd.DataFrame(
        {
            "sseqid": sseqid,
            "cseqid": qseqid,
            "sstart": sstart,
            "send": send,
            "qstart": qstart,
            "qend": qend,
        }
    )
    aln.to_csv("output/ALLELES.aln")


def save_required_ouputs(contigs):
    longest_contig = get_longest_contig(contigs)
    format_aln_file(longest_contig)
