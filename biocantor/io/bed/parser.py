"""
BED file parsing.
"""

from collections import defaultdict
from pathlib import Path
from typing import Iterable

from Bio import SeqIO

from biocantor.gene import CDSInterval
from biocantor.io.exc import DuplicateSequenceException
from biocantor.io.exc import InvalidInputError
from biocantor.io.models import AnnotationCollectionModel
from biocantor.io.parser import ParsedAnnotationRecord
from biocantor.location import Strand, CompoundInterval, SingleInterval


def parse_bed(bed: Path) -> Iterable[ParsedAnnotationRecord]:
    genes = defaultdict(list)
    with open(bed) as fh:
        for row in fh:
            if row.startswith("#") or row.startswith("track "):
                continue
            row = row.rstrip().split("\t")
            if len(row) == 3:
                tx = dict(
                    exon_starts=[int(row[1])],
                    exon_ends=[int(row[2])],
                    strand=Strand.PLUS.name,
                    sequence_name=row[0],
                )
            elif len(row) == 6:
                tx = dict(
                    exon_starts=[int(row[1])],
                    exon_ends=[int(row[2])],
                    strand=Strand.from_symbol(row[5]).name,
                    sequence_name=row[0],
                    transcript_symbol=row[3],
                )
            elif len(row) == 12:
                start = int(row[1])
                strand = Strand.from_symbol(row[5])
                thick_start = int(row[6])
                thick_end = int(row[7])
                block_sizes = list(map(int, row[10].split(",")))
                block_starts = list(map(int, row[11].split(",")))
                exon_starts = [start + x for x in block_starts]
                exon_ends = [start + x + y for x, y in zip(block_starts, block_sizes)]
                exon_interval = CompoundInterval(exon_starts, exon_ends, strand)
                cds_interval = SingleInterval(thick_start, thick_end, strand)
                cds_blocks = exon_interval.intersection(cds_interval)
                cds_starts = [x.start for x in cds_blocks.blocks]
                cds_ends = [x.end for x in cds_blocks.blocks]
                frames = CDSInterval.construct_frames_from_location(cds_blocks)
                tx = dict(
                    exon_starts=exon_starts,
                    exon_ends=exon_ends,
                    strand=strand.name,
                    cds_starts=cds_starts,
                    cds_ends=cds_ends,
                    cds_frames=[x.name for x in frames],
                    transcript_symbol=row[3],
                    sequence_name=row[0],
                )
            else:
                raise InvalidInputError(f"Found BED row of length {len(row)}. Must be 3, 6 or 12.")
            genes[tx["sequence_name"]].append(dict(transcripts=[tx], sequence_name=tx["sequence_name"]))
    for chrom, gene_dicts in genes.items():
        yield ParsedAnnotationRecord(
            AnnotationCollectionModel.Schema().load(dict(genes=gene_dicts, sequence_name=chrom))
        )


def parse_bed_fasta(bed: Path, fasta: Path) -> Iterable[ParsedAnnotationRecord]:
    seqrecords_dict = {}
    for seq_rec in SeqIO.parse(fasta, format="fasta"):
        if seq_rec.id in seqrecords_dict:
            raise DuplicateSequenceException(f"Sequence {seq_rec.id} found twice in FASTA file.")
        seqrecords_dict[seq_rec.id] = seq_rec
    for rec in parse_bed(bed):
        rec.seqrecord = seqrecords_dict[rec.annotation.sequence_name]
        yield rec
