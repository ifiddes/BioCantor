"""
FeatureIntervals have the ability to write to BED.

"""

import pytest
from collections import OrderedDict
from biocantor.gene.cds_frame import CDSFrame
from biocantor.io.bed import RGB
from biocantor.io.bed.parser import parse_bed, parse_bed_fasta
from biocantor.io.models import (
    TranscriptIntervalModel,
    FeatureIntervalModel,
    AnnotationCollectionModel,
)
from biocantor.location.strand import Strand
from biocantor.util.hashing import digest_object


class TestBedWriter:
    tx1 = dict(
        exon_starts=[2],
        exon_ends=[18],
        strand=Strand.PLUS.name,
        cds_starts=[5],
        cds_ends=[9],
        cds_frames=[CDSFrame.ZERO.name],
        sequence_name="chr1",
    )
    tx2 = dict(
        exon_starts=[2, 7, 12],
        exon_ends=[6, 10, 15],
        strand=Strand.PLUS.name,
        cds_starts=[4, 7, 12],
        cds_ends=[6, 10, 13],
        cds_frames=[CDSFrame.ZERO.name, CDSFrame.TWO.name, CDSFrame.TWO.name],
        transcript_interval_guid=digest_object(123),
        transcript_symbol="name",
    )
    feat1 = dict(interval_starts=[2], interval_ends=[5], strand=Strand.PLUS.name, sequence_name="chr10")
    feat2 = dict(interval_starts=[2, 7, 12], interval_ends=[6, 10, 15], strand=Strand.PLUS.name)
    feat3 = dict(interval_starts=[25], interval_ends=[30], strand=Strand.MINUS.name)

    @pytest.mark.parametrize(
        "tx,expected",
        [
            (tx1, ["chr1", "2", "18", "None", "0", "+", "5", "9", "0,0,0", "1", "16", "0"]),
            (
                tx2,
                [
                    "None",
                    "2",
                    "15",
                    "name",
                    "0",
                    "+",
                    "4",
                    "13",
                    "0,0,0",
                    "3",
                    "4,3,3",
                    "0,5,10",
                ],
            ),
        ],
    )
    def test_tx(self, tx, expected):
        model = TranscriptIntervalModel.Schema().load(tx)
        obj = model.to_transcript_interval()
        assert str(obj.to_bed12()) == "\t".join(expected)

    @pytest.mark.parametrize(
        "feat,expected",
        [
            (feat1, ["chr10", "2", "5", "None", "0", "+", "0", "0", "0,0,0", "1", "3", "0"]),
            (feat2, ["None", "2", "15", "None", "0", "+", "0", "0", "0,0,0", "3", "4,3,3", "0,5,10"]),
            (feat3, ["None", "25", "30", "None", "0", "-", "0", "0", "0,0,0", "1", "5", "0"]),
        ],
    )
    def test_feat(self, feat, expected):
        model = FeatureIntervalModel.Schema().load(feat)
        obj = model.to_feature_interval()
        assert str(obj.to_bed12()) == "\t".join(expected)

    @pytest.mark.parametrize(
        "tx,score,rgb,name,expected",
        [
            (
                tx2,
                10,
                RGB(128, 128, 128),
                "transcript_symbol",
                ["None", "2", "15", "name", "10", "+", "4", "13", "128,128,128", "3", "4,3,3", "0,5,10"],
            ),
            (  # if name is not an attribute, just pass it along
                tx2,
                10,
                RGB(128, 128, 128),
                "test",
                ["None", "2", "15", "test", "10", "+", "4", "13", "128,128,128", "3", "4,3,3", "0,5,10"],
            ),
        ],
    )
    def test_changed_metadata(self, tx, score, rgb, name, expected):
        model = TranscriptIntervalModel.Schema().load(tx)
        obj = model.to_transcript_interval()
        assert str(obj.to_bed12(score, rgb, name)) == "\t".join(expected)

    def test_parse_bed(self, tmp_path):
        tmp_bed = tmp_path / "tmp.bed"
        with open(tmp_bed, "w") as fh:
            for tx in [self.tx1, self.tx2]:
                model = TranscriptIntervalModel.Schema().load(tx)
                obj = model.to_transcript_interval()
                print(obj.to_bed12(), file=fh)
        parsed = list(parse_bed(tmp_bed))
        assert AnnotationCollectionModel.Schema().dump([x.annotation for x in parsed], many=True) == [
            OrderedDict(
                [
                    ("feature_collections", []),
                    (
                        "genes",
                        [
                            OrderedDict(
                                [
                                    (
                                        "transcripts",
                                        [
                                            OrderedDict(
                                                [
                                                    ("exon_starts", [2]),
                                                    ("exon_ends", [18]),
                                                    ("strand", "PLUS"),
                                                    ("cds_starts", [5]),
                                                    ("cds_ends", [9]),
                                                    ("cds_frames", ["ZERO"]),
                                                    ("qualifiers", None),
                                                    ("is_primary_tx", None),
                                                    ("transcript_id", None),
                                                    ("protein_id", None),
                                                    ("product", None),
                                                    ("transcript_symbol", "None"),
                                                    ("transcript_type", None),
                                                    ("sequence_name", "chr1"),
                                                    ("sequence_guid", None),
                                                    ("transcript_interval_guid", None),
                                                    ("transcript_guid", None),
                                                ]
                                            )
                                        ],
                                    ),
                                    ("gene_id", None),
                                    ("gene_symbol", None),
                                    ("gene_type", None),
                                    ("locus_tag", None),
                                    ("qualifiers", None),
                                    ("sequence_name", "chr1"),
                                    ("sequence_guid", None),
                                    ("gene_guid", None),
                                ]
                            )
                        ],
                    ),
                    ("variant_collections", []),
                    ("name", None),
                    ("id", None),
                    ("sequence_name", "chr1"),
                    ("sequence_guid", None),
                    ("sequence_path", None),
                    ("qualifiers", None),
                    ("start", None),
                    ("end", None),
                    ("completely_within", None),
                    ("parent_or_seq_chunk_parent", None),
                ]
            ),
            OrderedDict(
                [
                    ("feature_collections", []),
                    (
                        "genes",
                        [
                            OrderedDict(
                                [
                                    (
                                        "transcripts",
                                        [
                                            OrderedDict(
                                                [
                                                    ("exon_starts", [2, 7, 12]),
                                                    ("exon_ends", [6, 10, 15]),
                                                    ("strand", "PLUS"),
                                                    ("cds_starts", [4, 7, 12]),
                                                    ("cds_ends", [6, 10, 13]),
                                                    ("cds_frames", ["ZERO"]),
                                                    ("qualifiers", None),
                                                    ("is_primary_tx", None),
                                                    ("transcript_id", None),
                                                    ("protein_id", None),
                                                    ("product", None),
                                                    ("transcript_symbol", "name"),
                                                    ("transcript_type", None),
                                                    ("sequence_name", "None"),
                                                    ("sequence_guid", None),
                                                    ("transcript_interval_guid", None),
                                                    ("transcript_guid", None),
                                                ]
                                            )
                                        ],
                                    ),
                                    ("gene_id", None),
                                    ("gene_symbol", None),
                                    ("gene_type", None),
                                    ("locus_tag", None),
                                    ("qualifiers", None),
                                    ("sequence_name", "None"),
                                    ("sequence_guid", None),
                                    ("gene_guid", None),
                                ]
                            )
                        ],
                    ),
                    ("variant_collections", []),
                    ("name", None),
                    ("id", None),
                    ("sequence_name", "None"),
                    ("sequence_guid", None),
                    ("sequence_path", None),
                    ("qualifiers", None),
                    ("start", None),
                    ("end", None),
                    ("completely_within", None),
                    ("parent_or_seq_chunk_parent", None),
                ]
            ),
        ]

    def test_spliced(self, tmp_path):
        tmp_bed = tmp_path / "tmp.bed"
        with open(tmp_bed, "w") as fh:
            fh.write("chr1\t926007\t930198\tENSG00000187634|SAMD11|ce0f9f1\t0\t-\t926007\t930198\t0\t2\t6,44\t0,4147\n")
        recs = parse_bed(tmp_bed)
        _ = list(recs)[0].to_annotation_collection()

    def test_parse_bed_fasta(self, tmp_path):
        tmp_bed = tmp_path / "tmp.bed"
        tmp_fasta = tmp_path / "tmp.fasta"
        with open(tmp_bed, "w") as fh:
            model = TranscriptIntervalModel.Schema().load(self.tx1)
            obj = model.to_transcript_interval()
            print(obj.to_bed12(), file=fh)
        with open(tmp_fasta, "w") as fh:
            fh.write(">chr1\nATGCATGCATGAGAAGACCCGAGTAA\n")
        parsed = list(parse_bed_fasta(tmp_bed, tmp_fasta))
        assert len(parsed) == 1
        rec = parsed[0].to_annotation_collection()
        assert str(rec.genes[0].get_primary_transcript().get_spliced_sequence()) == "GCATGCATGAGAAGAC"
