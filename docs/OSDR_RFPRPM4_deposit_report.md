# Draft report to the NASA OSDR curators: GLDS-528 RFPRPM4 raw FASTQ files

**Status: draft, not sent.** Prepared on 17 October 2026 for Richard Barker to review and send.
Nothing in it has been sent to OSDR or to the study's depositors.

---

**Subject:** OSD-528 / GLDS-528: RFPRPM4 raw FASTQs contain the RFPRPM41 files followed by a second library

Dear OSDR curation team,

While re-analysing OSD-528 (*Mycobacterium marinum* under a 3D clinostat, an RPM 2.0 and static 1g;
DOI [10.26030/r3re-fd65](https://doi.org/10.26030/r3re-fd65)), I found that both raw FASTQ files for
sample RFPRPM4 begin with the complete raw FASTQ files of sample RFPRPM41. A second, different library
follows. The details below should let you check this quickly.

## What the files contain

The files were downloaded from the OSDR file API and each matched the deposited
`GLDS-528_rna-seq_raw_md5sum.txt`.

| File | Size (bytes) | md5 |
|---|---:|---|
| `GLDS-528_rna-seq_RFPRPM4_R1_raw.fastq.gz` | 1,553,788,471 | `bc83634f12641531804942debc4c445f` |
| `GLDS-528_rna-seq_RFPRPM4_R2_raw.fastq.gz` | 831,240,941 | `8c7d36e68d8ceaf81486cd43397c5c55` |
| `GLDS-528_rna-seq_RFPRPM41_R1_raw.fastq.gz` | 748,733,492 | `5f477dacb81f09b5932047d377b5518e` |
| `GLDS-528_rna-seq_RFPRPM41_R2_raw.fastq.gz` | 391,318,381 | `dfa32f5ed3bde08c62d02d986544c2e3` |

1. **The first part of each RFPRPM4 file is the whole RFPRPM41 file.** The md5 of the first
   748,733,492 bytes of RFPRPM4 R1 is `5f477dacb81f09b5932047d377b5518e`, which is the deposited md5 of
   RFPRPM41 R1. Likewise, the first 391,318,381 bytes of RFPRPM4 R2 have the RFPRPM41 R2 md5,
   `dfa32f5ed3bde08c62d02d986544c2e3`.
2. **A new gzip member starts exactly where the RFPRPM41 content ends:** the bytes `1f 8b` at offset
   748,733,492 in R1 and 391,318,381 in R2. The remainder is 805,054,979 bytes in R1 (md5
   `d71dd4afc4f80bbd86e5d69375ea3fd5`) and 439,922,560 bytes in R2 (md5 `9eec020b7322f6e1be1b04e8ddef42cd`).
3. **The remainder is a different library from the same run.** Its reads come from the same NextSeq run
   and flow cell (`NB551111:147:HT22VAFX2`) but carry index `CACAGCGGTC+NAATAGGAAT`. RFPRPM41's reads
   carry `AGCCTATGAT+NTAGGTGGTG`. The remainder index is not used by any other GLDS-528 sample. It holds
   15,304,202 read pairs; RFPRPM41 holds 14,370,842.

You can reproduce points 1 and 2 with standard tools, for example:

```bash
head -c 748733492 GLDS-528_rna-seq_RFPRPM4_R1_raw.fastq.gz | md5     # = RFPRPM41 R1 md5
tail -c +748733493 GLDS-528_rna-seq_RFPRPM4_R1_raw.fastq.gz | gunzip -c | head -1
```

## Why it matters

Anyone who reads the RFPRPM4 files from the start gets RFPRPM41's reads, so RFPRPM4 duplicates RFPRPM41
in every count matrix built that way. The RPM group then has two independent samples, not three. An
earlier version of my own re-analysis made exactly this mistake.

## Request

Could you check with the depositors whether the second library is the true RFPRPM4 library? If it is,
please replace the RFPRPM4 R1 and R2 files with only that library and update the md5 list.

In my re-analysis I treat the second library as an *inferred* RFPRPM4. It's excluded from the primary
analysis and used only in a sensitivity analysis until this is confirmed. The scripts that detect and
split the files are in
[`analysis/00_quantify_real_osdr_rnaseq.py`](../analysis/00_quantify_real_osdr_rnaseq.py)
(function `split_rpm4`) in https://github.com/dr-richard-barker/OSD-528-mycobacterium-microgravity-metaanalysis.

Kind regards,
Richard Barker
Collaborative Science Environment Inc., Madison, Wisconsin, USA
