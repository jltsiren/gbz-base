# The Graph Alignment Format (GAF)

This document describes the [vg](https://github.com/vgteam/vg) version of the Graph Alignment Format (GAF).
It is a superset of a subset of the [original GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md).
That format in turn is a superset of the [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md).
Sequence names and optional fields follow conventions set in the [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).
Difference strings are defined in the [minimap2 man page](https://lh3.github.io/minimap2/minimap2.html).
Paths are represented as [GFA](https://gfa-spec.github.io/GFA-spec/GFA1.html) walks.

## Overview

GAF is a tab-delimited file format for sequence alignments to bidirected sequence graphs.
The file is encoded in UTF-8.
Unless otherwise specified, all fields are restricted to 7-bit US-ASCII.

Each file consists of a number of header lines followed by a number of alignment lines.
Each line can be split into a number of fields separated by TAB (`\t`) characters.

## Header lines

**NOT IMPLEMENTED**

Header lines start with character `@`.

## Alignment lines

Each alignment line has 12 mandatory fields.
Missing values in fields 3 to 11 are indicated by character `*`.

|Field|Type  |Description|
|----:|:----:|:----------|
|1    |string|Query sequence name|
|2    |int   |Query sequence length|
|3    |int   |Query start (0-based; closed)|
|4    |int   |Query end (0-based; open)|
|5    |char  |Strand relative to the path; always `+`|
|6    |string|Target path represented as a GFA walk|
|7    |int   |Target path length|
|8    |int   |Start position on the target path (0-based; closed)|
|9    |int   |End position on the target path (0-based; open)|
|10   |int   |Number of matches|
|11   |int   |Number of matches, mismatches, insertions, and deletions|
|12   |int   |Mapping quality (0-255; 255 for missing)|

**Example:**
```txt
read1 	6 	0 	6 	+ 	>2>3>4 	12 	2 	8 	6 	6 	60 	cs:Z::6
read2 	7 	0 	7 	+ 	>2>5>6 	11 	1 	8 	7 	7 	60 	cs:Z::7
read3 	7 	0 	7 	* 	* 	* 	* 	* 	* 	* 	255 	cs:Z:+GATTACA
```

### Query sequence name

Query sequence names must follow SAM conventions.
A name may contain any printable ASCII characters in the range `[!-~]`, except `@`.
This allows distinguishing header lines from alignment lines.

### Target path

This version of GAF does not allow specifying the target path using stable rGFA coordinates or nodes (GFA segments) with string names.
Nodes must have positive integer identifiers.
Node identifier `0` cannot be used, as many graph implementations reserve it for technical purposes.

### Optional fields

Optional fields are stored in the SAM-style `TAG:TYPE:VALUE` format.
The tag is a two-character string matching `/[A-Za-z][A-Za-z0-9]/`.
No tag can appear more than once on the same line, and the order of the optional fields does not matter.

The following types are currently supported:

|Type|Description|
|:--:|:----------|
|`A` |Printable character in `[!-~]`|
|`Z` |String of printable characters and spaces (`[ !-~]*`)|
|`i` |Signed 64-bit integer|
|`f` |Double-precision floating point number|
|`b` |Boolean value, with `1` for true and `0` for false|

### Difference string

Difference strings represent an edit script that transforms the given interval of the target path to the given interval of the query sequence.
They are stored as an optional field `cs` of type `Z`.
We support a subset of the operations defined for minimap2 difference strings.

|Operation|Regex   |Description|
|:-------:|:------:|:----------|
|`:`      |`[0-9]+`|Number of matching bases|
|`*`      |`[ACGTN][ACGTN]`|Mismatch as (target base, query base)|
|`+`      |`[ACGTN]+`|Insertion as the unaligned query bases|
|`-`      |`[ACGTN]+`|Deletion as the unaligned target bases|

### Other defined optional fields

|Tag |Type|Description|
|:--:|:--:|:----------|
|`AS`|`i` |Alignment score|
|`bq`|`Z` |Base quality string; must have the same length as the query sequence|
|`fn`|`Z` |Name of the next fragment (for paired alignments; cannot be used with `fp`)|
|`fp`|`Z` |Name of the previous fragment (for paired alignments; cannot be used with `fn`)|
|`pd`|`b` |This alignment and its pair (specified by `fn` or `fp`) are properly paired|
|`fi`|`i` |Fragment identifier for a fragmented alignment (see below)|

## Conventions

### Primary alignments

**NOT IMPLEMENTED**

A primary alignment represents an alignment of the entire query sequence to a non-empty interval of a target path.
Query start (field 3) must be `0` and query end (field 4) must have the same value as query sequence length (field 2).
A difference string must be present to allow recovering the entire query sequence.

### Unaligned sequences

**NOT IMPLEMENTED**

An unaligned sequence is represented as an alignment of the entire query sequence to a missing interval of a missing target path.
Query start (field 3) must be `0` and query end (field 4) must have the same value as query sequence length (field 2).
A difference string must be present, with the entire query sequence as a single insertion, to allow recovering the sequence.

### Fragmented alignments

A fragmented alignment is a single alignment represented as number of alignment lines (e.g. corresponding to subpaths that are within a specific subgraph).
The fragments (alignment lines) correspond to non-overlapping intervals of the underlying alignment.
Each fragment represents an alignment of a non-empty query interval to a non-empty interval of a non-empty target path.

Query interval (fields 3 and 4), target path (fields 6 to 9), and the difference string must be specific to each fragment.
Other fields are inherited from the underlying alignment.
Fragments are identified by fragment indexes starting from `1`, stored as an optional field `fi` of type `i`.
