# FASTQC Results

## HER
- **HER21 R1**: total seq 61,247,419, seq length 60 bp, phred score good at the start of the sequence then drop toward the end with wiskers below 20 from bp 47 onwards, average quality phred score mostly between 33 and 39, shakey base content distribution for first 10 bp, bias in GC content, high duplication rate with only 29.05% remaining after deduplication, overrepresented sequences without identified source, adopter content low 0.7% towards the end of the reads
- **HER21 R2**: total seq 61,247,419, seq length 60 bp, phred score good at the start of the sequence then drop toward the end with wiskers below 20 from bp 38 onwards -> drop in quality score starts earlier than for R1, average quality phred score mostly between 33 and 39, shakey base content distribution for first 10 bp, bias in GC content, high duplication rate with only 31.51% remaining after deduplication, overrepresented sequences without identified source, adopter content low 0.7% towards the end of the reads
- **HER22 R1**: total seq 68,888,018, seq length 60 bp, phred score good at the start of the sequence then drop toward the end with wiskers below 20 from bp 48 onwards, average quality phred score mostly between 33 and 39 (majority at 33), shakey base content distribution for first 10 bp, bias in GC content, high duplication rate with only 35.19% remaining after deduplication, overrepresented sequences without identified source, adopter content low
- **HER22 R2**: total seq 68,888,018, seq length 60 bp, phred score good at the start of the sequence then drop toward the end with wiskers below 20 from bp **26** onwards, drop is drastic, average quality phred score mostly between 33 and 38 (majority at 33), shakey base content distribution for first 10 bp, bias in GC content, little amount of base N content present, high duplication rate with only 39.3% remaining after deduplication, overrepresented sequences without identified source (one sequence is NNNNN), adopter content low
- **HER23 R1**: total seq 52,010,599, seq length 60 bp, phred score good at the start of the sequence then drop toward the end but generally quite good, average quality phred score mostly between 33 and 39, shakey base content distribution for first 10 bp, bias in GC content, high duplication rate with only 32.4% remaining after deduplication, overrepresented sequences without identified source, adopter content low
- **HER23 R1**: total seq 52,010,599, seq length 60 bp, phred score good at the start of the sequence then drop toward the end with wiskers below 20 from bp 55, average quality phred score mostly between 33 and 39, shakey base content distribution for first 10 bp, bias in GC content, high duplication rate with only 36.17% remaining after deduplication, overrepresented sequences without identified source, adopter content low

**Summary**: drop in quality towards the end of the reads, two peaks in per sequence quality score, high dublication rate, bias in GC content, presence of overrepresented sequences

&rarr; reasons for the two peaks: adapter or primaer contamination - contaminated reads tend to have lower quality scores, low quality bases at reads end, multiple sequencing batches/ lanes, contaminant sequences from other organisms, heterogenous library preparation




