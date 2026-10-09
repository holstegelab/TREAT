"""Resolve BED chromosome names against the reference and all input BAMs."""

import pysam


def chromosomeSources(inBam, ref):
    sources = []
    with pysam.FastaFile(ref) as fasta:
        sources.append(('reference FASTA', ref, set(fasta.references)))
    for bam_path in inBam:
        with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as bam:
            sources.append(('BAM', bam_path, set(bam.references)))
    return sources


def validateChromosomes(bed, inBam, ref):
    required = set(bed)
    sources = chromosomeSources(inBam, ref)
    errors = []
    for kind, path, available in sources:
        missing = sorted(required - available)
        if missing:
            errors.append('%s %s is missing BED chromosome name(s): %s' %
                          (kind, path, ', '.join(missing)))
    if errors:
        raise ValueError('\n'.join(errors) +
                         '\nBED, BAM, and FASTA chromosome names must match exactly '
                         '(for example, 4 and chr4 are different). '
                         'TREAT preserves the names supplied in the BED.')


def resolveChromosomes(bed, inBam, ref, bed_path):
    sources = chromosomeSources(inBam, ref)
    common = set.intersection(*(available for _, _, available in sources))
    mapping = {}
    errors = []
    for chrom in bed:
        alternate = chrom[3:] if chrom.startswith('chr') else 'chr' + chrom
        if chrom in common:
            mapping[chrom] = chrom
        elif alternate in common:
            mapping[chrom] = alternate
        else:
            missing = ['%s %s (missing %s)' %
                       (kind, path, ', '.join(name for name in [chrom, alternate]
                                             if name not in available))
                       for kind, path, available in sources
                       if chrom not in available or alternate not in available]
            errors.append('Cannot resolve BED chromosome %s: neither %s nor %s '
                          'exists consistently in the FASTA and every BAM.\n%s' %
                          (chrom, chrom, alternate, '\n'.join(missing)))
    if errors:
        raise ValueError('\n'.join(errors) +
                         '\nBAM and FASTA chromosome names must agree. '
                         'Only a matching chr-prefix alternative can be resolved automatically.')
    resolved = {}
    for chrom, intervals in bed.items():
        target = mapping[chrom]
        resolved.setdefault(target, []).extend(
            [[start, end, '%s:%s-%s' % (target, start, end)]
             for start, end, _ in intervals])
    if any(chrom != target for chrom, target in mapping.items()):
        with open(bed_path, 'w') as output:
            for chrom, intervals in resolved.items():
                for start, end, _ in intervals:
                    output.write('%s\t%s\t%s\n' % (chrom, start, end))
        for chrom, target in mapping.items():
            if chrom != target:
                print('** BED chromosome name resolved: %s -> %s' % (chrom, target))
    return resolved
