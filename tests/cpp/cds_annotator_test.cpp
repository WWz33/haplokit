#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include "cds_annotator.h"

namespace {

const char* kGff =
    "c1\tsrc\tgene\t1000\t3000\t.\t+\t.\tID=gA\n"
    "c1\tsrc\tmRNA\t1000\t3000\t.\t+\t.\tID=gA.1;Parent=gA\n"
    "c1\tsrc\tCDS\t1000\t1300\t.\t+\t0\tID=gA.1.cds1;Parent=gA.1\n"
    "c1\tsrc\tCDS\t2000\t2200\t.\t+\t2\tID=gA.1.cds2;Parent=gA.1\n"
    // Minus-strand gene: two CDS spans, spliced CDS runs right-to-left.
    "c1\tsrc\tgene\t5000\t7500\t.\t-\t.\tID=gB\n"
    "c1\tsrc\tmRNA\t5000\t7500\t.\t-\t.\tID=gB.1;Parent=gB\n"
    "c1\tsrc\tCDS\t5000\t5300\t.\t-\t2\tID=gB.1.cds2;Parent=gB.1\n"
    "c1\tsrc\tCDS\t7000\t7300\t.\t-\t0\tID=gB.1.cds1;Parent=gB.1\n"
    // Overlapping gene on the plus strand sharing gB's span range.
    "c1\tsrc\tgene\t5000\t7500\t.\t+\t.\tID=gC\n"
    "c1\tsrc\tmRNA\t5000\t7500\t.\t+\t.\tID=gC.1;Parent=gC\n"
    "c1\tsrc\tCDS\t5100\t5200\t.\t+\t0\tID=gC.1.cds1;Parent=gC.1\n";

// gA.1 CDS (plus): span1 301nt + span2 201nt = 502nt.
// Layout plan (1-based CDS coordinates):
//   c.1-3    ATG   start codon
//   c.4-6    GGG   Ala (positions 1003-1005 genomic span1)
//   c.10-12  TAC   Tyr
//   c.499-501 TAA  stop (genomic 2197-2199, span2)
// gB.1 CDS (minus): spliced from 7300 down to 7000 (301nt) then 5300 down to 5000 (301nt), total 602nt.
//   c.1-3 comes from genomic 7300,7299,7298 (complemented).
// gC.1 CDS (plus): 101nt, c.1-3 = genomic 5100-5102.

std::string build_genome() {
    std::string genome(8000, 'A');
    const auto set = [&](int pos, const std::string& seq) {
        for (std::size_t i = 0; i < seq.size(); ++i) {
            genome[pos - 1 + i] = seq[i];
        }
    };
    // Fill the gA spans with an arbitrary base so only our forced codons matter.
    for (int pos = 1000; pos <= 1300; ++pos) genome[pos - 1] = 'C';
    for (int pos = 2000; pos <= 2200; ++pos) genome[pos - 1] = 'C';
    // c.1-3 ATG (span1 start)
    set(1000, "ATG");
    // c.4-6 CAG (Gln); a single c.4 C>T makes TAG, a stop_gained
    set(1003, "CAG");
    // c.10-12 TAC -> Tyr (codon 4)
    set(1009, "TAC");
    // stop at c.499-501: span1 is 301nt, so c.499 = span2 pos 2000+198-301-1 = ...
    // c.302 = span2 first base (2000). c.499 -> span2 index 499-302 = 197 -> 2197.
    set(2197, "TAA");
    // minus-strand gene gB: CDS runs 7300..7000 then 5300..5000, complemented.
    for (int pos = 7000; pos <= 7300; ++pos) genome[pos - 1] = 'G';  // complement 'C'
    for (int pos = 5000; pos <= 5300; ++pos) genome[pos - 1] = 'G';
    // c.1-3: genomic 7300,7299,7298 complemented. Make it ATG:
    // complement(A)=T so genomic should be T,A? genomic 7300='T'(->A), 7299='A'(->T), 7298='G'(->C)? need ATG: A T G
    // c.1 = comp(genome[7300]) = 'A' => genome[7300]='T'
    // c.2 = comp(genome[7299]) = 'T' => genome[7299]='A'
    // c.3 = comp(genome[7298]) = 'G' => genome[7298]='C'
    genome[7300 - 1] = 'T'; genome[7299 - 1] = 'A'; genome[7298 - 1] = 'C';
    // c.4-6 (genomic 7297..7295): want GCC -> Ala: genome 7297='C',7296='G',7295='G'
    genome[7297 - 1] = 'C'; genome[7296 - 1] = 'G'; genome[7295 - 1] = 'G';
    // Codon 200 covers c.598-600 (the last full codon of the 602nt CDS; c.601-602
    // are a trailing partial codon). c.598=comp(genome[5004]), c.599=comp(genome[5003]),
    // c.600=comp(genome[5002]). Want TAA (stop):
    // c.598='T' => genome[5004]='A'; c.599='A' => genome[5003]='T'; c.600='A' => genome[5002]='T'
    genome[5004 - 1] = 'A'; genome[5003 - 1] = 'T'; genome[5002 - 1] = 'T';
    // gC.1: span 5100-5200 plus.
    for (int pos = 5100; pos <= 5200; ++pos) genome[pos - 1] = 'T';
    set(5100, "ATG");
    return genome;
}

const char* complement(const std::string& in) {
    static std::string out;
    out.clear();
    for (auto it = in.rbegin(); it != in.rend(); ++it) {
        switch (*it) {
            case 'A': out += 'T'; break;
            case 'T': out += 'A'; break;
            case 'G': out += 'C'; break;
            case 'C': out += 'G'; break;
            default: out += 'N';
        }
    }
    return out.c_str();
}

bool write_fixtures(const std::string& genome, const char* gff_path, const char* fa_path) {
    {
        std::ofstream g(gff_path);
        g << kGff;
    }
    {
        std::ofstream f(fa_path);
        // Headers follow gffread -x output: ">ID loc:seqid:start-end,... attrs".
        // Only the first whitespace-delimited token is used as the transcript ID.
        // gA.1: span1 (1000-1300) + span2 (2000-2200)
        f << ">gA.1 loc:c1:1000-1300,2000-2200 ID=gA.1.cds1;Parent=gA.1\n"
          << genome.substr(1000 - 1, 301) << genome.substr(2000 - 1, 201) << "\n";
        // gB.1: minus strand, spliced = revcomp(7000-7300) + revcomp(5000-5300)
        f << ">gB.1 loc:c1:5000-5300,7000-7300 ID=gB.1.cds2;Parent=gB.1\n"
          << complement(genome.substr(7000 - 1, 301))
          << complement(genome.substr(5000 - 1, 301)) << "\n";
        // gC.1: plus span 5100-5200, plain header (no suffix) also accepted
        f << ">gC.1\n" << genome.substr(5100 - 1, 101) << "\n";
    }
    return true;
}

int failures = 0;

void check(bool ok, const std::string& what) {
    if (!ok) {
        std::cerr << "FAIL: " << what << "\n";
        ++failures;
    } else {
        std::cout << "ok: " << what << "\n";
    }
}

int run() {
    const std::string genome = build_genome();
    write_fixtures(genome, "data/cds_test_tmp.gff3", "data/cds_test_tmp.fa");

    haplokit::CdsAnnotator annotator;
    if (!annotator.load("data/cds_test_tmp.gff3", "data/cds_test_tmp.fa")) {
        std::cerr << "load failed: " << annotator.last_error() << "\n";
        return 1;
    }

    // --- plus strand, start codon: genomic 1000 = c.1 A, A>G start_lost ---
    {
        auto hits = annotator.annotate_snp("c1", 1000, "A", "G");
        check(hits.size() == 1, "start-codon snv hits gA.1 only");
        if (hits.size() == 1) {
            check(hits[0].transcript_id == "gA.1", "start snv transcript");
            check(hits[0].cds_pos == 1, "start snv cds_pos=1");
            check(hits[0].effect == "start_lost" && hits[0].protein_change == "p.Met1?",
                  "start snv effect+protein p.Met1?");
        }
    }

    // --- plus strand, synonymous: codon c.4-6 CAG (Gln); c.6 G->A gives CAA,
    // still Gln. c.6 = genomic 1005.
    {
        auto hits = annotator.annotate_snp("c1", 1005, "G", "A");
        check(hits.size() == 1 && hits[0].effect == "synonymous_variant" &&
                  hits[0].protein_change == "p.(=)",
              "synonymous CAG>CAA");
        if (!hits.empty()) check(hits[0].cds_pos == 6, "syn cds_pos=6");
    }

    // --- plus strand, missense: c.10-12 TAC (Tyr). c.10 T->C -> CAC (His). genomic 1009.
    {
        auto hits = annotator.annotate_snp("c1", 1009, "T", "C");
        check(hits.size() == 1, "missense hits");
        if (hits.size() == 1) {
            check(hits[0].cds_pos == 10 && hits[0].codon_index == 4, "missense cds_pos/codon_index");
            check(hits[0].effect == "missense_variant", "missense effect");
            check(hits[0].protein_change == "p.Tyr4His", "missense protein change");
        }
    }

    // --- plus strand, stop_gained: codon c.4-6 CAG (Gln); c.4 C->T -> TAG.
    {
        auto hits = annotator.annotate_snp("c1", 1003, "C", "T");
        check(hits.size() == 1 && hits[0].effect == "stop_gained" &&
                  hits[0].protein_change == "p.Gln2Ter",
              "stop_gained CAG>TAG");
    }

    // --- plus strand, stop_lost: c.499-501 TAA stop. c.499 T->C -> CAA (Gln). genomic 2197.
    {
        auto hits = annotator.annotate_snp("c1", 2197, "T", "C");
        check(hits.size() == 1 && hits[0].effect == "stop_lost" &&
                  hits[0].cds_pos == 499 &&
                  hits[0].protein_change == "p.Ter167GlnextTer?",
              "stop_lost TAA>CAA");
    }

    // --- plus strand, stop_retained: TAA>TAG at c.501 (genomic 2199).
    {
        auto hits = annotator.annotate_snp("c1", 2199, "A", "G");
        check(hits.size() == 1 && hits[0].effect == "stop_retained_variant" &&
                  hits[0].protein_change == "p.Ter167=",
              "stop_retained TAA>TAG");
    }

    // --- outside CDS (upstream of gA, between spans: 1500 intron) ---
    {
        auto hits = annotator.annotate_snp("c1", 1500, "A", "T");
        check(hits.empty(), "intron snv no hits");
    }

    // --- minus strand, start codon: genomic 7300 = c.1. genome[7300]='T', ref on
    // transcript strand = comp('T') = 'A'. Change genome T->C => transcript A->G, start_lost.
    {
        auto hits = annotator.annotate_snp("c1", 7300, "T", "C");
        check(hits.size() == 1, "minus start snv hits gB.1");
        if (hits.size() == 1) {
            check(hits[0].transcript_id == "gB.1", "minus start transcript");
            check(hits[0].cds_pos == 1 && hits[0].ref_base == 'A' && hits[0].alt_base == 'G',
                  "minus start bases complemented");
            check(hits[0].effect == "start_lost" && hits[0].protein_change == "p.Met1?",
                  "minus start_lost p.Met1?");
        }
    }

    // --- minus strand, span2 (genomic 5300 end of the second-in-genome-order span):
    // c.302 = comp(genome[5300]). genome[5300]='G' -> transcript 'C'.
    {
        auto hits = annotator.annotate_snp("c1", 5300, "G", "A");
        check(hits.size() == 1 && hits[0].cds_pos == 302, "minus span2 cds_pos=302");
        if (!hits.empty()) check(hits[0].ref_base == 'C' && hits[0].alt_base == 'T', "minus span2 bases");
    }

    // --- minus strand, stop codon c.598-600 (codon 200): stop_lost. c.600 is
    // genomic 5002, transcript base comp('T')='A'; A>G -> CAA (Gln), stop lost.
    {
        auto hits = annotator.annotate_snp("c1", 5002, "T", "A");
        check(hits.size() == 1 && hits[0].effect == "stop_lost" && hits[0].cds_pos == 600 &&
                  hits[0].protein_change == "p.Ter200TyrextTer?",
              "minus stop_lost at c.600");
    }

    // --- multi-gene overlap: genomic 5150 is inside gB.1 (minus) and gC.1 (plus) CDS ---
    {
        auto hits = annotator.annotate_snp("c1", 5150, "T", "G");
        check(hits.size() == 2, "overlap snv hits both transcripts");
        if (hits.size() == 2) {
            check(hits[0].transcript_id == "gB.1" && hits[1].transcript_id == "gC.1",
                  "overlap sorted by transcript id");
        }
    }

    // --- REF mismatch: pretend REF is wrong; annotator must skip ---
    {
        auto hits = annotator.annotate_snp("c1", 1000, "C", "G");
        check(hits.empty(), "ref mismatch skipped");
    }

    // --- indel / multi-base skipped ---
    {
        auto hits = annotator.annotate_snp("c1", 1000, "AT", "G");
        check(hits.empty(), "indel not annotated");
    }

    // --- missing transcript in FASTA: gC.1 exists; query a transcript absent from fasta
    // is not directly reachable via overlap, so cover the fasta-lookup miss via gC span
    // with a mismatching header by loading a second annotator.
    {
        haplokit::CdsAnnotator ann2;
        {
            std::ofstream f("data/cds_test_tmp2.fa");
            f << ">only_gA\n" << genome.substr(1000 - 1, 301) << genome.substr(2000 - 1, 201) << "\n";
        }
        check(ann2.load("data/cds_test_tmp.gff3", "data/cds_test_tmp2.fa"), "second fixture loads");
        auto hits = ann2.annotate_snp("c1", 1009, "T", "C");
        check(hits.empty(), "transcript absent from fasta skipped");
        std::remove("data/cds_test_tmp2.fa");
    }

    check(failures == 0, "all checks passed");

    std::remove("data/cds_test_tmp.gff3");
    std::remove("data/cds_test_tmp.fa");
    return failures == 0 ? 0 : 1;
}

}  // namespace

int main() {
    try {
        return run();
    } catch (const std::exception& e) {
        std::cerr << "exception: " << e.what() << "\n";
        return 1;
    }
}
