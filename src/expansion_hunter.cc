//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <boost/algorithm/string/join.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "third_party/spdlog/fmt/ostr.h"
#include "third_party/spdlog/spdlog.h"

#include "classification/alignment_summary.h"
#include "classification/mapping_classifier.h"
#include "classification/overlap_quantification.h"
#include "common/parameters.h"
#include "common/ref_genome.h"
#include "common/seq_operations.h"
#include "common/timestamp.h"
#include "genotyping/repeat_genotyper.h"
#include "graphs/gapless_aligner.h"
#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "graphs/graph_mapping_operations.h"
#include "graphs/path.h"
#include "include/bam_file.h"
#include "include/bam_index.h"
#include "include/irr_counting.h"
#include "include/json_output.h"
#include "include/read_group.h"
#include "include/region_findings.h"
#include "include/vcf_output.h"
#include "include/version.h"
#include "purity/purity.h"
#include "reads/aligned_reader.h"
#include "reads/read_operations.h"
#include "reads/read_pairs.h"
#include "region_analysis/region_analysis.h"
#include "rep_align/rep_align.h"
#include "stats/counts.h"

namespace spd = spdlog;

using std::array;
using std::cerr;
using std::endl;
using std::list;
using std::map;
using std::pair;
using std::string;
using std::unordered_set;
using std::vector;

// Returns the length of the first read in a BAM file.
size_t CalcReadLen(const string &bam_path) {
  // Open a BAM file for reading.
  samFile *file_ptr = sam_open(bam_path.c_str(), "r");
  if (!file_ptr) {
    throw std::runtime_error("Failed to read BAM file '" + bam_path + "'");
  }
  bam_hdr_t *header_ptr = sam_hdr_read(file_ptr);
  if (!header_ptr) {
    throw std::runtime_error("BamFile::Init: Failed to read BAM header: '" +
                             bam_path + "'");
  }

  enum { kSupplimentaryAlign = 0x800, kSecondaryAlign = 0x100 };

  size_t read_len = 99;
  bam1_t *align_ptr = bam_init1();
  int ret;
  while ((ret = sam_read1(file_ptr, header_ptr, align_ptr)) >= 0) {
    const bool is_supplimentary = align_ptr->core.flag & kSupplimentaryAlign;
    const bool is_secondary = align_ptr->core.flag & kSecondaryAlign;
    const bool is_primary_align = (!is_supplimentary) && (!is_secondary);
    if (is_primary_align) {
      read_len = align_ptr->core.l_qseq;
      break;
    }
  }

  if (ret < 0) {
    throw std::runtime_error("Failed to extract a read from BAM file");
  }

  bam_destroy1(align_ptr);
  bam_hdr_destroy(header_ptr);
  sam_close(file_ptr);

  return read_len;
}

void TestFlankingReads(const Parameters &parameters,
                       const vector<RepeatAlign> &flanking_repaligns) {
  int32_t num_left_flanking = 0;
  int32_t num_right_flanking = 0;
  for (const RepeatAlign &repalign : flanking_repaligns) {
    if (repalign.left_flank_len != 0) {
      assert(repalign.right_flank_len == 0);
      ++num_left_flanking;
    } else {
      assert(repalign.left_flank_len == 0 && repalign.right_flank_len != 0);
      ++num_right_flanking;
    }
  }

  const int32_t read_length = parameters.read_len();
  const double depth = parameters.depth();

  const double start_probability = depth / read_length;
  const int32_t num_possible_starts = read_length - 2;
  const double expected_num_flanking_reads =
      start_probability * num_possible_starts;

  assert(num_left_flanking <= num_possible_starts);
  assert(num_right_flanking <= num_possible_starts);

  ExpectedCountTest count_test((int32_t)expected_num_flanking_reads);

  cerr << "\t[Flanking reads: " << flanking_repaligns.size() << "]" << endl;
  cerr << "\t[Left-flanking reads: " << num_left_flanking
       << " p-value: " << count_test.Test(num_left_flanking) << "]" << endl;
  cerr << "\t[Right-flanking reads: " << num_right_flanking
       << " p-value: " << count_test.Test(num_right_flanking) << "]" << endl;
}

// Search for reads spanning the entire repeat sequence.
void FindShortRepeats(const Parameters &parameters, BamFile &bam_file,
                      const RepeatSpec &repeat_spec, AlignPairs &align_pairs,
                      vector<RepeatReadGroup> &read_groups,
                      vector<RepeatAlign> *flanking_repaligns) {
  const int unit_len = repeat_spec.units_shifts[0][0].length();

  map<int, vector<RepeatAlign>> size_spanning_repaligns;
  flanking_repaligns->clear();

  // Graph test.
  GraphSharedPtr graph_ptr =
      MakeStrGraph(repeat_spec.left_flank, repeat_spec.units_shifts[0][0],
                   repeat_spec.right_flank);
  const int32_t kmer_len = 14;
  StrandClassifier strand_classifier(graph_ptr, kmer_len);
  GaplessAligner aligner(graph_ptr, kmer_len);
  StrMappingClassifier mapping_classifier(0, 1, 2);
  // End graph test.

  // Align each read to the repeat.
  for (auto &kv : align_pairs) {
    AlignPair &frag = kv.second;

    for (Align &align : frag) {
      RepeatAlign rep_align;
      const bool aligns = AlignRead(parameters, repeat_spec, align.bases,
                                    align.quals, &rep_align);
      rep_align.read.name = align.name;

      cerr << " /" << endl;
      cerr << "/" << endl;
      {
        cerr << align.name << endl;
        const string read_seq = strand_classifier.IsForwardOriented(align.bases)
                                    ? align.bases
                                    : ReverseComplement(align.bases);
        list<GraphMapping> mappings = aligner.GetBestAlignment(read_seq);
        GraphMapping canonical_mapping =
            mapping_classifier.GetCanonicalMapping(mappings);
        cerr << read_seq << " "
             << mapping_classifier.Classify(canonical_mapping) << endl;
        cerr << canonical_mapping << endl;
      }

      cerr << align.bases << " ";

      if (!aligns) {
        cerr << "not_flanking_or_spanning" << endl;
        cerr << "\\ " << endl;
        cerr << " \\" << endl;
        continue;
      }

      if (rep_align.type == RepeatAlign::Type::kSpanning) {
        cerr << "spanning" << endl;
      } else {
        assert(rep_align.type == RepeatAlign::Type::kFlanking);
        cerr << "flanking" << endl;
      }
      cerr << "\\ " << endl;
      cerr << " \\" << endl;

      // Not pretty, but lets downstream code know that this is not an
      // IRR.
      align.status = kFlankingRead;

      if (rep_align.type == RepeatAlign::Type::kSpanning) {
        size_spanning_repaligns[rep_align.size].push_back(rep_align);
      } else {
        assert(rep_align.type == RepeatAlign::Type::kFlanking);
        flanking_repaligns->push_back(rep_align);
      }
    }
  }

  for (const auto &size_repaligns : size_spanning_repaligns) {
    RepeatReadGroup read_group;
    read_group.read_type = ReadType::kSpanning;
    read_group.size = size_repaligns.first;
    read_group.num_supporting_reads = size_repaligns.second.size();
    read_group.rep_aligns = size_repaligns.second;
    read_groups.push_back(read_group);
  }

  DistributeFlankingReads(parameters, repeat_spec, &read_groups,
                          flanking_repaligns);

  const double haplotype_depth = parameters.depth() / 2;
  CoalesceFlankingReads(repeat_spec, read_groups, flanking_repaligns,
                        parameters.read_len(), haplotype_depth,
                        repeat_spec.units[0].length(), repeat_spec.units_shifts,
                        parameters.min_baseq(), parameters.min_wp());
}

// Caches alignments from extended target and off-target regions. For BAM files,
// mates of reads in the relevant regions are cached too.
void CacheAligns(BamFile *bam_file, const RepeatSpec &repeat_spec,
                 AlignPairs &align_pairs,
                 unordered_set<string> &ontarget_frag_names,
                 int extension_len) {
  align_pairs.clear();
  Region extended_target_region =
      repeat_spec.target_region.Extend(extension_len);

  // Cache on-target reads.
  CacheReadsFromRegion(extended_target_region, kCacheAll,
                       repeat_spec.units_shifts, 0.9, &(*bam_file),
                       &align_pairs);

  // Save names of fragments that overlap target locus before caching reads
  // from confusion regions.
  ontarget_frag_names.clear();
  for (const auto &kv : align_pairs) {
    const AlignPair &frag = kv.second;
    const string name = !frag[0].name.empty() ? frag[0].name : frag[1].name;
    ontarget_frag_names.insert(name);
  }

  cerr << TimeStamp() << ",\t[Found " << ontarget_frag_names.size()
       << " reads in target locus]" << endl;

  // Cache aligned off-target reads.
  for (const Region &confusion_region : repeat_spec.offtarget_regions) {
    Region extended_confusion_region = confusion_region.Extend(extension_len);
    CacheReadsFromRegion(extended_confusion_region, kCacheAll,
                         repeat_spec.units_shifts, 0.9, &(*bam_file),
                         &align_pairs);
  }

  // Filling-in missing mates by jumping around the BAM.
  cerr << TimeStamp() << ",\t[Filling in mates]" << endl;
  if ((*bam_file).format() == BamFile::kBamFile) {
    FillinMates(*bam_file, align_pairs, repeat_spec.units_shifts, 0.9,
                ontarget_frag_names);
  } else {
    cerr << TimeStamp() << ",\t[Skipping filling in mates]" << endl;
  }
  cerr << TimeStamp() << ",\t[Done filling in mates]" << endl;
}

bool IsFlannkingGroup(const RepeatReadGroup &read_group) {
  return read_group.read_type == ReadType::kFlanking;
}

bool FindLongRepeats(const Parameters &parameters,
                     const RepeatSpec &repeat_spec, BamFile &bam_file,
                     unordered_set<string> &ontarget_frag_names,
                     AlignPairs &align_pairs, RegionFindings &region_findings) {
  // Count the number of anchored IRRs.
  region_findings.num_anchored_irrs = 0;
  const int extension_len = parameters.region_extension_len();
  Region target_nhood = repeat_spec.target_region.Extend(extension_len);

  // Anchored IRRs and IRR pairs will be stored here.
  vector<RepeatAlign> irr_rep_aligns;
  CountAnchoredIrrs(bam_file, parameters, target_nhood, repeat_spec,
                    ontarget_frag_names, align_pairs,
                    region_findings.num_anchored_irrs, repeat_spec.units_shifts,
                    &irr_rep_aligns);

  cerr << TimeStamp() << ",\t[Found " << region_findings.num_anchored_irrs
       << " anchored IRRs]" << endl
       << TimeStamp() << ",\t[Cached " << align_pairs.size() << " reads]"
       << endl;

  // Stores the total IRR count.
  region_findings.num_irrs = region_findings.num_anchored_irrs;

  // Look for IRR pairs only if anchored IRRs are found.
  if (region_findings.num_anchored_irrs) {
    if (!repeat_spec.is_common_unit()) {
      cerr << TimeStamp() << ",\t[Counting aligned IRR pairs]" << endl;
      map<string, int> numIrrConfRegion;
      region_findings.num_irrs +=
          CountAlignedIrr(bam_file, parameters, align_pairs, numIrrConfRegion,
                          repeat_spec.units_shifts, &irr_rep_aligns);

      // Record paired IRR counts from each confusion region.
      region_findings.offtarget_irr_counts.clear();
      for (const Region &confusionRegion : repeat_spec.offtarget_regions) {
        Region confusionNhood = confusionRegion.Extend(extension_len);
        region_findings.offtarget_irr_counts.push_back(
            numIrrConfRegion[confusionNhood.ToString()]);
      }

      region_findings.num_unaligned_irrs = 0;

      if (!parameters.skip_unaligned()) {
        cerr << TimeStamp() << ",\t[Counting unaligned IRRs]" << endl;
        CountUnalignedIrrs(bam_file, parameters,
                           region_findings.num_unaligned_irrs,
                           repeat_spec.units_shifts, &irr_rep_aligns);
        region_findings.num_irrs += region_findings.num_unaligned_irrs;
      } else {
        cerr << TimeStamp() << ",\t[Skipping unaligned IRRs]" << endl;
      }
    }

    const size_t unit_len = repeat_spec.units[0].length();

    RepeatReadGroup read_group;
    read_group.read_type = ReadType::kInrepeat;

    read_group.size =
        (int)(std::ceil(parameters.read_len() / (double)unit_len));
    read_group.num_supporting_reads = region_findings.num_irrs;
    read_group.rep_aligns = irr_rep_aligns;

    region_findings.read_groups.push_back(read_group);

    // If there is evidence for a long repeat allele we assume that flanking
    // reads came from and so any previously-found alleles whose size was
    // estimated from flanking reads are no longer valid.

    vector<RepeatReadGroup>::const_iterator flanking_allele_it =
        std::find_if(region_findings.read_groups.begin(),
                     region_findings.read_groups.end(), IsFlannkingGroup);

    if (flanking_allele_it != region_findings.read_groups.end()) {
      region_findings.flanking_repaligns.insert(
          region_findings.flanking_repaligns.end(),
          flanking_allele_it->rep_aligns.begin(),
          flanking_allele_it->rep_aligns.end());
    }

    region_findings.read_groups.erase(
        std::remove_if(region_findings.read_groups.begin(),
                       region_findings.read_groups.end(), IsFlannkingGroup),
        region_findings.read_groups.end());
  }
}

static bool CompareRegionFindings(const RegionFindings &r1,
                                  const RegionFindings &r2) {
  return r1.region < r2.region;
}

// For each repeat region, calculate sizes of repeats that are (a) shorter and
// (b) longer than the read length; reconsile alleles and output results to
// appropriate files.

void EstimateRepeatSizes(const Parameters &parameters,
                         const map<string, RepeatSpec> &repeat_specs,
                         BamFile *bam_file, Outputs *outputs) {
  // boost::property_tree::ptree ptree_root;
  vector<RegionFindings> sample_findings;

  // Analyze repeats one by one.
  for (auto &kv : repeat_specs) {
    const RepeatSpec &repeat_spec = kv.second;

    string repeat_header = repeat_spec.target_region.ToString();
    repeat_header += " " + boost::algorithm::join(repeat_spec.units, "/");

    GenotypeType genotype_type = GenotypeType::kDiploid;
    const string chrom = repeat_spec.target_region.chrom();

    const bool is_female_chrom_y =
        parameters.sex() == Sex::kFemale && (chrom == "chrY" || chrom == "Y");

    if (is_female_chrom_y) {
      cerr << TimeStamp() << ",[Skipping " << repeat_header
           << " because the sample is female]" << endl;
      continue;
    }

    cerr << TimeStamp() << ",[Analyzing " << repeat_header << "]" << endl;

    const bool is_sex_chrom =
        chrom == "chrX" || chrom == "X" || chrom == "chrY" || chrom == "Y";

    if (parameters.sex() == Sex::kMale && is_sex_chrom) {
      genotype_type = GenotypeType::kHaploid;
    }

    RegionFindings region_findings;
    region_findings.region_id = repeat_spec.repeat_id;
    region_findings.region = repeat_spec.target_region;
    region_findings.offtarget_irr_counts =
        vector<int>(repeat_spec.offtarget_regions.size(), 0);

    cerr << TimeStamp() << ",\t[Caching reads]" << endl;
    AlignPairs align_pairs;
    unordered_set<string> ontarget_frag_names;
    CacheAligns(&(*bam_file), repeat_spec, align_pairs, ontarget_frag_names,
                parameters.region_extension_len());
    if (align_pairs.empty()) {
      cerr << TimeStamp() << ",\t[Found no on-target or off-target reads]"
           << endl;
      continue;
    }

    cerr << TimeStamp() << ",\t[Estimating short repeat sizes]" << endl;
    FindShortRepeats(parameters, *bam_file, repeat_spec, align_pairs,
                     region_findings.read_groups,
                     &region_findings.flanking_repaligns);

    cerr << TimeStamp() << ",\t[Estimating long repeat sizes]" << endl;
    region_findings.num_anchored_irrs = 0;
    region_findings.num_unaligned_irrs = 0;
    region_findings.num_irrs = 0;

    FindLongRepeats(parameters, repeat_spec, *bam_file, ontarget_frag_names,
                    align_pairs, region_findings);

    map<int, int> flanking_size_counts;
    map<int, int> spanning_size_counts;

    for (const auto &align : region_findings.flanking_repaligns) {
      flanking_size_counts[align.size] += 1;
    }

    const int num_units_in_read = (int)(std::ceil(
        parameters.read_len() / (double)repeat_spec.units[0].length()));
    int num_supporting_irr_reads = -1;

    vector<RepeatAllele> haplotype_candidates;
    // Add count of in-repeat reads to flanking.
    for (const auto &read_group : region_findings.read_groups) {
      if (read_group.read_type == ReadType::kSpanning) {
        spanning_size_counts[read_group.size] +=
            read_group.num_supporting_reads;
        haplotype_candidates.push_back(
            RepeatAllele(read_group.size, read_group.num_supporting_reads,
                         ReadType::kSpanning));
      } else if (read_group.read_type == ReadType::kInrepeat) {
        flanking_size_counts[num_units_in_read] +=
            read_group.num_supporting_reads;
        haplotype_candidates.push_back(
            RepeatAllele(num_units_in_read, read_group.num_supporting_reads,
                         ReadType::kInrepeat));
      } else if (read_group.read_type == ReadType::kFlanking) {
        haplotype_candidates.push_back(
            RepeatAllele(read_group.size, read_group.num_supporting_reads,
                         ReadType::kFlanking));
        for (const auto &align : read_group.rep_aligns) {
          flanking_size_counts[align.size] += 1;
        }
      } else {
        throw std::logic_error("Do not know how to deal with " +
                               kReadTypeToString.at(read_group.read_type) +
                               " alleles");
      }
    }

    cerr << TimeStamp() << ",\t[Flanking:";
    for (const auto &kv : flanking_size_counts) {
      cerr << " (" << kv.first << ", " << kv.second << ")";
    }
    cerr << "]" << endl;

    cerr << TimeStamp() << ",\t[Spanning:";
    for (const auto &kv : spanning_size_counts) {
      cerr << " (" << kv.first << ", " << kv.second << ")";
    }
    cerr << "]" << endl;

    cerr << TimeStamp() << ",\t[Haplotype candidates:";
    for (const auto &candiate : haplotype_candidates) {
      cerr << " " << candiate.size_;
    }
    cerr << "]" << endl;

    if (haplotype_candidates.empty()) {
      cerr
          << TimeStamp()
          << ",\t[Skipping this region because no informative reads were found]"
          << endl;
      continue;
    }

    const int unit_len = repeat_spec.units[0].length();
    double kPropCorrectMolecules = 0.97;
    if (unit_len <= 2) {
      kPropCorrectMolecules = 0.70;
    }
    const double hap_depth = parameters.depth() / 2;
    const int max_num_units_in_read =
        (int)(std::ceil(parameters.read_len() / (double)unit_len));

    GenotypeRepeat(parameters, repeat_spec, max_num_units_in_read,
                   kPropCorrectMolecules, hap_depth, parameters.read_len(),
                   haplotype_candidates, flanking_size_counts,
                   spanning_size_counts, genotype_type,
                   region_findings.genotype);

    // End genotyping
    sample_findings.push_back(region_findings);

    TestFlankingReads(parameters, region_findings.flanking_repaligns);

    OutputRepeatAligns(parameters, repeat_spec, region_findings.read_groups,
                       region_findings.flanking_repaligns, &(*outputs).log());
  }

  std::sort(sample_findings.begin(), sample_findings.end(),
            CompareRegionFindings);
  WriteJson(parameters, repeat_specs, sample_findings, outputs->json());
  WriteVcf(parameters, repeat_specs, sample_findings, outputs->vcf());
  cerr << TimeStamp() << ",[All done]" << endl;
}

void ReorientReads(const GraphSharedPtr &graph_ptr, int32_t kmer_len,
                   vector<reads::ReadPtr> &read_ptrs) {
  StrandClassifier strand_classifier(graph_ptr, kmer_len);
  for (reads::ReadPtr &read_ptr : read_ptrs) {
    reads::ReorientRead(strand_classifier, *read_ptr);
  }
}

void AlignReadsToGraph(const GraphSharedPtr &graph_ptr, int32_t kmer_len,
                       vector<reads::ReadPtr> &read_ptrs) {
  GaplessAligner aligner(graph_ptr, kmer_len);
  StrMappingClassifier mapping_classifier(0, 1, 2);
  int32_t str_unit_len = graph_ptr->NodeSeq(1).length();
  StrOverlapQuantifier str_overlap_quantifier(0, 1, 2, str_unit_len);

  int32_t num_reads_aligned = 0;
  int32_t num_reads_passed_filter = 0;
  for (reads::ReadPtr &read_ptr : read_ptrs) {
    list<GraphMapping> mappings = aligner.GetBestAlignment(read_ptr->Bases());
    ++num_reads_aligned;

    if (mappings.empty()) {
      continue;
    }

    const GraphMapping canonical_mapping =
        mapping_classifier.GetCanonicalMapping(mappings);

    const int32_t num_matches = canonical_mapping.NumMatches();
    const int32_t reference_span = canonical_mapping.ReferenceSpan();
    const double prop_matches =
        static_cast<double>(num_matches) / reference_span;
    if (prop_matches <= 0.8) {
      continue;
    }
    ++num_reads_passed_filter;

    read_ptr->SetCanonicalMapping(canonical_mapping);

    const MappingType mapping_type =
        mapping_classifier.Classify(canonical_mapping);
    read_ptr->SetCanonicalMappingType(mapping_type);

    const int32_t num_str_units_spanned =
        str_overlap_quantifier.NumUnitsOverlapped(canonical_mapping);
    read_ptr->SetNumStrUnitsSpanned(num_str_units_spanned);
  }

  spd::get("console")->info("{} out of {} reads passed filter",
                            num_reads_passed_filter, num_reads_aligned);
}

void OutputGraphAlignments(const RepeatSpec &repeat_spec,
                           const vector<reads::ReadPtr> &read_ptrs,
                           std::ostream &out) {
  out << repeat_spec.repeat_id << ":" << endl;
  const int32_t indentation_size = 2;
  const string spacer(indentation_size, ' ');
  for (const reads::ReadPtr &read_ptr : read_ptrs) {
    if (read_ptr->CanonicalMappingType() == MappingType::kUnmapped) {
      continue;
    }
    out << spacer << "- name: " << read_ptr->FragmentId() << endl;
    out << spacer << spacer << "type: " << read_ptr->CanonicalMappingType()
        << endl;
    out << spacer << spacer
        << "graph_cigar: " << read_ptr->CanonicalMapping().GetCigarString()
        << endl;
    out << spacer << spacer << "alignment: |" << endl;
    out << EncodeGraphMapping(read_ptr->CanonicalMapping(),
                              3 * indentation_size)
        << endl;
  }
}

int32_t ComputeFlankingHaplotypeCandidate(
    const map<int32_t, int32_t> &flanking_size_counts,
    const map<int32_t, int32_t> &spanning_size_counts) {
  int32_t longest_spanning = 0;
  for (const auto &size_count : spanning_size_counts) {
    if (size_count.second > longest_spanning) {
      longest_spanning = size_count.second;
    }
  }

  int32_t longest_flanking = 0;
  for (const auto &size_count : flanking_size_counts) {
    if (size_count.second > longest_flanking &&
        size_count.second > longest_spanning) {
      longest_flanking = size_count.second;
    }
  }

  return longest_flanking;
}

vector<RepeatAllele> GenerateHaplotypeCandidates(
    const vector<reads::ReadPtr> &read_ptrs,
    map<int32_t, int32_t> &flanking_size_counts,
    map<int32_t, int32_t> &spanning_size_counts) {
  SummarizeAlignments(read_ptrs, flanking_size_counts, spanning_size_counts);

  vector<RepeatAllele> haplotype_candidates;

  for (const auto &size_count : spanning_size_counts) {
    haplotype_candidates.push_back(
        RepeatAllele(size_count.first, size_count.second, ReadType::kSpanning));
  }

  const int32_t flanking_haplotype_candidate =
      ComputeFlankingHaplotypeCandidate(flanking_size_counts,
                                        spanning_size_counts);
  if (flanking_haplotype_candidate != 0) {
    haplotype_candidates.push_back(
        RepeatAllele(flanking_haplotype_candidate, 0, ReadType::kFlanking));
  }

  return haplotype_candidates;
}

void RunGenotyper(const Parameters &parameters, const RepeatSpec &repeat_spec,
                  const vector<RepeatAllele> &haplotype_candidates,
                  const map<int32_t, int32_t> &flanking_size_counts,
                  const map<int32_t, int32_t> &spanning_size_counts,
                  GenotypeType &genotype_type, RepeatGenotype &genotype) {
  const string &repeat_unit = repeat_spec.units_shifts[0][0];
  const int32_t unit_len = repeat_unit.length();
  const double prop_correct_molecules = unit_len > 2 ? 0.97 : 0.70;
  const double hap_depth = parameters.depth() / 2;
  const int max_num_units_in_read =
      (int)(std::ceil(parameters.read_len() / (double)unit_len));

  spd::get("console")->warn("max_num_units_in_read: {}", max_num_units_in_read);
  spd::get("console")->warn("prop_correct_molecules: {}",
                            prop_correct_molecules);
  spd::get("console")->warn("hap_depth: {}", hap_depth);
  spd::get("console")->warn("parameters.read_len(): {}", parameters.read_len());

  GenotypeShortRepeat(max_num_units_in_read, prop_correct_molecules, hap_depth,
                      parameters.read_len(), haplotype_candidates,
                      flanking_size_counts, spanning_size_counts, genotype_type,
                      genotype);
}

string EncodeCounts(const map<int32_t, int32_t> &size_counts) {
  string encoding;
  for (const auto &size_count : size_counts) {
    encoding += "(" + std::to_string(size_count.first) + ", " +
                std::to_string(size_count.second) + ") ";
  }
  return encoding;
}

int main(int argc, char *argv[]) {
  auto console = spd::stderr_color_mt("console");
  spd::set_pattern("%Y-%m-%dT%H:%M:%S,[%v]");

  try {
    Parameters parameters;
    console->info("Starting {}", kProgramVersion);

    if (!parameters.Load(argc, argv)) {
      return 1;
    }

    console->info("Analyzing sample {}", parameters.sample_name());

    Outputs outputs(parameters.vcf_path(), parameters.json_path(),
                    parameters.log_path());

    map<string, RepeatSpec> repeat_specs;
    LoadRepeatSpecs(parameters.repeat_specs_path(), parameters.genome_path(),
                    parameters.min_wp(), repeat_specs);

    const int read_len = CalcReadLen(parameters.bam_path());
    parameters.set_read_len(read_len);

    reads::AlignedReader aligned_reader(parameters.bam_path(),
                                        parameters.genome_path());

    reads::ReadPairs read_pairs;
    for (const auto &repeat_id_spec : repeat_specs) {
      const RepeatSpec &repeat_spec = repeat_id_spec.second;
      GenotypeType genotype_type = GenotypeType::kDiploid;
      const string chrom = repeat_spec.target_region.chrom();

      const bool is_female_chrom_y =
          parameters.sex() == Sex::kFemale && (chrom == "chrY" || chrom == "Y");

      if (is_female_chrom_y) {
        continue;
      }

      const bool is_sex_chrom =
          chrom == "chrX" || chrom == "X" || chrom == "chrY" || chrom == "Y";

      if (parameters.sex() == Sex::kMale && is_sex_chrom) {
        genotype_type = GenotypeType::kHaploid;
      }

      ExtractReads(repeat_spec, parameters.region_extension_len(),
                   aligned_reader, read_pairs);
      console->info("Extracted {} reads from {}", read_pairs.NumReads(),
                    repeat_spec.repeat_id);
      vector<reads::ReadPtr> read_ptrs;
      read_pairs.GetReads(read_ptrs);

      console->info("Constructing STR graph");
      GraphSharedPtr graph_ptr =
          MakeStrGraph(repeat_spec.left_flank, repeat_spec.units_shifts[0][0],
                       repeat_spec.right_flank);

      console->info("Reorienting reads");
      const int32_t kmer_len = 14;
      ReorientReads(graph_ptr, kmer_len, read_ptrs);
      console->info("Aligning reads to graph");
      AlignReadsToGraph(graph_ptr, kmer_len, read_ptrs);
      console->info("Writing alignments to log file");
      OutputGraphAlignments(repeat_spec, read_ptrs, outputs.log());

      console->info("Generating haplotype candidates");
      map<int32_t, int32_t> flanking_size_counts;
      map<int32_t, int32_t> spanning_size_counts;

      vector<RepeatAllele> haplotype_candidates = GenerateHaplotypeCandidates(
          read_ptrs, flanking_size_counts, spanning_size_counts);

      console->info("Flanking counts: {}", EncodeCounts(flanking_size_counts));
      console->info("Spanning counts: {}", EncodeCounts(spanning_size_counts));

      string candidates_encoding;
      for (const RepeatAllele &repeat_allele : haplotype_candidates) {
        candidates_encoding += std::to_string(repeat_allele.size_) + " ";
      }
      console->info("Haplotype candidates: {}", candidates_encoding);

      console->info("Genotyping the repeat");
      RepeatGenotype genotype;
      RunGenotyper(parameters, repeat_spec, haplotype_candidates,
                   flanking_size_counts, spanning_size_counts, genotype_type,
                   genotype);

      string genotype_encoding;
      for (const RepeatAllele &allele : genotype) {
        if (!genotype_encoding.empty()) {
          genotype_encoding += "/";
        }
        genotype_encoding += std::to_string(allele.size_);
      }

      console->info("{}\t{}\t{}", parameters.sample_name(),
                    repeat_spec.repeat_id, genotype_encoding);
    }

    /*
      if (!parameters.depth_is_set()) {
        cerr << TimeStamp() << ",[Calculating depth]" << endl;
        const double depth =
            bam_file.CalcMedianDepth(parameters, parameters.read_len());
        parameters.set_depth(depth);
      }

      cerr << TimeStamp() << ",[Read length: " << parameters.read_len() << "]"
           << endl;
      cerr << TimeStamp() << ",[Depth: " << parameters.depth() << "]" << endl;

      if (parameters.depth() < parameters.kSmallestPossibleDepth) {
        throw std::runtime_error("Estimated depth of " +
                                 lexical_cast<string>(parameters.depth()) +
                                 " is too low for a meaningful inference of "
                                 "repeat sizes");
      }

      EstimateRepeatSizes(parameters, repeat_specs, &bam_file, &outputs); */
  } catch (const std::exception &e) {
    console->error(e.what());
    return 1;
  }

  return 0;
}
