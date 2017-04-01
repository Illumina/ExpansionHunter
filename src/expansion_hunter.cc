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

#include <boost/lexical_cast.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::lexical_cast;
#include <boost/algorithm/string/join.hpp>

#include <cassert>
#include <iostream>
using std::cerr;
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <unordered_set>
using std::unordered_set;
#include <sstream>
#include <stdexcept>
#include <utility>
using std::pair;
#include "genotyping/genotyping.h"

#include "include/allele.h"
#include "include/bam_file.h"
#include "include/bam_index.h"
#include "include/irr_counting.h"
#include "include/parameters.h"
#include "include/ref_genome.h"
#include "include/repeat_length.h"
#include "include/version.h"
#include "purity/purity.h"
#include "rep_align/rep_align.h"

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

// Search for reads spanning the entire repeat sequence.
void FindShortRepeats(const Parameters &parameters, BamFile &bam_file,
                      const RepeatSpec &repeat_spec, AlignPairs &align_pairs,
                      vector<Allele> &alleles,
                      vector<RepeatAlign> *flanking_repaligns) {
  const size_t unit_len = repeat_spec.units_shifts[0][0].length();

  map<size_t, vector<RepeatAlign>> size_spanning_repaligns;
  flanking_repaligns->clear();

  // Align each read to the repeat.
  for (auto &kv : align_pairs) {
    AlignPair &frag = kv.second;

    for (Align &align : frag) {
      RepeatAlign rep_align;
      const bool aligns = AlignRead(
          parameters.min_baseq(), parameters.min_wp(), repeat_spec.units_shifts,
          repeat_spec.left_flank, repeat_spec.right_flank, align.bases,
          align.quals, &rep_align);
      rep_align.name = align.name;

      if (!aligns) {
        continue;
      }

      // Not pretty, but lets downstream code know that this is not an IRR.
      align.status = kFlankingRead;

      if (rep_align.type == kSpanning) {
        size_spanning_repaligns[rep_align.size].push_back(rep_align);
      } else {
        assert(rep_align.type == kFlanking);
        flanking_repaligns->push_back(rep_align);
      }
    }
  }

  for (const auto &size_repaligns : size_spanning_repaligns) {
    Allele allele;
    allele.type = kSpanningAllele;
    allele.size = size_repaligns.first;
    allele.size_ci_lower = size_repaligns.first;
    allele.size_ci_upper = size_repaligns.first;
    allele.num_supporting_reads = size_repaligns.second.size();
    allele.rep_aligns = size_repaligns.second;
    alleles.push_back(allele);
  }

  DistributeFlankingReads(parameters, repeat_spec, &alleles,
                          flanking_repaligns);

  const double haplotype_depth = parameters.depth() / 2;
  CoalesceFlankingReads(repeat_spec, alleles, flanking_repaligns,
                        parameters.read_len(), haplotype_depth,
                        repeat_spec.units[0].length(), repeat_spec.units_shifts,
                        parameters.min_baseq(), parameters.min_wp());
}

// Caches alignments from extended target and off-target regions. For BAM files,
// mates of reads in the relevant regions are cached too.
void CacheAligns(BamFile *bam_file, const RepeatSpec &repeat_spec,
                 AlignPairs &align_pairs,
                 unordered_set<string> &ontarget_frag_names,
                 size_t extension_len) {
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

  cerr << "\t[Found " << ontarget_frag_names.size() << " reads in target locus]"
       << endl;

  // Cache aligned off-target reads.
  for (const Region &confusion_region : repeat_spec.offtarget_regions) {
    Region extended_confusion_region = confusion_region.Extend(extension_len);
    CacheReadsFromRegion(extended_confusion_region, kCacheAll,
                         repeat_spec.units_shifts, 0.9, &(*bam_file),
                         &align_pairs);
  }

  // Filling-in missing mates by jumping around the BAM.
  cerr << "\t[Filling in mates]" << endl;
  if ((*bam_file).format() == BamFile::kBamFile) {
    FillinMates(*bam_file, align_pairs);
  } else {
    cerr << "\t[Skipping filling in mates]" << endl;
  }
  cerr << "\t[Done filling in mates]" << endl;
}

bool is_flannking_allele(const Allele &allele) {
  return allele.type == kFlankingAllele;
}

bool FindLongRepeats(const Parameters &parameters,
                     const RepeatSpec &repeat_spec, BamFile &bam_file,
                     vector<Allele> &alleles, size_t &num_irrs,
                     size_t &num_unaligned_irrs, size_t &num_anchored_irrs,
                     vector<size_t> &off_target_irr_counts,
                     unordered_set<string> &ontarget_frag_names,
                     AlignPairs &align_pairs,
                     vector<RepeatAlign> *flanking_repaligns) {
  // Count the number of anchored IRRs.
  num_anchored_irrs = 0;
  const size_t extension_len = parameters.region_extension_len();
  Region target_nhood = repeat_spec.target_region.Extend(extension_len);

  // Anchored IRRs and IRR pairs will be stored here.
  vector<RepeatAlign> irr_rep_aligns;
  CountAnchoredIrrs(bam_file, parameters, target_nhood, repeat_spec,
                    ontarget_frag_names, align_pairs, num_anchored_irrs,
                    repeat_spec.units_shifts, &irr_rep_aligns);

  cerr << "\t[Found " << num_anchored_irrs << " anchored IRRs]" << endl
       << "\t[Cached " << align_pairs.size() << " reads]" << endl;

  // Stores the total IRR count.
  num_irrs = num_anchored_irrs;

  // Look for IRR pairs only if anchored IRRs are found.
  if (num_anchored_irrs) {
    if (!repeat_spec.is_common_unit()) {
      cerr << "\t[Counting aligned IRR pairs]" << endl;
      map<string, size_t> numIrrConfRegion;
      num_irrs +=
          CountAlignedIrr(bam_file, parameters, align_pairs, numIrrConfRegion,
                          repeat_spec.units_shifts, &irr_rep_aligns);

      // Record paired IRR counts from each confusion region.
      off_target_irr_counts.clear();
      for (const Region &confusionRegion : repeat_spec.offtarget_regions) {
        Region confusionNhood = confusionRegion.Extend(extension_len);
        off_target_irr_counts.push_back(
            numIrrConfRegion[confusionNhood.AsString()]);
      }

      num_unaligned_irrs = 0;

      if (!parameters.skip_unaligned()) {
        cerr << "\t[Counting unaligned IRRs]" << endl;
        CountUnalignedIrrs(bam_file, parameters, num_unaligned_irrs,
                           repeat_spec.units_shifts, &irr_rep_aligns);
        num_irrs += num_unaligned_irrs;
      } else {
        cerr << "\t[Skipping unaligned IRRs]" << endl;
      }
    }

    const size_t unit_len = repeat_spec.units[0].length();

    cerr << "\t[Estimating repeat length from IRRs]" << endl;
    Allele allele;
    allele.type = kInRepeatAllele;

    const double haplotype_depth = parameters.depth() / 2;
    EstimateRepeatLen(num_irrs, parameters.read_len(), haplotype_depth,
                      allele.size, allele.size_ci_lower, allele.size_ci_upper);

    allele.size /= unit_len;
    allele.size_ci_lower /= unit_len;
    allele.size_ci_upper /= unit_len;
    allele.num_supporting_reads = num_irrs;
    allele.rep_aligns = irr_rep_aligns;

    alleles.push_back(allele);

    // If there is evidence for a long repeat allele we assume that flanking
    // reads came from and so any previously-found alleles whose size was
    // estimated from flanking reads are no longer valid.

    vector<Allele>::const_iterator flanking_allele_it =
        std::find_if(alleles.begin(), alleles.end(), is_flannking_allele);

    if (flanking_allele_it != alleles.end()) {
      flanking_repaligns->insert(flanking_repaligns->end(),
                                 flanking_allele_it->rep_aligns.begin(),
                                 flanking_allele_it->rep_aligns.end());
    }

    alleles.erase(
        std::remove_if(alleles.begin(), alleles.end(), is_flannking_allele),
        alleles.end());
  }
}

// For each repeat region, calculate sizes of repeats that are (a) shorter and
// (b) longer than the read length; reconsile alleles and output results to
// appropriate files.

void EstimateRepeatSizes(const Parameters &parameters,
                         const map<string, RepeatSpec> &repeat_specs,
                         BamFile *bam_file, Outputs *outputs) {
  boost::property_tree::ptree ptree_root;

  // Analyze repeats one by one.
  for (auto &kv : repeat_specs) {
    const RepeatSpec &repeat_spec = kv.second;

    string repeat_header = repeat_spec.target_region.AsString();
    repeat_header += " " + boost::algorithm::join(repeat_spec.units, "/");
    cerr << "[Analyzing " << repeat_header << "]" << endl;

    cerr << "\t[Caching reads]" << endl;
    AlignPairs align_pairs;
    unordered_set<string> ontarget_frag_names;
    CacheAligns(&(*bam_file), repeat_spec, align_pairs, ontarget_frag_names,
                parameters.region_extension_len());
    if (align_pairs.empty()) {
      cerr << "\t[Found no on-target or off-target reads]" << endl;
      continue;
    }

    cerr << "\t[Estimating short repeat sizes]" << endl;
    vector<Allele> alleles;
    vector<RepeatAlign> flanking_repaligns;
    FindShortRepeats(parameters, *bam_file, repeat_spec, align_pairs, alleles,
                     &flanking_repaligns);

    cerr << "\t[Estimating long repeat sizes]" << endl;
    vector<size_t> offtarget_irr_counts;
    size_t num_anchored_irrs = 0;
    size_t num_unaligned_irrs = 0;
    size_t num_irrs = 0;

    FindLongRepeats(parameters, repeat_spec, *bam_file, alleles, num_irrs,
                    num_unaligned_irrs, num_anchored_irrs, offtarget_irr_counts,
                    ontarget_frag_names, align_pairs, &flanking_repaligns);

    map<int, int> flanking_size_counts;
    map<int, int> spanning_size_counts;

    for (const auto &align : flanking_repaligns) {
      flanking_size_counts[align.size] += 1;
    }

    vector<int> haplotype_candidates;
    // Add count of in-repeat reads to flanking.
    for (const auto &allele : alleles) {
      cerr << allele.readtypeToStr.at(allele.type) << endl;
      if (allele.type == kSpanningAllele) {
        spanning_size_counts[allele.size] += allele.num_supporting_reads;
        haplotype_candidates.push_back(allele.size);
      } else if (allele.type == kInRepeatAllele) {
        const int num_units_in_read = (int)(std::ceil(
            parameters.read_len() / (double)repeat_spec.units[0].length()));
        const int bounded_num_irrs =
            allele.num_supporting_reads <= 5 ? allele.num_supporting_reads : 5;
        flanking_size_counts[num_units_in_read] += bounded_num_irrs;
        haplotype_candidates.push_back(num_units_in_read);
      } else if (allele.type == kFlankingAllele) {
        haplotype_candidates.push_back(allele.size);
        for (const auto &align : allele.rep_aligns) {
          flanking_size_counts[align.size] += 1;
        }
      } else {
        throw std::logic_error("Do not know how to deal with " +
                               allele.readtypeToStr.at(allele.type) +
                               " alleles");
      }
    }

    cerr << "Flanking:" << endl;
    for (const auto &kv : flanking_size_counts) {
      cerr << "\t" << kv.first << " -- " << kv.second << endl;
    }

    cerr << "Spanning:" << endl;
    for (const auto &kv : spanning_size_counts) {
      cerr << "\t" << kv.first << " -- " << kv.second << endl;
    }

    cerr << "Haplotype candidates: ";
    for (const auto &size : haplotype_candidates) {
      cerr << size << " ";
    }
    cerr << endl;

    const int unit_len = repeat_spec.units[0].length();
    double kPropCorrectMolecules = 0.97;
    if (unit_len <= 2) {
      kPropCorrectMolecules = 0.70;
    }
    const double hap_depth = parameters.depth() / 2;
    const int max_num_units_in_read =
        (int)(std::ceil(parameters.read_len() / (double)unit_len));

    GenotypeType genotype_type = GenotypeType::kDiploid;
    const string chrom = repeat_spec.target_region.chrom();
    if (parameters.sex() == Sex::kMale && (chrom == "chrX" || chrom == "X")) {
      genotype_type = GenotypeType::kHaploid;
    }

    vector<int> genotype = genotypeOneUnitStr(
        max_num_units_in_read, kPropCorrectMolecules, hap_depth,
        parameters.read_len(), haplotype_candidates, flanking_size_counts,
        spanning_size_counts, genotype_type);

    // End genotyping

    boost::property_tree::ptree region_node;
    AsPtree(region_node, alleles, repeat_spec, num_irrs, num_unaligned_irrs,
            num_anchored_irrs, offtarget_irr_counts, genotype);
    ptree_root.put_child(repeat_spec.repeat_id, region_node);

    OutputRepeatAligns(parameters, repeat_spec, alleles, flanking_repaligns,
                       &(*outputs).log());
  }

  boost::property_tree::ptree bam_stats_node;
  bam_stats_node.put<size_t>("ReadLength", parameters.read_len());
  bam_stats_node.put<float>("MedianDepth", parameters.depth());
  ptree_root.put_child("BamStats", bam_stats_node);

  boost::property_tree::json_parser::write_json((*outputs).json(), ptree_root);
  DumpVcf(parameters, repeat_specs, ptree_root, *outputs);
  cerr << "[All done]" << endl;
}

int main(int argc, char *argv[]) {
  try {
    Parameters parameters;
    cerr << kProgramVersion << endl;

    if (!parameters.Load(argc, argv)) {
      return 1;
    }

    Outputs outputs(parameters.vcf_path(), parameters.json_path(),
                    parameters.log_path());

    map<string, RepeatSpec> repeat_specs;
    if (!LoadRepeatSpecs(parameters.repeat_specs_path(),
                         parameters.genome_path(), parameters.min_wp(),
                         &repeat_specs)) {
      throw std::invalid_argument(
          "Failed to load repeat table from disease specs in '" +
          parameters.repeat_specs_path() + "'");
    }

    const size_t read_len = CalcReadLen(parameters.bam_path());
    parameters.set_read_len(read_len);

    BamFile bam_file;
    bam_file.Init(parameters.bam_path(), parameters.genome_path());

    cerr << "[Sample " << parameters.sample_name() << "]" << endl;

    if (!parameters.depth_is_set()) {
      cerr << "[Calculating depth]" << endl;
      const double depth =
          bam_file.CalcMedianDepth(parameters, parameters.read_len());
      parameters.set_depth(depth);
    }

    cerr << "[Read length: " << parameters.read_len() << "]" << endl;
    cerr << "[Depth: " << parameters.depth() << "]" << endl;

    if (parameters.depth() < parameters.kSmallestPossibleDepth) {
      throw std::runtime_error("Estimated depth of " +
                               lexical_cast<string>(parameters.depth()) +
                               " is too low for a meaningful inference of "
                               "repeat sizes");
    }

    EstimateRepeatSizes(parameters, repeat_specs, &bam_file, &outputs);
  } catch (const std::exception &e) {
    cerr << e.what() << endl;
    return 1;
  }

  return 0;
}
