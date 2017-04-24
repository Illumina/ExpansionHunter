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
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/parameters.h"
#include "common/ref_genome.h"
#include "genotyping/genotyping.h"
#include "include/bam_file.h"
#include "include/bam_index.h"
#include "include/irr_counting.h"
#include "include/json_output.h"
#include "include/read_group.h"
#include "include/region_findings.h"
#include "include/vcf_output.h"
#include "include/version.h"
#include "purity/purity.h"
#include "rep_align/rep_align.h"

using std::unordered_set;
using std::map;
using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::pair;
using std::array;

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
                      vector<RepeatReadGroup> &read_groups,
                      vector<RepeatAlign> *flanking_repaligns) {
  const int unit_len = repeat_spec.units_shifts[0][0].length();

  map<int, vector<RepeatAlign>> size_spanning_repaligns;
  flanking_repaligns->clear();

  // Align each read to the repeat.
  for (auto &kv : align_pairs) {
    AlignPair &frag = kv.second;

    for (Align &align : frag) {
      RepeatAlign rep_align;
      const bool aligns = AlignRead(parameters, repeat_spec, align.bases,
                                    align.quals, &rep_align);
      rep_align.read.name = align.name;

      if (!aligns) {
        continue;
      }

      // Not pretty, but lets downstream code know that this is not an IRR.
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
    read_group.supported_by = RepeatReadGroup::SupportType::kSpanning;
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

bool is_flannking_allele(const RepeatReadGroup &repeat) {
  return repeat.supported_by == RepeatReadGroup::SupportType::kFlanking;
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

  cerr << "\t[Found " << region_findings.num_anchored_irrs << " anchored IRRs]"
       << endl
       << "\t[Cached " << align_pairs.size() << " reads]" << endl;

  // Stores the total IRR count.
  region_findings.num_irrs = region_findings.num_anchored_irrs;

  // Look for IRR pairs only if anchored IRRs are found.
  if (region_findings.num_anchored_irrs) {
    if (!repeat_spec.is_common_unit()) {
      cerr << "\t[Counting aligned IRR pairs]" << endl;
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
        cerr << "\t[Counting unaligned IRRs]" << endl;
        CountUnalignedIrrs(bam_file, parameters,
                           region_findings.num_unaligned_irrs,
                           repeat_spec.units_shifts, &irr_rep_aligns);
        region_findings.num_irrs += region_findings.num_unaligned_irrs;
      } else {
        cerr << "\t[Skipping unaligned IRRs]" << endl;
      }
    }

    const size_t unit_len = repeat_spec.units[0].length();

    RepeatReadGroup read_group;
    read_group.supported_by = RepeatReadGroup::SupportType::kInrepeat;

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
                     region_findings.read_groups.end(), is_flannking_allele);

    if (flanking_allele_it != region_findings.read_groups.end()) {
      region_findings.flanking_repaligns.insert(
          region_findings.flanking_repaligns.end(),
          flanking_allele_it->rep_aligns.begin(),
          flanking_allele_it->rep_aligns.end());
    }

    region_findings.read_groups.erase(
        std::remove_if(region_findings.read_groups.begin(),
                       region_findings.read_groups.end(), is_flannking_allele),
        region_findings.read_groups.end());
  }
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
      cerr << "[Skipping " << repeat_header << " because the sample is female]"
           << endl;
      continue;
    }

    cerr << "[Analyzing " << repeat_header << "]" << endl;

    const bool is_sex_chrom =
        chrom == "chrX" || chrom == "X" || chrom == "chrY" || chrom == "Y";

    if (parameters.sex() == Sex::kMale && is_sex_chrom) {
      genotype_type = GenotypeType::kHaploid;
    }

    RegionFindings region_findings;
    region_findings.region_id = repeat_spec.repeat_id;
    region_findings.offtarget_irr_counts =
        vector<int>(repeat_spec.offtarget_regions.size(), 0);

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
    FindShortRepeats(parameters, *bam_file, repeat_spec, align_pairs,
                     region_findings.read_groups,
                     &region_findings.flanking_repaligns);

    cerr << "\t[Estimating long repeat sizes]" << endl;
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
    // int long_repeat_size = -1;
    int num_supporting_irr_reads = -1;

    vector<int> haplotype_candidates;
    // Add count of in-repeat reads to flanking.
    for (const auto &repeat : region_findings.read_groups) {
      cerr << repeat.readtypeToStr.at(repeat.supported_by) << endl;
      if (repeat.supported_by == RepeatReadGroup::SupportType::kSpanning) {
        spanning_size_counts[repeat.size] += repeat.num_supporting_reads;
        haplotype_candidates.push_back(repeat.size);
      } else if (repeat.supported_by ==
                 RepeatReadGroup::SupportType::kInrepeat) {
        flanking_size_counts[num_units_in_read] += repeat.num_supporting_reads;
        haplotype_candidates.push_back(num_units_in_read);
      } else if (repeat.supported_by ==
                 RepeatReadGroup::SupportType::kFlanking) {
        haplotype_candidates.push_back(repeat.size);
        for (const auto &align : repeat.rep_aligns) {
          flanking_size_counts[align.size] += 1;
        }
      } else {
        throw std::logic_error("Do not know how to deal with " +
                               repeat.readtypeToStr.at(repeat.supported_by) +
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

    GenotypeRepeat(max_num_units_in_read, kPropCorrectMolecules, hap_depth,
                   parameters.read_len(), haplotype_candidates,
                   flanking_size_counts, spanning_size_counts, genotype_type,
                   region_findings.genotype, region_findings.genotype_ci,
                   region_findings.genotype_support);

    //
    // vector<int> allele_sizes = region_findings.genotype.ExtractAlleleSizes();
    // for (int i = 0; i < allele_sizes.size(); ++i) {
    //  if (allele_sizes[i] == num_units_in_read) {
    // assert(long_repeat_size != -1);
    // assert(num_supporting_irr_reads != -1);
    // region_findings.genotype[i] = long_repeat_size;
    // region_findings.genotype_support[i].set_num_inrepeat(
    //    num_supporting_irr_reads);
    //}
    //}
    //

    // End genotyping
    sample_findings.push_back(region_findings);

    OutputRepeatAligns(parameters, repeat_spec, region_findings.read_groups,
                       region_findings.flanking_repaligns, &(*outputs).log());
  }

  WriteJson(parameters, repeat_specs, sample_findings, outputs->json());
  WriteVcf(parameters, repeat_specs, sample_findings, outputs->vcf());
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

    const int read_len = CalcReadLen(parameters.bam_path());
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
