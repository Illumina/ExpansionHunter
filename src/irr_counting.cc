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

#include "include/irr_counting.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "common/timestamp.h"
#include "purity/purity.h"
#include "rep_align/rep_align.h"

using std::string;
using std::endl;
using std::cerr;
using std::unordered_set;
using std::vector;
using std::map;
using std::array;

// Check if two alignments are same.
static bool SameAlign(Align &al1, Align &al2) {
  return (al1.name == al2.name && al1.mate_pos == al2.mate_pos &&
          al1.flag == al2.flag && al1.bases == al2.bases);
}

void CacheReadsFromRegion(const Region &region, const WhatToCache whatToCache,
                          const vector<vector<string>> &units_shifts,
                          double min_wp_score, BamFile *bam_file,
                          AlignPairs *align_pairs) {
  // Jump to the target region or the unaligned reads if the chromosome name
  // is "*" then jump to unaligned reads.
  if (region.chrom() == "*") {
    if (!(*bam_file).JumpToUnaligned()) {
      cerr << TimeStamp() << ",\t[Warning: there appears to be no aligned reads]" << endl;
      return;
    }
  } else {
    if (!(*bam_file).SetRegionToRange(region)) {
      throw std::runtime_error("Failed to jump to " + region.ToString());
    }
  }

  Align align;

  while ((*bam_file).GetRead(align)) {
    const AlignPairs::iterator it = (*align_pairs).find(align.name);

    if (it == (*align_pairs).end()) {
      Align dummy_align;

      align.region = region.ToString();
      array<Align, 2> pair = {align, dummy_align};
      if (!align.IsFirstMate()) {
        std::reverse(pair.begin(), pair.end());
      }
      (*align_pairs)[align.name] = pair;
    } else {
      // The same alignment might be encountered twice if two confusion
      // regions are near each other. Such duplicate alignments are skipped
      // with a warning. Two alignments for the same mate, however, are not
      // permitted.

      align.region = region.ToString();
      AlignPair &frag = (*align_pairs)[align.name];

      if (align.IsFirstMate()) {
        if (frag[0].name.empty()) {
          frag[0] = align;
        } else {
          if (!SameAlign(frag[0], align)) {
            cerr << TimeStamp()
                 << ",\t[WARNING: There are multiple first mates named \""
                 << frag[0].name << "\"]" << endl;
          }
        }
      } else {
        if (frag[1].name.empty()) {
          frag[1] = align;
        } else {
          if (!SameAlign(frag[1], align)) {
            cerr << TimeStamp()
                 << ",\t[WARNING: There are multiple second mates named \""
                 << frag[1].name << "\"]" << endl;
          }
        }
      }

      if (whatToCache == kCacheIrr) {
        // If both mates are cached. Remove both mates from cache
        // unless one of them is in-repeat.
        if (!frag[0].name.empty() && !frag[1].name.empty()) {
          double score1 =
              MatchRepeatRc(units_shifts, frag[0].bases, frag[0].quals);
          double score2 =
              MatchRepeatRc(units_shifts, frag[1].bases, frag[1].quals);
          score1 /= frag[0].bases.length();
          score2 /= frag[1].bases.length();

          if (score1 < min_wp_score && score2 < min_wp_score) {
            (*align_pairs).erase(it);
          }
        }
      }
    }
  }
}

bool CheckAnchoredIrrs(const BamFile &bam_file, const Parameters &parameters,
                       const Region &target_neighborhood,
                       const RepeatSpec &repeat_spec, const Align &read_align,
                       const Align &mate_align,
                       const vector<vector<string>> &units_shifts) {
  const int min_mapq = parameters.min_anchor_mapq();
  // Check if the read has low mapping quality and it is an off-target
  // anchor; such reads are reported but not included into the calculation.
  if (read_align.mapq < min_mapq) {
    if (mate_align.IsMapped() && mate_align.mapq >= min_mapq) {
      Region mateRegion;
      mate_align.GetReadRegion(mateRegion, bam_file.ref_vec());

      if (!mateRegion.Overlaps(target_neighborhood)) {
        cerr << TimeStamp() << ",\t[Discarding IRR " << read_align.name
             << " read" << (read_align.IsFirstMate() ? "1" : "2")
             << read_align.pos << " MAPQ " << read_align.mapq
             << ") because anchoring mate " << mate_align.name << " read"
             << (mate_align.IsFirstMate() ? "1" : "2") << mate_align.pos
             << " MAPQ " << mate_align.mapq << ") not on target ("
             << target_neighborhood << ")]" << endl;
        return false;
      }
    }
  }

  double wp = MatchRepeatRc(units_shifts, read_align.bases, read_align.quals);
  if (!read_align.bases.empty()) {
    wp /= read_align.bases.length();
  }
  const bool is_irr = wp >= parameters.min_wp();

  // no repeat above the length/score threshold was detected in the read
  if (!is_irr) {
    return false;
  }

  if (mate_align.mapq < min_mapq) {
    return false;
  }

  return true;
}

/*****************************************************************************/

void FillinMates(BamFile &bam_file, AlignPairs &align_pairs,
                 const vector<vector<string>> &units_shifts,
                 double min_wp_score,
                 const unordered_set<string> &ontarget_frag_names) {
  for (AlignPairs::iterator it = align_pairs.begin(); it != align_pairs.end();
       ++it) {
    AlignPair &frag = it->second;

    // At least one read must always be filled out.
    assert(!frag[0].name.empty() || !frag[1].name.empty());
    // Try to fill in the missing mate.
    if (frag[0].name.empty() || frag[1].name.empty()) {
      // Get references for exisitng Align and the one that
      // needs to be filled in and then process them below. This
      // is done to avoid code duplication.
      Align &existing_al = frag[0].name.empty() ? frag[1] : frag[0];
      Align &missing_al = frag[0].name.empty() ? frag[0] : frag[1];

      // Do not recover nearby mates.
      if ((existing_al.chrom_id == existing_al.mate_chrom_id) &&
          (std::abs(existing_al.pos - existing_al.mate_pos) < 1000)) {
        continue;
      }

      // Do not recover mates of off-target reads that are not IRRs.
      if (ontarget_frag_names.find(existing_al.name) ==
          ontarget_frag_names.end()) {
        double score =
            MatchRepeatRc(units_shifts, existing_al.bases, existing_al.quals);
        score /= existing_al.bases.length();

        if (score < min_wp_score) {
          continue;
        }
      }

      if (bam_file.GetAlignedMate(existing_al, missing_al)) {
        // region typically stores the position of the region from
        // which a read was cached. Since this read was not cached
        // from anywhere, we store its position in the region field.
        const vector<string> &refVec = bam_file.ref_vec();
        assert(missing_al.chrom_id < refVec.size());
        missing_al.region = Region(refVec[missing_al.chrom_id],
                                   missing_al.pos + 1, missing_al.pos + 2)
                                .ToString();
      } else {
        missing_al = Align();
        missing_al.region = Region("chr-1", 0, 0).ToString();
      }
    }
  }
}

// TODO(edolzhenko): Move cache creation out of the function.
bool CountUnalignedIrrs(BamFile &bam_file, const Parameters &parameters,
                        int &numUnalignedIRRs,
                        const vector<vector<string>> &units_shifts,
                        vector<RepeatAlign> *irr_rep_aligns) {
  numUnalignedIRRs = 0;

  AlignPairs align_pairs;
  Region unalignedRegion("*", 0, 0, "");

  cerr << TimeStamp() << ",\t[Caching unaligned IRRs]" << endl;
  CacheReadsFromRegion(unalignedRegion, kCacheIrr, units_shifts,
                       parameters.min_wp(), &bam_file, &align_pairs);

  cerr << TimeStamp() << ",\t[Done; cached " << align_pairs.size()
       << " unaligned fragments containing at least one IRR read]" << endl;

  for (AlignPairs::const_iterator it = align_pairs.begin();
       it != align_pairs.end(); ++it) {
    const AlignPair &frag = it->second;

    double topScore1 =
        MatchRepeatRc(units_shifts, frag[0].bases, frag[0].quals);
    double topScore2 =
        MatchRepeatRc(units_shifts, frag[1].bases, frag[1].quals);
    if (!frag[0].bases.empty()) {
      topScore1 /= frag[0].bases.length();
    }
    if (!frag[1].bases.empty()) {
      topScore2 /= frag[1].bases.length();
    }

    const bool is_irr1 = topScore1 >= parameters.min_wp();
    const bool is_irr2 = topScore2 >= parameters.min_wp();

    assert(!(frag[0].name.empty() && is_irr1));
    assert(!(frag[1].name.empty() && is_irr2));

    if (is_irr1 && is_irr2) {
      numUnalignedIRRs += 2;

      RepeatAlign rep_align;
      rep_align.read.name = it->first;
      rep_align.read.bases = frag[0].bases;
      rep_align.read.quals = frag[0].quals;
      rep_align.left_flank_len = 0;
      rep_align.right_flank_len = 0;
      rep_align.type = RepeatAlign::Type::kUnalignedIrrPair;
      const int unit_len = units_shifts[0][0].length();
      rep_align.size = rep_align.read.bases.length() / unit_len;
      rep_align.mate.bases = frag[1].bases;
      rep_align.mate.quals = frag[1].quals;
      irr_rep_aligns->push_back(rep_align);
    } else if (is_irr1 || is_irr2) {
      ++numUnalignedIRRs;
      const Align &irr = is_irr1 ? frag[0] : frag[1];
      const Align &mate = is_irr1 ? frag[1] : frag[0];

      RepeatAlign rep_align;
      rep_align.read.name = it->first;
      rep_align.read.bases = irr.bases;
      rep_align.read.quals = irr.quals;
      rep_align.left_flank_len = 0;
      rep_align.right_flank_len = 0;
      rep_align.type = RepeatAlign::Type::kUnalignedIrrSingleton;
      const int unit_len = units_shifts[0][0].length();
      rep_align.size = rep_align.read.bases.length() / unit_len;
      rep_align.mate.bases = mate.bases;
      rep_align.mate.quals = mate.quals;
      irr_rep_aligns->push_back(rep_align);
    }
  }
}

int CountAlignedIrr(const BamFile &bam_file, const Parameters &parameters,
                    const AlignPairs &align_pairs,
                    map<string, int> &numIrrConfRegion,
                    const vector<vector<string>> &units_shifts,
                    vector<RepeatAlign> *irr_rep_aligns) {
  int irr_count = 0;
  bool isFwdKmer = false;
  for (AlignPairs::const_iterator it = align_pairs.begin();
       it != align_pairs.end(); ++it) {
    const AlignPair &frag = it->second;
    double first_top_score =
        MatchRepeatRc(units_shifts, frag[0].bases, frag[0].quals);
    if (!frag[0].bases.empty()) {
      first_top_score /= frag[0].bases.length();
    }
    const bool is_first_irr = first_top_score >= parameters.min_wp();

    double second_top_score =
        MatchRepeatRc(units_shifts, frag[1].bases, frag[1].quals);
    if (!frag[1].bases.empty()) {
      second_top_score /= frag[1].bases.length();
    }
    const bool is_second_irr = second_top_score >= parameters.min_wp();

    // Update status and increase count if both mates are IRR.
    if (is_first_irr && is_second_irr) {
      assert(!frag[0].name.empty() && !frag[1].name.empty());
      irr_count += 2;
      // Increase the count for the corresponding confusion region.
      numIrrConfRegion[frag[0].region] += 1;
      numIrrConfRegion[frag[1].region] += 1;

      RepeatAlign rep_align;
      rep_align.read.name = it->first;
      rep_align.read.bases = frag[0].bases;
      rep_align.read.quals = frag[0].quals;
      rep_align.left_flank_len = 0;
      rep_align.right_flank_len = 0;
      rep_align.type = RepeatAlign::Type::kAlignedIrrPair;
      const int unit_len = units_shifts[0][0].length();
      rep_align.size = rep_align.read.bases.length() / unit_len;
      rep_align.mate.bases = frag[1].bases;
      rep_align.mate.quals = frag[1].quals;
      irr_rep_aligns->push_back(rep_align);
    }
  }

  return irr_count;
}

void CountAnchoredIrrs(const BamFile &bam_file, const Parameters &parameters,
                       const Region &targetNhood, const RepeatSpec &repeat_spec,
                       const unordered_set<string> &ontarget_frag_names,
                       AlignPairs &align_pairs, int &numAnchoredIrrs,
                       const vector<vector<string>> &units_shifts,
                       vector<RepeatAlign> *anchored_irrs) {
  numAnchoredIrrs = 0;

  // Check fragments from the target locus.
  for (unordered_set<string>::const_iterator it = ontarget_frag_names.begin();
       it != ontarget_frag_names.end(); ++it) {
    const AlignPair frag = align_pairs[*it];
    const Align &al1 = frag[0];
    const Align &al2 = frag[1];

    // Counts fragments for which both mates are present.
    if (!al1.name.empty() && !al2.name.empty()) {
      assert(al1.name == al2.name);

      const bool is_mate1_anchored_irr =
          (al1.status != kFlankingRead) &&
          CheckAnchoredIrrs(bam_file, parameters, targetNhood, repeat_spec, al1,
                            al2, units_shifts);

      const bool is_mate2_anchored_irr =
          (al2.status != kFlankingRead) &&
          CheckAnchoredIrrs(bam_file, parameters, targetNhood, repeat_spec, al2,
                            al1, units_shifts);

      if (is_mate1_anchored_irr || is_mate2_anchored_irr) {
        const Align &irr = is_mate1_anchored_irr ? al1 : al2;
        const Align &anchor = is_mate1_anchored_irr ? al2 : al1;

        ++numAnchoredIrrs;
        RepeatAlign rep_align;
        rep_align.read.name = irr.name;
        rep_align.read.bases = irr.bases;
        rep_align.read.quals = irr.quals;
        rep_align.left_flank_len = 0;
        rep_align.right_flank_len = 0;
        rep_align.type = RepeatAlign::Type::kAnchored;
        const int unit_len = units_shifts[0][0].length();
        rep_align.size = rep_align.read.bases.length() / unit_len;
        rep_align.mate.bases = anchor.bases;
        rep_align.mate.quals = anchor.quals;
        anchored_irrs->push_back(rep_align);
      }
    }
  }
}
