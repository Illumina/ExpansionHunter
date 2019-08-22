//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include "graphalign/DagAlignerAffine.hh"

#include "gtest/gtest.h"

using std::string;

using namespace graphalign;
using namespace graphalign::dagAligner;

template <typename T> std::string toString(const T& obj)
{
    std::stringstream ss;
    ss << obj;
    return ss.str();
}

TEST(SimpleAlignment, Short_to_short)
{
    DagAligner<true> aligner({ 5, -4 }, 0, -8);

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 151, 151 } }), std::vector<int>({ 0 }));

    const string query = "tgCccgcCCcCCCCcccC";

    const string reference = "TGCAGTCCCGCCCCGTCCC";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score secondBestScore = 0;
    std::vector<Cigar> cigars;
    const Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);

    EXPECT_EQ(std::size_t(3), cigars.size());
    EXPECT_EQ("0[3=3X3=1X4=1X1D3=]", toString(cigars.at(2)));
    EXPECT_EQ("0[3=3X3=1X4=1D1X3=]", toString(cigars.at(1)));
    EXPECT_EQ("0[3=3X3=1D4=2X3=]", toString(cigars.at(0)));
    EXPECT_EQ(37, bestScore);
}

TEST(SimpleAlignment, Long_to_Long)
{
    DagAligner<true> aligner({ 5, -4 }, 0, -8);

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 151, 151 } }), std::vector<int>({ 0 }));

    const string query = "TCTCGCCCCGCCCCTCAGGCGGCCTCCCTGCtgtgCCCCGCCCCGGCCcCGCCCCgCCCCcCCCCCcCCaCgCCCCCCcCccCcCCCCgCCCC"
                         "CCCCctCcCCCCccctCCCCtccCCtgCccgcCCcCCCCcccC";
    const string reference = "CCGCCCCGCCCCCGTCTCGCCCCGCCCCTCAGGCGGCCTCCCTGCTGTGCCCCGCCCCGGCCTCGCCACGCCCCTACCTCACCACGCCC"
                             "CCCGCATCGCCACGCCCCCCGCATCGCCACGCCTCCCTTACCATGCAGTCCCGCCCCGTCCC";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score secondBestScore = 0;
    std::vector<Cigar> cigars;
    const Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(3), cigars.size());

    EXPECT_EQ(
        "0[14D48=1X4=1X6=2X2=1X1=1X11=1X1=2X1=1X2=1X8=1X1=1X2=1X2=1X1=1X6=1X1=1X2=1X3=3X3=1X4=1X1D3=]",
        toString(cigars.at(2)));
    EXPECT_EQ(
        "0[14D48=1X4=1X6=2X2=1X1=1X11=1X1=2X1=1X2=1X8=1X1=1X2=1X2=1X1=1X6=1X1=1X2=1X3=3X3=1X4=1D1X3=]",
        toString(cigars.at(1)));
    EXPECT_EQ(
        "0[14D48=1X4=1X6=2X2=1X1=1X11=1X1=2X1=1X2=1X8=1X1=1X2=1X2=1X1=1X6=1X1=1X2=1X3=3X3=1D4=2X3=]",
        toString(cigars.at(0)));
    EXPECT_EQ(344, bestScore);
}

TEST(SimpleAlignment, AAAC_to_AGC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 3, 3 } }), std::vector<int>({ 0 }));

    const string query = "AAAC";
    const string reference = "AGC";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score secondBestScore = 0;

    std::vector<Cigar> cigars;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(3), cigars.size());
    EXPECT_EQ("0[1S1=1X1=]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=1I1X1=]", toString(cigars.at(1)));
    EXPECT_EQ("0[1=1X1I1=]", toString(cigars.at(2)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1S1=1X1=]", toString(cigar));
}

TEST(SimpleAlignment, ATGC_to_AGC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "ATGC";
    const string reference = "AGC";

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 3, 3 } }), std::vector<int>({ 0 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=1I2=]", toString(cigar));
}

TEST(SimpleAlignment, TAACTTTTGGG_to_TGCTTTTAA)
{
    // query:     TAACTTTTGGG
    //            |x |||||xx
    // reference: TG-CTTTTAA-

    const string query = "TAACTTTTGGG";
    const string reference = "TGCTTTTAA";

    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 9, 9 } }), std::vector<int>({ 0 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    Score secondBestScore = 0;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=1I1X5=3S]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=1X1I5=3S]", toString(cigars.at(1)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);

    EXPECT_EQ("0[1=1I1X5=3S]", toString(cigar));
}

TEST(SimpleAlignment, TCACGGAGA_to_TACGAGAG)
{
    // TCACGGAGA
    // | ||| |||
    // T-ACG-AGAG

    const string query = "TCACGGAGA";
    const string reference = "TACGAGAG";

    DagAligner<true> aligner({ 5, -4 }, 0, -8);

    EdgeMap edges(std::vector<std::pair<int, int>>({ { 8, 8 } }), std::vector<int>({ 0 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    Score secondBestScore = 0;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=1I2=1I4=]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=1I3=1I3=]", toString(cigars.at(1)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=1I2=1I4=]", toString(cigar));
}

TEST(ForkAlignment, AAAC_to_AAC_fork_AAC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAAC";
    const string reference = "AAC"
                             "AAC";

    EdgeMap edges(std::vector<std::pair<int, int>>({ { -1, 3 }, { 6, 6 } }), std::vector<int>({ 0, 1 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    Score secondBestScore = 0;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(6), cigars.size());
    EXPECT_EQ("0[1S3=]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=1I2=]", toString(cigars.at(1)));
    EXPECT_EQ("0[2=1I1=]", toString(cigars.at(2)));
    EXPECT_EQ("1[1S3=]", toString(cigars.at(3)));
    EXPECT_EQ("1[1=1I2=]", toString(cigars.at(4)));
    EXPECT_EQ("1[2=1I1=]", toString(cigars.at(5)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1S3=]", toString(cigar));
}

TEST(ForkAlignment, AAAC_to_AGC_fork_AAC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAAC";
    const string reference = "AGC"
                             "AAC";
    EdgeMap edges(std::vector<std::pair<int, int>>({ { -1, 3 }, { 6, 6 } }), std::vector<int>({ 0, 1 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    Score secondBestScore = 0;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(3), cigars.size());
    EXPECT_EQ("1[1S3=]", toString(cigars.at(0)));
    EXPECT_EQ("1[1=1I2=]", toString(cigars.at(1)));
    EXPECT_EQ("1[2=1I1=]", toString(cigars.at(2)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("1[1S3=]", toString(cigar));
}

TEST(Fork, AAC_to_AGC_AAC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAC";
    const string reference = "AGC"
                             "AAC";

    EdgeMap edges(std::vector<std::pair<int, int>>({ { -1, 3 }, { 6, 6 } }), std::vector<int>({ 0, 1 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("1[3=]", toString(cigar));
}

/*
 *    A
 *   / \
 *  A   C
 *   \ /
 *    G
 */
TEST(ForkJoin1Base, AAC_to_AGC_AAC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAC";
    const string reference = "A"
                             "A"
                             "G"
                             "C";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 }, { 0, 2 }, { 1, 3 }, { 2, 3 }, { 4, 4 } }),
        std::vector<int>({ 0, 1, 2, 3 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=]1[1=]3[1=]", toString(cigar));
}

/*
 *    AA
 *   /  \
 *  A    C
 *   \  /
 *    AG
 */
TEST(ForkJoin2Base, AAGC_to_AAGC_AAAC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAGC";
    const string reference = "A"
                             "AA"
                             "AG"
                             "C";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 }, { 0, 3 }, { 2, 5 }, { 4, 5 }, { 6, 6 } }),
        std::vector<int>({ 0, 1, 2, 3 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=]2[2=]3[1=]", toString(cigar));
}

/*
 *    AA
 *   /  \
 *  A    C
 *   \  /
 *    AC
 */
TEST(ForkJoin2Base, AAGC_to_AAAC_AACC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAGC";
    const string reference = "A"
                             "AA"
                             "AC"
                             "C";
    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 }, { 0, 3 }, { 2, 5 }, { 4, 5 }, { 6, 6 } }),
        std::vector<int>({ 0, 1, 2, 3 }));
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    Score secondBestScore = 0;
    Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=]1[1=1X]3[1=]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=]2[1=1X]3[1=]", toString(cigars.at(1)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=]2[1=1X]3[1=]", toString(cigar));
}

/*
 *    AA
 *   /  \
 *  A    C
 *   \  /
 *    AC
 */
TEST(ForkJoin2Base, AAC_to_AAAC_AACC)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAC";
    const string reference = "A"
                             "AA"
                             "AC"
                             "C";
    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 }, { 0, 3 }, { 2, 5 }, { 4, 5 }, { 6, 6 } }),
        std::vector<int>({ 0, 1, 2, 3 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[1=]2[2=]", toString(cigar));
}

/*
 *  AAAA
 *      \
 *       T
 *      /
 *  GACC
 */
TEST(JoinStartAtOffset, ACCT_to_AAAAT_GACCT)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "ACCT";
    const string reference = "AAAA"
                             "GACC"
                             "T";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { -1, 4 }, { 3, 8 }, { 7, 8 }, { 9, 9 } }), std::vector<int>({ 0, 1, 2 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("1[1D3=]2[1=]", toString(cigar));
}

/*
 *  AA
 *    \
 *     T
 *    /
 *   C
 */
TEST(JoinStartAtOffset, AAT_to_AAT_CT)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAT";
    const string reference = "AA"
                             "C"
                             "T";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { -1, 2 }, { 1, 3 }, { 2, 3 }, { 4, 4 } }), std::vector<int>({ 0, 1, 2 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    EXPECT_EQ(
        "Aligner(AffineAlignMatrix(AlignMatrix(\n"
        "[0\t-2\t-4\t-6]\n"
        "[-2\t1\t-1\t-3]\n"
        "[-4\t-1\t2\t0]\n"
        "[-2\t-1\t-3\t-5]\n"
        "[-4\t-3\t0\t3]\n"
        ")))",
        toString(aligner));

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[2=]2[1=]", toString(cigar));
}

/*
 *   C
 *    \
 *     T
 *    /
 *  AA
 */
TEST(JoinStartAtOffset, AAT_to_CT_AAT)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAT";
    const string reference = "C"
                             "AA"
                             "T";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { -1, 1 }, { 0, 3 }, { 2, 3 }, { 4, 4 } }), std::vector<int>({ 0, 1, 2 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    EXPECT_EQ(
        "Aligner(AffineAlignMatrix(AlignMatrix(\n"
        "[0\t-2\t-4\t-6]\n"
        "[-2\t-1\t-3\t-5]\n"
        "[-2\t1\t-1\t-3]\n"
        "[-4\t-1\t2\t0]\n"
        "[-4\t-3\t0\t3]\n"
        ")))",
        toString(aligner));

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("1[2=]2[1=]", toString(cigar));
}

/*
 *           AA
 *          /  \
 *  TCGTGTAA    CCCCCCCCTTTTT
 *          \  /
 *           GC
 */
TEST(ForkJoinLong, AAGCCCCCCCCCTTTTT_to_TCGTGTAACCCCCCCCTTTTT_TCGTGTGCCCCCCCCCTTTTT)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAGCCCCCCCCCTTTTT";
    const string reference = "TCGTGTAA"
                             "AA"
                             "GC"
                             "CCCCCCCCTTTTT";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 7, 8 }, { 7, 10 }, { 9, 12 }, { 11, 12 }, { 25, 25 } }),
        std::vector<int>({ 0, 1, 2, 3 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[6D2=]2[2=]3[13=]", toString(cigar));
}

/*
 *           AAA
 *          /   \
 *  TCGTGTAA     CCCCCCCCTTTTT
 *          \   /
 *           GC
 */
TEST(ForkJoinLong, AAGCCCCCCCCCTTTTT_to_TCGTGTAAAAACCCCCCCCTTTTT_TCGTGTAAGCCCCCCCCCTTTTT)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "AAGCCCCCCCCCTTTTT";
    const string reference = "TCGTGTAA"
                             "AAA"
                             "GC"
                             "CCCCCCCCTTTTT";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 7, 8 }, { 7, 11 }, { 10, 13 }, { 12, 13 }, { 26, 26 } }),
        std::vector<int>({ 0, 1, 2, 3 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ("0[6D2=]2[2=]3[13=]", toString(cigar));
}

/*
 *  A---   TT
 *      \ /
 *       T
 *      / \
 *  GACC   C
 */
class SimpleGraphTest : public testing::Test
{
protected:
    string reference;
    EdgeMap edges;
    DagAligner<true> aligner;
    Score bestScore = 0;
    Score secondBestScore = 0;

    SimpleGraphTest()
        : reference("A"
                    "GACC"
                    "T"
                    "TT"
                    "C")
        , edges(
              std::vector<std::pair<int, int>>({ { -1, 1 }, { 0, 5 }, { 4, 5 }, { 5, 6 }, { 5, 8 }, { 9, 9 } }),
              std::vector<int>({ 0, 1, 2, 3, 4 }))
        , aligner({ 1, -1 }, 0, -2)
    {
    }
};

TEST_F(SimpleGraphTest, OffEnd)
{
    string query = "ATCTG";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=]2[1=1I]3[1=1S]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=]2[1=]3[1X1=1S]", toString(cigars.at(1)));

    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(0, bestScore);
    EXPECT_EQ("0[1=]2[1=1I]3[1=1S]", toString(cigar));
}
//
TEST_F(SimpleGraphTest, QueryAllN)
{
    string query = "N";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    bestScore = aligner.backtrackAllPaths<true>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=]", toString(cigars.at(0)));
    EXPECT_EQ("1[1=]", toString(cigars.at(1)));

    Cigar cigar = aligner.backtrackBestPath<true>(edges, bestScore, secondBestScore);
    EXPECT_EQ(1, bestScore);
    EXPECT_EQ("0[1=]", toString(cigar));

    cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(1, bestScore);
    EXPECT_EQ("0[1=]", toString(cigar));

    query = "NNNNNNNNNN";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    cigars.clear();
    bestScore = aligner.backtrackAllPaths<true>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(1), cigars.size());
    EXPECT_EQ("1[4=]2[1=]3[2=3S]", toString(cigars.at(0)));

    cigar = aligner.backtrackBestPath<true>(edges, bestScore, secondBestScore);
    EXPECT_EQ(7, bestScore);
    EXPECT_EQ("1[4=]2[1=]3[2=3S]", toString(cigar));

    cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(1, bestScore);
    EXPECT_EQ("1[3S4=]2[1=]3[2=]", toString(cigar));
}

TEST_F(SimpleGraphTest, QuerySomeN)
{
    string query = "GANCNC";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(6, bestScore);
    EXPECT_EQ("1[4=]2[1=]4[1=]", toString(cigar));
}

TEST_F(SimpleGraphTest, EmptyQuery)
{
    string query = "";
    EXPECT_ANY_THROW(aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges));
}

TEST_F(SimpleGraphTest, BadQualities)
{
    string query = "gACc";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(4, bestScore);
    EXPECT_EQ("1[4=]", toString(cigar));
}

/*
 *     _
 *    / \
 *  G-TCC-AAAAA
 */
TEST(RepeatExpansion, SimpleRepeat)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string reference = "G"
                             "TCC"
                             "TCC"
                             "TCC"
                             "AAAAA";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 },
                                           { 3, 4 },
                                           { 6, 7 },
                                           { 3, 10 },
                                           { 6, 10 },
                                           { 9, 10 },
                                           { reference.length(), reference.length() } }),
        std::vector<int>({ 0, 1, 2, 3, 4 }));

    string query = "TCCTCCAA";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);
    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(6, bestScore);
    EXPECT_EQ("0[1D]1[3=]2[3=]4[2=]", toString(cigar));

    query = "GTCTCCCCAA";
    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    std::vector<Cigar> cigars;
    bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);
    EXPECT_EQ(std::size_t(2), cigars.size());
    EXPECT_EQ("0[1=]1[1=1D1=]2[3=]3[1D2=]4[2=]", toString(cigars.at(0)));
    EXPECT_EQ("0[1=]1[2=1D]2[3=]3[1D2=]4[2=]", toString(cigars.at(1)));

    cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(6, bestScore);
    EXPECT_EQ("0[1=]1[1=1D1=]2[3=]3[1D2=]4[2=]", toString(cigar));
}

/*
 *  |\     |\
 *  G--TCC-C
 */
TEST(RepeatExpansion, HomoPolymers)
{
    DagAligner<true> aligner({ 1, -1 }, 0, -2);

    const string query = "GGTCCGC";
    const string reference = "G"
                             "G"
                             "G"
                             "TCC"
                             "C"
                             "C";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 },
                                           { 1, 2 },
                                           { 0, 3 },
                                           { 1, 3 },
                                           { 2, 3 },
                                           { 5, 6 },
                                           { 6, 7 },
                                           { reference.length(), reference.length() } }),
        std::vector<int>({ 0, 1, 2, 3, 4, 5 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score bestScore = 0;
    Score secondBestScore = 0;
    Cigar cigar = aligner.backtrackBestPath<false>(edges, bestScore, secondBestScore);
    EXPECT_EQ(5, bestScore);
    EXPECT_EQ("0[1=]1[1=]3[3=]4[1X]5[1=]", toString(cigar));
}

/*
 *              ---
 *             /   \
 * C-GCC-GACAAC-GAC-CTTCCTGAACT
 *   \ /        \ /
 *    -          -
 */
TEST(RepeatExpansion, two_repeats)
{
    DagAligner<true> aligner({ 5, -4 }, 0, -8);

    const string query = "CgCCGCCA";
    const string reference = "C"
                             "GCC"
                             "GCC"
                             "GCC"
                             "GCC"
                             "GCC"
                             "GCC"
                             "GACAAC"
                             "GAC"
                             "GAC"
                             "GAC"
                             "GAC"
                             "CTTCCTGAACT";

    EdgeMap edges(
        std::vector<std::pair<int, int>>({ { 0, 1 },
                                           { 3, 4 },
                                           { 6, 7 },
                                           { 9, 10 },
                                           { 12, 13 },
                                           { 15, 16 },
                                           { 0, 19 },
                                           { 3, 19 },
                                           { 6, 19 },
                                           { 9, 19 },
                                           { 12, 19 },
                                           { 15, 19 },
                                           { 18, 19 },
                                           { 24, 25 },
                                           { 27, 28 },
                                           { 30, 31 },
                                           { 33, 34 },
                                           { 27, 37 },
                                           { 30, 37 },
                                           { 33, 37 },
                                           { 36, 37 },
                                           { 24, 37 },
                                           { reference.length(), reference.length() } }),
        std::vector<int>({ 6, 7, 8, 9, 10, 11, 12, 5, 1, 2, 3, 4, 0 }));

    aligner.align(query.begin(), query.end(), reference.begin(), reference.end(), edges);

    Score secondBestScore = 0;
    std::vector<Cigar> cigars;
    const Score bestScore = aligner.backtrackAllPaths<false>(edges, cigars, secondBestScore);

    EXPECT_EQ(32, bestScore);
    EXPECT_EQ("6[1=]7[3=]8[3=]5[1D1=]", toString(cigars.at(0)));
}
