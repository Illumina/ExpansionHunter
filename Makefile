CC=/illumina/thirdparty/gcc/gcc-4.9.2/bin/g++ -O3 -std=c++11 -static

BOOST_BASE=/illumina/thirdparty/boost/boost_1_54_0_python2.7/

BOOST_INC_DIR=${BOOST_BASE}/include
BOOST_LIB_DIR=${BOOST_BASE}/lib
BOOST_LIB_SPEC=-lboost_program_options -lboost_filesystem -lboost_system -lboost_regex -lboost_date_time -lpthread

SAMTOOLS_BASE=~/popdata/progs/htslib/
SAMTOOLS_INC_DIR=${SAMTOOLS_BASE}
SAMTOOLS_LIB_DIR=${SAMTOOLS_BASE}
# actually -lhts -lz BUT -lz already specified separately
SAMTOOLS_LIB_SPEC=-lhts

COMPR_LIB_SPEC=-lz

SRC_DIR=src
INC_DIR=.
LIB_DIR=src
INSTALL_DIR=bin
BUILD_DIR=local-build

LIBS=${BUILD_DIR}/genomic_region.o ${BUILD_DIR}/repeat_spec.o ${BUILD_DIR}/allele.o \
     ${BUILD_DIR}/parameters.o ${BUILD_DIR}/ref_genome.o ${BUILD_DIR}/bam_file.o     \
     ${BUILD_DIR}/bam_index.o ${BUILD_DIR}/read_alignment.o   \
     ${BUILD_DIR}/irr_counting.o ${BUILD_DIR}/purity.o ${BUILD_DIR}/repeat_length.o \
     ${BUILD_DIR}/purity.o ${BUILD_DIR}/rep_align.o

all: ${INSTALL_DIR} ${INSTALL_DIR}/ExpansionHunter

${BUILD_DIR}:
	mkdir -p ${BUILD_DIR}

${BUILD_DIR}/%.o: purity/%.cc
	${CC} -c -o $@ $< -I ${INC_DIR}

${BUILD_DIR}/%.o: rep_align/%.cc
	${CC} -c -o $@ $< -I ${INC_DIR} -I ${BOOST_INC_DIR}

${BUILD_DIR}/%.o: ${LIB_DIR}/%.cc ${BUILD_DIR}
	${CC} -c -o $@ $< -I ${SAMTOOLS_INC_DIR} -I ${BOOST_INC_DIR} -I ${INC_DIR}

${INSTALL_DIR}:
	mkdir -p ${INSTALL_DIR}

${INSTALL_DIR}/ExpansionHunter: ${SRC_DIR}/expansion_hunter.cc ${LIBS}
	${CC} -o $@ $^ -I ${SAMTOOLS_INC_DIR} -I ${BOOST_INC_DIR} -I ${INC_DIR} -L ${SAMTOOLS_LIB_DIR} -L ${BOOST_LIB_DIR} ${SAMTOOLS_LIB_SPEC} ${COMPR_LIB_SPEC} ${BOOST_LIB_SPEC}
