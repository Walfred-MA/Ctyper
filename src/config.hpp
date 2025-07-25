#ifndef config_hpp
#define config_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include <atomic>
typedef unsigned int uint;
typedef unsigned __int128 u128;
typedef unsigned long long ull;
typedef unsigned short  uint16;
typedef unsigned char  uint8;
typedef std::atomic<unsigned long long> ull_atom;

constexpr int background_prime = 1000003;
constexpr const char* vctyper = "v1.2.0";
#define klen 31

#define large_prime 2147483647
#define prime10M 9999991
#define MAX_UINT16 65535
#define MAX_COUNT 16383 //2^14 -1
#define MAX_LINE 2000000
#define MAX_UINT32 4294967295
#define MAX_UINT31 2147483647
#define Comb2( size ) (size + 1) * size / 2
#define MIN( A , B ) ( A <= B ) ? A : B
#define MAX( A , B ) ( A >= B ) ? A : B
#define MAXABS( A , B ) ( abs(A) >= abs(B) ) ? A : B
#define spair std::pair<std::string,std::string>
typedef float FLOAT_T;
typedef Eigen::MatrixXf Matrix_T ;
typedef Eigen::VectorXf Vector_T ;
constexpr ull REVNUM = ~0ULL ^ (3ULL << 62); //2^62-1


#define FLOAT_T float
#define FIXCOL 6
#define THREECOLLENGTH 17

#define genomesize_mean 6320012150.0
#define genomesize_male 6270605410.0
#define genomesize_female 6320012150.0

extern bool optioncorr;

#define errorcutoff1 15
#define errorcutoff2 30
#define varcutoff 30

#define sufficient 1000
#define corrstartpoint 0.1

#define windowmerge 15
#define minwindowcutoff 3
#define windowunit 30

#define maxiteration 2000

#define ABS(x) ((unsigned long long)(((x) >= 0) ? (x) : -(x)))

static constexpr std::array<uint16_t, 128> reciprocals = [] {
    std::array<uint16_t, 128> r{};
    for (int d = 1; d < 128; ++d) {
        r[d] = static_cast<uint16_t>((1ULL << 16) / d);
    }
    return r;
}();



static const uint8_t RECIPROCALS[84] = {200, 196, 192, 189, 185, 182, 179, 175, 172, 169, 167, 164, 161, 159, 156, 154, 152, 149, 147, 145, 143, 141, 139, 137, 135, 133, 132, 130, 128, 127, 125, 123, 122, 120, 119, 118, 116, 115, 114, 112, 111, 110, 109, 108, 106, 105, 104, 103, 102, 101, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67};
#define RECIPROCAL100 33
//
#define APPLY_OPTION_CORR(optioncorr, flag, count_f)                     \
    do {                                                                 \
        if ((optioncorr)  && (flag & 1)   ) {                             \
            (count_f) *= 0.01 * ((flag & 0xFF00) >> 8);      \
        }                                                                \
    } while (0)

#define APPLY_OPTION_CORR2(optioncorr, flag, count_f, new_weight)                         \
    do {                                                                                 \
        if ((optioncorr)  && (flag & 1)  ) {                                              \
            (count_f) *= 0.01 * ((flag & 0xFF00) >> 8);                                \
            (new_weight) = 0.01 * ((flag & 0x000000FE)) * 0.5;                       \
        }                                                                                \
    } while (0)

#define GET_CRAM_REF_PATH() \
    ((!refpath.empty()) ? refpath : \
    (std::getenv("SAMTOOLS_REF_PATH") ? std::getenv("SAMTOOLS_REF_PATH") : \
    (std::getenv("REF_PATH") ? std::getenv("REF_PATH") : "NOT SET (HTS_DETAULT)")))


using hotspot = std::pair<int,int>;

#define DEPTH_REGION {"chr1:9995000-10005000", "chr1:37995000-38005000", "chr1:65995000-66005000", "chr1:93995000-94005000", "chr1:121995000-122005000", "chr1:149995000-150005000", "chr1:177995000-178005000", "chr1:205995000-206005000", "chr1:233995000-234005000", "chr2:9995000-10005000", "chr2:37995000-38005000", "chr2:65995000-66005000", "chr2:93995000-94005000", "chr2:121995000-122005000", "chr2:149995000-150005000", "chr2:177995000-178005000", "chr2:205995000-206005000", "chr3:9995000-10005000", "chr3:37995000-38005000", "chr3:65995000-66005000", "chr3:93995000-94005000", "chr3:121995000-122005000", "chr3:149995000-150005000", "chr3:177995000-178005000", "chr4:9995000-10005000", "chr4:37995000-38005000", "chr4:65995000-66005000", "chr4:93995000-94005000", "chr4:121995000-122005000", "chr4:149995000-150005000", "chr4:177995000-178005000", "chr5:9995000-10005000", "chr5:37995000-38005000", "chr5:65995000-66005000", "chr5:93995000-94005000", "chr5:121995000-122005000", "chr5:149995000-150005000", "chr6:9995000-10005000", "chr6:37995000-38005000", "chr6:65995000-66005000", "chr6:93995000-94005000", "chr6:121995000-122005000", "chr6:149995000-150005000", "chr7:9995000-10005000", "chr7:37995000-38005000", "chr7:65995000-66005000", "chr7:93995000-94005000", "chr7:121995000-122005000", "chr8:9995000-10005000", "chr8:37995000-38005000", "chr8:65995000-66005000", "chr8:93995000-94005000", "chr8:121995000-122005000", "chr9:9995000-10005000", "chr9:37995000-38005000", "chr9:65995000-66005000", "chr9:93995000-94005000", "chr9:121995000-122005000", "chr10:9995000-10005000", "chr10:37995000-38005000", "chr10:65995000-66005000", "chr10:93995000-94005000", "chr10:121995000-122005000", "chr11:9995000-10005000", "chr11:37995000-38005000", "chr11:65995000-66005000", "chr11:93995000-94005000", "chr11:121995000-122005000", "chr12:9995000-10005000", "chr12:37995000-38005000", "chr12:65995000-66005000", "chr12:93995000-94005000", "chr12:121995000-122005000", "chr13:9995000-10005000", "chr13:37995000-38005000", "chr13:65995000-66005000", "chr13:93995000-94005000", "chr14:9995000-10005000", "chr14:37995000-38005000", "chr14:65995000-66005000", "chr14:93995000-94005000", "chr15:9995000-10005000", "chr15:37995000-38005000", "chr15:65995000-66005000", "chr16:9995000-10005000", "chr16:37995000-38005000", "chr16:65995000-66005000", "chr17:9995000-10005000", "chr17:37995000-38005000", "chr17:65995000-66005000", "chr18:9995000-10005000", "chr18:37995000-38005000", "chr18:65995000-66005000", "chr19:9995000-10005000", "chr19:37995000-38005000", "chr20:9995000-10005000", "chr20:37995000-38005000", "chr21:9995000-10005000", "chr22:9995000-10005000", "chr22:37995000-38005000"};



#endif /* config_hpp */
