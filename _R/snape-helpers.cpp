#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>



//[[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

template<typename T> using vector2d = std::vector<std::vector<T>>;
template<typename T> using vector3d = std::vector<std::vector<std::vector<T>>>;




// Get allele counts from gSYNC count outputs
//[[Rcpp::export]]
List split_count_strings(const std::vector<std::string>& strings) {

    List out(strings.size());

    // Calculate number of sub-strings separated by ':':
    size_t nc = 0;
    for (const char& c : strings[0]) if (c == ':') nc++;
    nc++;

    std::string tmp0;
    std::vector<int> tmp1(nc);

    for (size_t i = 0; i < strings.size(); i++) {
        const std::string& s(strings[i]);
        std::istringstream iss(s);
        size_t j = 0;
        while (std::getline(iss, tmp0, ':')) {
            tmp1[j] = std::stoi(tmp0);
            j++;
        }
        out[i] = tmp1;
    }

    return out;
}


//' Get allele frequencies from string outputs from SNAPE
//[[Rcpp::export]]
std::vector<double> get_snape_af(const std::vector<std::string>& strings) {

    std::vector<double> out;
    out.reserve(strings.size());

    // Calculate number of sub-strings separated by ':':
    size_t nc = 0;
    for (const char& c : strings[0]) if (c == ':') nc++;
    nc++;

    std::string tmp0;
    std::vector<std::string> tmp1(nc);

    for (size_t i = 0; i < strings.size(); i++) {
        const std::string& s(strings[i]);
        std::istringstream iss(s);
        size_t j = 0;
        while (std::getline(iss, tmp0, ':')) {
            tmp1[j] = tmp0;
            j++;
        }
        out.push_back(std::stod(tmp1.back()));
    }

    return out;
}




//' This function calculates the number of alleles at each locus
//[[Rcpp::export]]
std::vector<int> count_alleles(const vector3d<unsigned>& allele_counts) {

 size_t n_pools = allele_counts.size();
 size_t n_alleles = allele_counts[0].size();
 if (n_pools == 0) return std::vector<int>();
 if (n_alleles == 0) return std::vector<int>();

 std::vector<int> out;
 out.reserve(n_alleles);

 for (size_t j = 0; j < n_alleles; j++) {
     out.push_back(0);
     int& n_nonzeros(out.back());
     for (size_t k = 0; k < 6; k++) {
         unsigned sum_k = 0;
         for (size_t i = 0; i < n_pools; i++) {
             sum_k += allele_counts[i][j][k];
         }
         if (sum_k > 0) n_nonzeros++;
     }
 }

 return out;
}









inline int one_which_allele(const std::vector<unsigned>& cnt,
                            const double& frq) {

    double cnt_sum = static_cast<double>(std::accumulate(cnt.begin(),
                                                         cnt.end(), 0));
    double cnt_frq;
    std::vector<double> cnt_frq_diffs;
    cnt_frq_diffs.reserve(cnt.size());
    for (const unsigned& x : cnt) {
        cnt_frq = static_cast<double>(x) / cnt_sum;
        cnt_frq_diffs.push_back(std::abs(cnt_frq - frq));
    }

    int mdi = std::min_element(cnt_frq_diffs.begin(), cnt_frq_diffs.end()) -
        cnt_frq_diffs.begin();
    // Because output should be 1-based indices.
    mdi++;
    return mdi;

}




//' This function determines which allele the frequencies must be referring to.
//' This is a sanity check to make sure they don't vary between pools.
//'
//' Only recommended for biallelic loci!
//'
//[[Rcpp::export]]
IntegerMatrix which_allele(const vector3d<unsigned>& allele_counts,
                           NumericMatrix allele_freqs) {

 size_t n_pools = allele_counts.size();
 size_t n_alleles = allele_counts[0].size();

 IntegerMatrix out(n_alleles, n_pools);

 if (n_pools == 0) return out;
 if (n_alleles == 0) return out;
 if (allele_freqs.ncol() != n_pools) {
     stop("allele_freqs.ncol() != n_pools");
 }
 if (allele_freqs.nrow() != n_alleles) {
     stop("allele_freqs.nrow() != n_alleles");
 }

 int tmp;

 for (size_t i = 0; i < n_pools; i++) {
     for (size_t j = 0; j < n_alleles; j++) {
         tmp = one_which_allele(allele_counts[i][j],
                                allele_freqs(j, i));
         out(j, i) = tmp;
     }
 }

 return out;
}



//' Below function filters positions so that they're all at least `min_dist`
//' from each other, and it does this within sequences.
//[[Rcpp::export]]
LogicalVector min_dist_filter(CharacterVector seq,
                              IntegerVector pos,
                              const int& min_dist) {
    size_t n_pos = pos.size();
    if (seq.size() != n_pos) stop("seq.size() != pos.size()");
    LogicalVector out(n_pos, 0);
    if (n_pos == 0) return out;
    out[0] = true;
    int lpos = pos[0];
    String lseq = seq[0];
    size_t i = 0;
    while (true) {
        while (i < n_pos && lseq == seq[i] && (pos[i] - lpos) < min_dist) {
            i++;
        }
        if (i >= n_pos) break;
        out[i] = true;
        lpos = pos[i];
        if (lseq != seq[i]) lseq = seq[i];
        i++;
    }
    return out;
}
