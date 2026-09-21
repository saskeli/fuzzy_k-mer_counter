#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <seqan3/io/sequence_file/input.hpp>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "include/mer.hpp"

void help(char const* arg) {
  std::cout
      << "Counts number of occurrences of {ACGTN}^k -mers in fasta file(s).\n"
      << "Reads fasta file paths and mers from parameters, batch files\n"
      << "or standard input (if no paths are given) for piping.\n\n"
      << "Usage: " << arg << " [fasta path] [options] [mers ...]\n\n"
      << "fasta_file     Path to fasta file (possibly compressed).\n"
      << "               Needs to be placed before any mer parameters.\n"
      << "<mers ..  >    List of mers to search for.\n"
      << "-b <file_path> Reads file containing one fasta_file and\n"
      << "               at least one mer per line. Overrides <fasta_path>\n"
      << "               and <mers ...>.\n"
      << "-r             Do not include reverse complements for\n"
      << "               non-palindromic mers.\n"
      << "-s             Report at most one match per fasta record.\n"
      << "-t N           Set the maximum number of threads.\n"
      << "               Only runs multithreaded on batch jobs.\n"
      << "-p <path>      Prefix path to add to (each) fasta_file.\n"
      << "-h             Display this message and exit." << std::endl;
}

struct options {
  std::filesystem::path prefix;
  int threads;
  bool reverse_complement;
  bool single_match_per_read;
};

template <bool single_match, class seq_t, class mer_v_t, class count_t>
void count_matches(seq_t& fin, mer_v_t& mers, count_t& counts) {
  std::vector<bool> seq_count;
  if constexpr (single_match) {
    seq_count = std::vector<bool>(mers.size());
  }
  for (auto& record : fin) {
    auto seqan_seq = record.sequence();
    uint64_t len = 0;
    uint64_t sq = 0;
    for (auto base : seqan_seq) {
      sq = (sq << 2) | fkc::internal::nuc_to_v[seqan3::to_char(base)];
      ++len;
      for (size_t i = 0; i < mers.size(); ++i) {
        if (len >= mers[i].length()) {
          if constexpr (single_match) {
            seq_count[i] = seq_count[i] || mers[i].match(sq);
          } else {
            counts[i] += mers[i].match(sq);
          }
        }
      }
    }
    if constexpr (single_match) {
      for (size_t i = 0; i < mers.size(); ++i) {
        counts[i] += seq_count[i];
        seq_count[i] = false;
      }
    }
  }
}

template <class mer_v_t, class count_v_t>
void output_results(const std::filesystem::path& fasta_s, mer_v_t& mers, count_v_t& counts,
                    bool do_rc) {
  std::cout << fasta_s << " mer counts:\n";
  for (size_t i = 0; i < mers.size(); ++i) {
    fkc::mer mer = mers[i];
    uint64_t count = counts[i];
    if (do_rc && i + 1 < mers.size() && mers[i].rc() == mers[i + 1]) {
      count += counts[++i];
    }
    mer.print(std::cout) << "\t" << count << "\n";
  }
  std::cout << std::flush;
}

template <class mers_t>
void single_run(const std::filesystem::path& fasta_path, const mers_t& mers,
                const options& opt) {
  seqan3::sequence_file_input fin{fasta_path};
  std::vector<uint64_t> counts(mers.size());
  if (opt.single_match_per_read) {
    count_matches<true>(fin, mers, counts);
  } else {
    count_matches<false>(fin, mers, counts);
  }
#pragma omp critical
  {
    output_results(fasta_path, mers, counts, opt.reverse_complement);
  }
}

template <class in_t>
void run_batch(in_t& in, const options& opt) {
#ifdef _OPENMP
  omp_set_num_threads(opt.threads);
#endif

#pragma omp parallel
  {
    std::string line;
    while (true) {
      bool ok = true;
#pragma omp critical
      {
        if (in.bad() || in.eof()) {
          ok = false;
        } else {
          std::getline(in, line);
        }
      }
      if (not ok) {
        break;
      }
      std::istringstream iss(line);
      std::string fasta_path;
      std::string mer_s;
      std::vector<fkc::mer> mers;
      iss >> fasta_path;
      while (not iss.bad() && not iss.eof()) {
        iss >> mer_s;
        mers.push_back({mer_s.c_str()});
        if (opt.reverse_complement) {
          fkc::mer rc = mers.back().rc();
          if (rc != mers.back()) {
            mers.push_back(rc);
          }
        }
      }
      single_run(opt.prefix / fasta_path, mers, opt);
    }
  }
}

int main(int argc, char const* argv[]) {
  options opt{"", 1, true, false};
#ifdef _OPENMP
  opt.threads = omp_get_max_threads();
#endif
  int fasta_i = 0;
  int batch_i = 0;
  std::vector<fkc::mer> mers;

  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-h") == 0) {
      help(argv[0]);
      return 0;
    } else if (std::strcmp(argv[i], "-r") == 0) {
      opt.reverse_complement = false;
    } else if (std::strcmp(argv[i], "-s") == 0) {
      opt.single_match_per_read = true;
    } else if (std::strcmp(argv[i], "-t") == 0) {
      opt.threads = std::stol(argv[++i]);
    } else if (std::strcmp(argv[i], "-b") == 0) {
      batch_i = ++i;
    } else if (std::strcmp(argv[i], "-p") == 0) {
      opt.prefix = argv[++i];
    } else if (fasta_i == 0) {
      fasta_i = i;
    } else {
      mers.push_back({argv[i]});
      fkc::mer rc = mers.back().rc();
      if (opt.reverse_complement && mers.back() != rc) {
        mers.push_back(rc);
      }
    }
  }

  if (batch_i) {
    std::cerr << "Reading batch data from " << argv[batch_i] << std::endl;
    std::ifstream in_file(argv[batch_i]);
    run_batch(in_file, opt);
  } else if (fasta_i != 0) {
    single_run(opt.prefix / argv[fasta_i], mers, opt);
  } else {
    std::cerr << "Reading batch data from standard input..." << std::endl;
    run_batch(std::cin, opt);
  }

  return 0;
}
