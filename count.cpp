#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <seqan3/io/sequence_file/input.hpp>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "include/mer.hpp"
#include "include/seq.hpp"

void help(char const* arg, int t) {
  std::cout
      << "Counts number of occurrences of {ACGTN}^k -mers in fasta file(s).\n"
      << "Reads fasta file paths and mers from parameters, batch files\n"
      << "or standard input (if no paths are given) for piping.\n\n"
      << "For batch operation, maximum k is 32, unless run with the [-L]\n"
      << "option for max k = 64. For non-batch operation max k is 64.\n\n"
      << "Usage: " << arg << " -h\n"
      << "       " << arg << " [fasta path] [-b file_path] [-r] [-s] [-t N]\n"
      << "       [-L] [-p path] [mers ...]\n\n"
      << "Options:\n"
      << "fasta_file     Path to fasta file (possibly compressed).\n"
      << "               Needs to be placed before any mer parameters.\n"
      << "mers ...       List of mers to search for.\n"
      << "-b file_path   Reads file containing one fasta_file and\n"
      << "               at least one mer per line. Overrides [fasta_path]\n"
      << "               and [mers ...].\n"
      << "-r             Do not include reverse complements for\n"
      << "               non-palindromic mers.\n"
      << "-s             Report at most one match per fasta record.\n"
      << "-t N           Set the maximum number of threads. Default " << t 
      << ".\n"
      << "               Only runs multithreaded on batch jobs.\n"
      << "-p path        Prefix path to add to (each) fasta_file.\n"
      << "-L             Extend max k to 64, by sacrificing some performance.\n"
      << "-h             Display this message and exit." << std::endl;
}

struct options {
  std::filesystem::path prefix;
  int threads;
  bool reverse_complement;
  bool single_match_per_read;
  bool large_k;
};

template <bool single_match, class seq_t, class rec_t, class mer_v_t,
          class count_t>
void count_matches(rec_t& fin, mer_v_t& mers, count_t& counts) {
  std::vector<bool> seq_count;
  if constexpr (single_match) {
    seq_count = std::vector<bool>(mers.size());
  }
  for (auto& record : fin) {
    auto seqan_seq = record.sequence();
    seq_t seq;
    for (auto base : seqan_seq) {
      seq.append(seqan3::to_char(base));
      for (size_t i = 0; i < mers.size(); ++i) {
        if (seq.length() >= mers[i].length()) {
          if constexpr (single_match) {
            seq_count[i] = seq_count[i] || mers[i].match(seq);
          } else {
            counts[i] += mers[i].match(seq);
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

template <bool reverse_complement, class mer_v_t, class count_v_t>
void output_results(const std::filesystem::path& fasta_s, mer_v_t& mers,
                    count_v_t& counts) {
  std::cout << fasta_s << " mer counts:\n";
  for (size_t i = 0; i < mers.size(); ++i) {
    auto mer = mers[i];
    uint64_t count = counts[i];
    if constexpr (reverse_complement) {
      if (i + 1 < mers.size() && mers[i].rc() == mers[i + 1]) {
        count += counts[++i];
      }
    }
    mer.print(std::cout) << "\t" << count << "\n";
  }
  std::cout << std::flush;
}

template <bool single_match, bool reverse_complement, class mers_t>
void single_run(const std::filesystem::path& fasta_path, const mers_t& mers) {
  typedef typename std::conditional<
      std::is_same<typename mers_t::value_type, fkc::mer>::value, fkc::sequence<1>,
      fkc::sequence<2>>::type seq_t;
  seqan3::sequence_file_input fin{fasta_path};
  std::vector<uint64_t> counts(mers.size());
  count_matches<single_match, seq_t>(fin, mers, counts);
#pragma omp critical
  {
    output_results<reverse_complement>(fasta_path, mers, counts);
  }
}

template <bool long_mers, class in_t>
void run_batch(in_t& in, const options& opt) {
  typedef typename std::conditional<long_mers, fkc::long_mer, fkc::mer>::type mer_t;
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
      std::vector<mer_t> mers;
      iss >> fasta_path;
      while (not iss.bad() && not iss.eof()) {
        iss >> mer_s;
        mers.push_back({mer_s.c_str()});
        if (opt.reverse_complement) {
          mer_t rc = mers.back().rc();
          if (rc != mers.back()) {
            mers.push_back(rc);
          }
        }
      }
      if (opt.single_match_per_read) {
        if (opt.reverse_complement) {
          single_run<true, true>(opt.prefix / fasta_path, mers);
        } else {
          single_run<true, false>(opt.prefix / fasta_path, mers);
        }
      } else {
        if (opt.reverse_complement) {
          single_run<false, true>(opt.prefix / fasta_path, mers);
        } else {
          single_run<false, false>(opt.prefix / fasta_path, mers);
        }
      }
    }
  }
}

int main(int argc, char const* argv[]) {
  options opt{"", 1, true, false, false};
#ifdef _OPENMP
  opt.threads = omp_get_max_threads();
#endif
  int fasta_i = 0;
  int batch_i = 0;
  std::vector<fkc::long_mer> mers;

  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-h") == 0) {
      help(argv[0], opt.threads);
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
    } else if (std::strcmp(argv[i], "-L") == 0) {
      opt.large_k = true;
    } else if (fasta_i == 0) {
      fasta_i = i;
    } else {
      mers.push_back({argv[i]});
      fkc::long_mer rc = mers.back().rc();
      if (opt.reverse_complement && mers.back() != rc) {
        mers.push_back(rc);
      }
    }
  }

  if (batch_i) {
    std::cerr << "Reading batch data from " << argv[batch_i] << std::endl;
    std::ifstream in_file(argv[batch_i]);
    if (opt.large_k) {
      run_batch<true>(in_file, opt);
    } else {
      run_batch<false>(in_file, opt);
    }
    
  } else if (fasta_i != 0) {
    if (opt.single_match_per_read) {
      if (opt.reverse_complement) {
        single_run<true, true>(opt.prefix / argv[fasta_i], mers);
      } else {
        single_run<true, true>(opt.prefix / argv[fasta_i], mers);
      }
    } else {
      if (opt.reverse_complement) {
        single_run<false, true>(opt.prefix / argv[fasta_i], mers);
      } else {
        single_run<false, true>(opt.prefix / argv[fasta_i], mers);
      }
    }
    
  } else {
    std::cerr << "Reading batch data from standard input..." << std::endl;
    if (opt.large_k) {
      run_batch<true>(std::cin, opt);
    } else {
      run_batch<false>(std::cin, opt);
    }
    
  }

  return 0;
}
