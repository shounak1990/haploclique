/* Copyright 2012-2014 Tobias Marschall and Armin Töpfer
 *
 * This file is part of HaploClique.
 *
 * HaploClique is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HaploClique is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HaploClique.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cassert>
#include <iomanip>
#include <ctime>

#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <boost/tokenizer.hpp>
#include <bamtools/api/BamReader.h>
#include <bamtools/api/BamWriter.h>

#include "HistogramBasedDistribution.h"
#include "SortedBamReader.h"
#include "Histogram.h"
#include "ThreadPool.h"

using namespace std;
namespace po = boost::program_options;

int phred_base = 33;
long long skipped_inside_clipping = 0;

typedef boost::unordered_map<std::string, HistogramBasedDistribution*> readgroupwise_distributions_t;

void usage(const char* name, const po::options_description& options_desc) {
	cerr << "Usage: " << name << " [options] <reference.fasta(.gz)> <input.bam>" << endl;
	cerr << endl;
	cerr << "Reads a sorted BAM file (also parsing the BWA-specific XA tags)" << endl;
	cerr << "and outputs a list of plausible alignments (one per line) of read pairs" << endl;
	cerr << "to stdout in the following format:" << endl;
	cerr << endl;
	cerr << "<read-name> <pair-nr> <read-group> <phred-sum1> <chromosome1> <start1> <end1> <strand1> <cigar1> <seq1> <qualities1> <phred-sum2> <chromosome2> <start2> <end2> <strand2> <cigar2> <seq2> <qualities2> <aln-pair-prob> <aln-pair-prob-inslength>" << endl;
	cerr << endl;
	cerr << "Where <aln-pair-prob> is the probability that this alignment pair is " << endl;
	cerr << "correct estimated based on the alignment quality while the probabilities " << endl;
	cerr << "given in <aln-pair-prob-inslength> are estimates based on alignment quality AND" << endl;
	cerr << "insert length." << endl;
	cerr << endl;
	cerr << "<startX> and <endX> coordinates are 1-based and inclusive (i.e. closed intervals)" << endl;
	cerr << "and give the region on the reference the respective read was aligned to." << endl;
	cerr << endl;
	cerr << "In case of single end alignments, the format is changed to:" << endl;
	cerr << "<read-name> <read-nr> <read-group> <phred-sum> <chromosome> <start> <end> <strand> <cigar> <seq> <qualities> <aln-prob>" << endl;
	cerr << endl;
	cerr << options_desc << endl;
	exit(1);
}

// TODO: some duplicate code shared with bam-to-alignment-priors and add-score-tags-to-bam

bool is_clipped(const BamTools::BamAlignment& alignment) {
	vector<BamTools::CigarOp>::const_iterator it = alignment.CigarData.begin();
	for (;it!=alignment.CigarData.end(); ++it) {
		switch (it->Type) {
		case 'S':
		case 'H':
		case 'P':
			return true;
		default:
			continue;
		}
	}
	return false;
}

bool is_unique(const BamTools::BamAlignment& alignment) {
	if (alignment.MapQuality == 0) return false;
	if (!alignment.IsPrimaryAlignment()) return false;
	uint32_t x0_tag = -1;
	uint32_t x1_tag = -1;
	if (alignment.GetTag("X0",x0_tag)) {
		if (x0_tag>1) return false;
	}
	if (alignment.GetTag("X1",x1_tag)) {
		if (x1_tag>0) return false;
	}
	return true;
}

/** Returns probability "distribution" according to MAPQ fields. Returns 0 if MAPQ is unavailible. */
vector<double>* mapq_probabilities(const vector<BamTools::BamAlignment*>& alignments) {
	typedef vector<BamTools::BamAlignment*>::const_iterator aln_iter_t;
	for (aln_iter_t it=alignments.begin(); it!=alignments.end(); ++it) {
		if ((*it)->MapQuality == 255) return 0;
	}
	double probsum = 0.0;
	vector<double>* probs = new vector<double>();
	for (aln_iter_t it=alignments.begin(); it!=alignments.end(); ++it) {
		if ((*it)->MapQuality == 0) {
			probs->push_back(0.0);
		} else {
			double p = 1.0 - pow(10.0,(*it)->MapQuality/-10.0);
			probs->push_back(p);
			probsum += p;
		}
	}
	if (probsum > 1.0) {
		for (size_t i=0; i<probs->size(); ++i) {
			probs->at(i) /= probsum;
		}
	}
	return probs;
}

typedef struct aln_pair_t {
	const BamTools::BamAlignment* aln1;
	const BamTools::BamAlignment* aln2;
	double pair_prob;
	aln_pair_t(const BamTools::BamAlignment* aln1, const BamTools::BamAlignment* aln2, double pair_prob) : aln1(aln1), aln2(aln2), pair_prob(pair_prob) {}
} aln_pair_t;

/** Given an output file (a plain file object) and a list of alignments belonging
to the same read pair, outputs one line per alignment pair. */
void process_read_allpairs(const BamTools::RefVector& ref_vector, const vector<BamTools::BamAlignment*>& alignments1, const vector<BamTools::BamAlignment*>& alignments2, const HistogramBasedDistribution& ild) {
	vector<double>* posteriors1;
	vector<double>* posteriors2;
	posteriors1 = mapq_probabilities(alignments1);
	if (posteriors1 == 0) return; //Fehler ausgeben, wenn keine MAPQ vorhanden ist
	posteriors2 = mapq_probabilities(alignments2);
	if (posteriors2 == 0) {
		delete posteriors1;
		return; //Fehler ausgeben, wenn keine MAPQ vorhanden ist
	}
	assert(posteriors1->size() == alignments1.size());
	assert(posteriors2->size() == alignments2.size());
	// contains an entry for every pair that should be printed
	vector<aln_pair_t> all_pairs;
	double p_sum = 0.0;
	for (size_t i=0; i<alignments1.size(); ++i) {
		const BamTools::BamAlignment& read_aln1 = *alignments1[i];
		assert(read_aln1.CigarData.size() > 0);
		double posterior1 = posteriors1->at(i);
		assert(read_aln1.IsMapped());
		for (size_t j=0; j<alignments2.size(); ++j) {
			const BamTools::BamAlignment& read_aln2 = *alignments2[j];
			assert(read_aln2.CigarData.size() > 0);
			double posterior2 = posteriors2->at(j);
			assert(read_aln2.IsMapped());
			// skip if strandedness is the same
			if (read_aln1.IsReverseStrand() == read_aln2.IsReverseStrand()) continue;
			// skip if both reads map to different chromosomes
			if (read_aln1.RefID != read_aln2.RefID) continue;
			double p = posterior1*posterior2;
			
			const BamTools::BamAlignment& left = (read_aln1.Position<=read_aln2.Position)?read_aln1:read_aln2;
			const BamTools::BamAlignment& right = (read_aln1.Position<=read_aln2.Position)?read_aln2:read_aln1;
			int insert_length = right.Position - left.GetEndPosition();
			p_sum += p;
			all_pairs.push_back(aln_pair_t(&left, &right, p));
		}
	}
	string rg = "";
	if (!alignments1[0]->GetTag("RG", rg)) {
		rg = "-1";
	}
	if (p_sum > 0.0) {
		int n = 0;
		for (size_t i=0; i<all_pairs.size(); ++i) {
			aln_pair_t& p = all_pairs[i];
			// cout << p.aln1->Name << " " << (n++) << " " << rg << " "
			//      << phred_sum(*(p.aln1)) << " " << ref_vector[p.aln1->RefID].RefName << " " << p.aln1->Position+1 << " " << p.aln1->GetEndPosition() << " " << (p.aln1->IsReverseStrand()?'-':'+') << " " << p.aln1->CigarData << " " << p.aln1->QueryBases << " " << p.aln1->Qualities << " "
			//      << phred_sum(*(p.aln2)) << " " << ref_vector[p.aln2->RefID].RefName << " " << p.aln2->Position+1 << " " << p.aln2->GetEndPosition() << " " << (p.aln2->IsReverseStrand()?'-':'+') << " " << p.aln2->CigarData << " " << p.aln2->QueryBases << " " << p.aln2->Qualities << " "
		 //         << setprecision(16) << p.pair_prob/p_sum << endl;
		}
	}
	delete posteriors1;
	delete posteriors2;
}

void process_read_single_end(const BamTools::RefVector& ref_vector, const vector<BamTools::BamAlignment*>& alignments1) {
	vector<double>* posteriors;

	posteriors = mapq_probabilities(alignments1);
	if (posteriors == 0) return;

	assert(posteriors->size() == alignments1.size());
	// contains an entry for every pair that should be printed
	vector<aln_pair_t> all_pairs;
	double p_sum = 0.0;
	for (size_t i=0; i<alignments1.size(); ++i) {
		const BamTools::BamAlignment& aln = *alignments1[i];
		assert(aln.CigarData.size() > 0);
		assert(aln.IsMapped());
		p_sum += posteriors->at(i);
	}
	string rg = "";
	if (!alignments1[0]->GetTag("RG", rg)) {
		rg = "-1";
	}
	if (p_sum > 0.0) {
		for (size_t i=0; i<alignments1.size(); ++i) {
			const BamTools::BamAlignment& aln = *alignments1[i];
			// out_stream << aln.Name << " " << i << " " << rg << " "
			// 	<< phred_sum(aln) << " " << ref_vector[aln.RefID].RefName << " " << aln.Position+1 << " " << aln.GetEndPosition() << " " << (aln.IsReverseStrand()?'-':'+') << " " << aln.CigarData << " " << aln.QueryBases << " " << aln.Qualities << " "
			// 	<< setprecision(16) << posteriors->at(i)/p_sum << endl;
		}
	}
	delete posteriors;
}

void read_readgroup_list(const string& filename, readgroupwise_distributions_t* result) {
	ifstream is(filename.c_str());
	if (is.fail()) {
		cerr << "Error: could not open \"" << filename << "\"." << endl;
		exit(1);
	}
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
	boost::char_separator<char> whitespace_separator(" \t");
	string line;
	int linenr = 0;
	while (getline(is,line)) {
		linenr += 1;
		tokenizer_t tokenizer(line, whitespace_separator);
		vector<string> tokens(tokenizer.begin(), tokenizer.end());
		if (tokens.size() != 2) {
			ostringstream oss;
			oss << "Error parsing read group list. Offending line: " << linenr;
			throw std::runtime_error(oss.str());
		}
		const string& readgroup = tokens[0];
		readgroupwise_distributions_t::const_iterator it = result->find(readgroup);
		if (it != result->end()) {
			ostringstream oss;
			oss << "Duplicate readgroup \"" << readgroup << "\" in file \"" << filename << "\"." << endl;
			throw std::runtime_error(oss.str());
		}
		HistogramBasedDistribution* distribution = new HistogramBasedDistribution(tokens[1]);
		(*result)[readgroup] = distribution;
	}
}

int main(int argc, char* argv[]) {
	// PARAMETERS
	bool dont_skip_non_xa = false;
	int max_span;
	int distribution_estimation_count;
	string insert_length_dist_filename = "";
	string rgwise_insert_length_dist_filename = "";
	string mean_sd_filename = "";
	int bad_score_threshold;
	double discard_threshold;


	po::options_description options_desc("Allowed options");
	options_desc.add_options()
		("dont_skip_non_xa,x", po::value<bool>(&dont_skip_non_xa)->zero_tokens(), "Do not skip reads for which other alignments exist (i.e. X0+X1>1, but no XA tag is present).")
		("phred_base,p", po::value<int>(&phred_base)->default_value(33), "Value to substract from ASCII code to get the PHRED quality.")
		("bad_alignment_threshold,b", po::value<int>(&bad_score_threshold)->default_value(1000), "Issue a warning when AS tag is above this value.")
		("max_span,s", po::value<int>(&max_span)->default_value(50000), "Maximal internal segment. Read pairs with larger internal segment will be ignored.")
		("discard_reads,d", po::value<double>(&discard_threshold)->default_value(0.0), "Discard \"concordant\" alignments within the given number of standard deviations (default: disabled).")
		("insert_size_dist,i", po::value<string>(&insert_length_dist_filename)->default_value(""), "Filename of known internal segment size distribution. If not given, this distribution is estimated.")
		("rg_insert_size_dist,r", po::value<string>(&rgwise_insert_length_dist_filename)->default_value(""), "Filename of read-group-wise known internal segment size distributions. Expects two-column text file: <readgroup> <distribution-filename>.")
		("dist_est_count,n", po::value<int>(&distribution_estimation_count)->default_value(5000000), "Number of uniquely mapping reads that are to be used to estimate internal segment size distribution.")
		("mean_and_sd,m", po::value<string>(&mean_sd_filename), "Write (robustly estimated) mean and standard deviation of main peak if internal segment size distribution to given filename.")
	;

	if (argc<3) {
		usage(argv[0], options_desc);
	}
	string reference_filename(argv[argc-2]);
	string bam_input_filename(argv[argc-1]);
	argc -= 2;

	po::variables_map options;
	try {
		po::store(po::parse_command_line(argc, argv, options_desc), options);
		po::notify(options);
	} catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	}

	if ((discard_threshold != 0.0) && (insert_length_dist_filename.size()>0)) {
		cerr << "Options -i an -d cannot be used together." << endl;
		return 1;
	}

	if ((insert_length_dist_filename.size() > 0) && (rgwise_insert_length_dist_filename.size() > 0)) {
		cerr << "Options -i and -r cannot be used at the same time." << endl;
		return 1;
	}

	clock_t clock_start = clock();

	// Create insert size distribution. Either read it from file (if given) or
	// estimate based on input
	HistogramBasedDistribution* insert_length_distribution = 0;
	readgroupwise_distributions_t* readgroupwise_distributions = 0;
	long long skipped_by_xa = 0;
	auto_ptr<BamReader> bam_reader(0);

	try {
		cerr << "Estimating internal segment size distribution" << endl;
		bam_reader = auto_ptr<BamReader>(new SortedBamReader(bam_input_filename));
		const BamTools::RefVector& bam_ref_data = bam_reader->getReferenceData();
		Histogram histogram;
		int i = 0;
		while (bam_reader->hasNext()) {
			bam_reader->advance();
			bool single_end = 0;
			//TODO: change in SortedBamReader
			// single_end = bam_reader->isSingleEnd();
			if (single_end) {
				if (bam_reader->isUnmapped()) continue;
				if (!dont_skip_non_xa && bam_reader->hasMultipleMappingsButNoXA()) {
					skipped_by_xa += 1;
					continue;
				}
				vector<BamTools::BamAlignment*>* alignments = bam_reader->releaseAlignments().release();
				process_read_single_end(bam_ref_data, *alignments);
			} else {
				if (bam_reader->isFirstUnmapped() || bam_reader->isSecondUnmapped()) continue;
				if (!dont_skip_non_xa && bam_reader->hasMultipleMappingsButNoXA()) {
					skipped_by_xa += 1;
					continue;
				}
				vector<BamTools::BamAlignment*>* alignments1 = bam_reader->releaseAlignmentsFirst().release();
				vector<BamTools::BamAlignment*>* alignments2 = bam_reader->releaseAlignmentsSecond().release();
				
				assert(alignments1 != 0);
				// compute prior probabilities and prepare output
				process_read_allpairs(bam_ref_data, *alignments1, *alignments2, *insert_length_distribution);
				
				if (bam_reader->isFirstUnmapped() || bam_reader->isSecondUnmapped()) continue;
				if (bam_reader->hasMultipleMappingsButNoXA()) continue;
				assert(alignments1->size() > 0);
				assert(alignments2->size() > 0);
				assert(alignments1->at(0) != 0);
				assert(alignments2->at(0) != 0);
				if ((alignments1->size() == 1) && (alignments2->size() == 1))
				if (alignments1->at(0)->IsReverseStrand() == alignments2->at(0)->IsReverseStrand()) continue;
				if (alignments1->at(0)->RefID != alignments2->at(0)->RefID) continue;
				if (is_clipped(*alignments1->at(0)) || is_clipped(*alignments2->at(0))) continue;
				if (!is_unique(*alignments1->at(0)) || !is_unique(*alignments2->at(0))) continue;
				int insert_size = 0;
				if (alignments1->at(0)->Position <= alignments2->at(0)->Position) {
					insert_size = alignments2->at(0)->Position - alignments1->at(0)->GetEndPosition();
				} else {
					insert_size = alignments1->at(0)->Position - alignments2->at(0)->GetEndPosition();
				}
				histogram.add(insert_size);
			}
			i += 1;
		}
		if (i < 1) {
			cerr << "Error: Too few reads estimate internal segment size distribution (" << i << ")." << endl;
			cerr << "Wont proceed with less than 1 uniquely mappable read pairs." << endl;
			return 1;
		}
		if (i < distribution_estimation_count) {
			cerr << "Warning: fewer uniquely mappable read (" << i << ") than asked for (" << distribution_estimation_count << ", see option -n)." << endl;
		}
		double mean;
		double sd;
		histogram.computeMeanAndStddev(&mean, &sd);
		cerr << "Main peak of internal segment length distribution: mean " << mean << ", sd " << sd << endl;
		if (mean_sd_filename.compare("") != 0) {
			ofstream ofs(mean_sd_filename.c_str());
			ofs << mean << " " << sd << endl;
			ofs.close();
		}
		std::auto_ptr<HistogramBasedDistribution> d = histogram.toDistribution(20);
		assert(d.get() != 0);
		insert_length_distribution = d.release();
	} catch(exception& e) {
		cerr << "Error: " << e.what() << "\n";
		return 1;
	}


	if (skipped_by_xa > 0) {
		cerr << "Skipped " << skipped_by_xa << " ambiguously mapped reads for which no XA tag was present (to prevent this, use option -x)." << endl;
	}
	if (bam_reader->getSkippedDuplicates() > 0) {
		cerr << "Skipped " << bam_reader->getSkippedDuplicates() << " duplicate alignments." << endl;
	}
	if (skipped_inside_clipping > 0) {
		cerr << "Skipped " << skipped_inside_clipping << " reads with alignments soft-clipped on the inside of an alignment pair!" << endl;
	}

	if (insert_length_distribution!=0) delete insert_length_distribution;
	if (readgroupwise_distributions != 0) {
		readgroupwise_distributions_t::const_iterator it = readgroupwise_distributions->begin();
		for (; it != readgroupwise_distributions->end(); ++it) delete it->second;
		delete readgroupwise_distributions;
	}
	double cpu_time = (double)(clock() - clock_start) / CLOCKS_PER_SEC;
	cerr << "Total CPU time: " << cpu_time << endl;
	return 0;
}
