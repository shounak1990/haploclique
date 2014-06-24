/* Copyright 2012,2014 Tobias Marschall
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

#include <sstream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "ShortDnaSequence.h"

#include "SortedBamReader.h"
#include "BamHelper.h"

using namespace std;

SortedBamReader::SortedBamReader(const std::string& filename, bool expand_xa) {
	this->counter = 0;
	this->expand_xa = expand_xa;
	this->progress_messages_os = 0;
	this->progress_message_frequency = -1;
	this->skipped_duplicates = 0;
	if (!bam_reader.Open(filename)) {
		ostringstream oss;
		oss << "Could not read BAM input from \"" << filename << "\".";
		throw std::runtime_error(oss.str());
	}
	readAllReads();
}

SortedBamReader::~SortedBamReader() {
	for (size_t i=0; i<read_list.size(); ++i) {
		free_record(i);
	}
}

void SortedBamReader::readAllReads() {
	while (true) {
		auto_ptr<BamTools::BamAlignment> read_aln(new BamTools::BamAlignment());
		if (!bam_reader.GetNextAlignment(*read_aln)){
			break;
		}
		// Check whether read name is already known
		read_name_map_t::const_iterator name_map_it = read_name_map.find(read_aln->Name);
		int read_idx = -1;
		// if not, create new record
		if (name_map_it == read_name_map.end()) {
			read_list.push_back(read_record_t());
			read_idx = read_list.size() - 1;
			read_name_map[read_aln->Name] = read_idx;
			read_list[read_idx].single_end = (!read_aln->IsFirstMate() && !read_aln->IsSecondMate());
		} else {
			read_idx = read_name_map[read_aln->Name];
		}
		read_record_t& record = read_list[read_idx];
		// Test whether flags are consistent with single/paired end status of record
		if (record.single_end) {
			if (read_aln->IsFirstMate() || read_aln->IsSecondMate()) {
				ostringstream oss;
				oss << "Invalid flags for read \"" << read_aln->Name << "\": expected single-end read.";
				throw std::runtime_error(oss.str());
			}
		} else {
			if (read_aln->IsFirstMate() == read_aln->IsSecondMate()) {
				ostringstream oss;
				oss << "Invalid flags for read \"" << read_aln->Name << "\": expected paired-end read.";
				throw std::runtime_error(oss.str());
			}
		}
		// check XA tag
		uint32_t x0 = 0;
		uint32_t x1 = 0;
		if (read_aln->GetTag("X0", x0) && read_aln->GetTag("X1", x1)) {
			string xa = "";
			if (!read_aln->GetTag("XA", xa) && (x0+x1>1)) {
				record.multiplyMappedButNoXA = true;
			}
		}
		counter += 1;
		// Print progress message if requested
		if ((progress_messages_os != 0) && (counter % progress_message_frequency == 0)) {
			*progress_messages_os << "Having processed " << counter << " read alignments" << endl;
		}
		if (read_aln->IsFailedQC() || read_aln->IsDuplicate()) {
			// cerr << "Warning: FailedQC or duplicated flag for \"" << read_aln->Name << "\" --> skipped." << endl;
			continue;
		}
		// Test whether we need to update position field
		if (read_aln->IsPrimaryAlignment() && ((record.position == -1) || (read_aln->Position < record.position))) {
			record.position = read_aln->Position;
		}
		// Store alignment in read record
		if (record.single_end || read_aln->IsFirstMate()) {
			BamTools::BamAlignment* a = read_aln.release();
			record.alignments1.push_back(a);
			record.coordinates1.insert(BamHelper::alignment_coordinate_t(*a));
			if (expand_xa) {
				BamHelper::expandXA(bam_reader, *a, &(record.alignments1), &(record.coordinates1), &skipped_duplicates);
			}
			if (a->IsMapped()) record.unmapped1 = false;
		} else {
			BamTools::BamAlignment* a = read_aln.release();
			record.alignments2.push_back(a);
			record.coordinates2.insert(BamHelper::alignment_coordinate_t(*a));
			if (expand_xa) {
				BamHelper::expandXA(bam_reader, *a, &(record.alignments2), &(record.coordinates2), &skipped_duplicates);
			}
			if (a->IsMapped()) record.unmapped2 = false;
		}
	}
	// We don't need the name to index map anymore.
	read_name_map.clear();
	// sort read records by position
	sort(read_list.begin(), read_list.end(), read_record_comparator_t());
	// set index to 0
	next_index = 0;
}

void SortedBamReader::free_record(size_t i) {
	assert(i < read_list.size());
	read_record_t& record = read_list[i];
	for (size_t j=0; j<record.alignments1.size(); ++j) {
		delete record.alignments1[j];
	}
	for (size_t j=0; j<record.alignments2.size(); ++j) {
		delete record.alignments2[j];
	}
	record.alignments1.clear();
	record.alignments2.clear();
}

void SortedBamReader::advance() {
	assert(next_index + 1 < (int)read_list.size());
	free_record(next_index);
	next_index += 1;
}

bool SortedBamReader::hasNext() const {
	return next_index + 1 < (int)read_list.size();
}

bool SortedBamReader::isSingleEnd() const {
	return read_list[next_index].single_end;
}

const std::string& SortedBamReader::getReadName() const {
	return read_list[next_index].alignments1[0]->Name;
}

const std::vector<BamTools::BamAlignment*>& SortedBamReader::getAlignments() const {
	assert(read_list[next_index].single_end);
	return read_list[next_index].alignments1;
}

const vector<BamTools::BamAlignment*>& SortedBamReader::getAlignmentsFirst() const {
	assert(!read_list[next_index].single_end);
	return read_list[next_index].alignments1;
}

const vector<BamTools::BamAlignment*>& SortedBamReader::getAlignmentsSecond() const {
	assert(!read_list[next_index].single_end);
	return read_list[next_index].alignments2;
}

auto_ptr<vector<BamTools::BamAlignment*> >  SortedBamReader::releaseAlignments() {
	assert(read_list[next_index].single_end);
	auto_ptr<vector<BamTools::BamAlignment*> > result(new vector<BamTools::BamAlignment*>(read_list[next_index].alignments1));
	read_list[next_index].alignments1.clear();
	return result;
}

auto_ptr<vector<BamTools::BamAlignment*> >  SortedBamReader::releaseAlignmentsFirst() {
	assert(!read_list[next_index].single_end);
	auto_ptr<vector<BamTools::BamAlignment*> > result(new vector<BamTools::BamAlignment*>(read_list[next_index].alignments1));
	read_list[next_index].alignments1.clear();
	return result;
}

auto_ptr<vector<BamTools::BamAlignment*> > SortedBamReader::releaseAlignmentsSecond() {
	assert(!read_list[next_index].single_end);
	auto_ptr<vector<BamTools::BamAlignment*> > result(new vector<BamTools::BamAlignment*>(read_list[next_index].alignments2));
	read_list[next_index].alignments2.clear();
	return result;
}

bool SortedBamReader::isUnmapped() const {
	assert(read_list[next_index].single_end);
	return read_list[next_index].unmapped1;
}

bool SortedBamReader::isFirstUnmapped() const {
	assert(!read_list[next_index].single_end);
	return read_list[next_index].unmapped1;
}

bool SortedBamReader::isSecondUnmapped() const {
	assert(!read_list[next_index].single_end);
	return read_list[next_index].unmapped2;
}

bool SortedBamReader::hasMultipleMappingsButNoXA() const {
	return read_list[next_index].multiplyMappedButNoXA;
}

void SortedBamReader::enableProgressMessages(std::ostream& os, int frequency) {
	assert(frequency > 0);
	this->progress_messages_os = &os;
	this->progress_message_frequency = frequency;
}

long long SortedBamReader::getSkippedDuplicates() const {
	return skipped_duplicates;
}

