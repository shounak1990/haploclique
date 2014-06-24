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

#ifndef SORTEDBAMREADER_H_
#define SORTEDBAMREADER_H_

#include <iostream>
#include <vector>
#include <map>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <bamtools/api/BamReader.h>

#include "BamReader.h"
#include "BamHelper.h"

/** Reads a sorted BAM file and returns groups of alignments belonging to
 *  the same read. If present, interprets XA tags as written by BWA.
 */
class SortedBamReader : public BamReader {
private:
	BamTools::BamReader bam_reader;

	typedef struct read_record_t {
		bool single_end;
		bool unmapped1;
		bool unmapped2;
		bool multiplyMappedButNoXA;
		int position; // position of the leftmost primary alignment
		std::vector<BamTools::BamAlignment*> alignments1;
		std::vector<BamTools::BamAlignment*> alignments2;
		boost::unordered_set<BamHelper::alignment_coordinate_t> coordinates1;
		boost::unordered_set<BamHelper::alignment_coordinate_t> coordinates2;
		read_record_t() : single_end(true), unmapped1(true), unmapped2(true), multiplyMappedButNoXA(false), position(-1) {}
		~read_record_t() {}
	} read_record_t;
	typedef struct read_record_comparator_t {
		bool operator()(const read_record_t& r1, const read_record_t& r2) { return r1.position > r2.position; }
	} read_record_comparator_t;

	// list of all reads
	std::vector<read_record_t> read_list;
	// map: read name --> list index
	typedef boost::unordered_map<std::string,int> read_name_map_t;
	read_name_map_t read_name_map;
	// count all processed alignments
	long long counter;
	std::ostream* progress_messages_os;
	int progress_message_frequency;
	
	int next_index;
	bool expand_xa;
	
	long long skipped_duplicates;

	void free_record(size_t i);
	void readAllReads();

public:
	/** If paired is true, then the input BAM is expected to be paired end data. In this case, the methods
	 *  getAlignmentsFirst and getAlignmentsSecond must be used (instead of getAlignments) to get the
	 *  alignments of first/second read end. */
	SortedBamReader(const std::string& filename, bool expand_xa = true);
	virtual ~SortedBamReader();

	/** Returns true if there is another read pair left in the input file. */
	virtual bool hasNext() const;

	/** Reads next read pair. Before calling this method, hasNext() should
	 *  be called. */
	virtual void advance();

	/** Returns the name of the currently processed read. */
	virtual const std::string& getReadName() const;

	/** Returns a reference to alignments current read. May only be called when
	 *  paired = false was given at construction time.
	 *  The returned reference is valid until the next call to advance(). */
	virtual const std::vector<BamTools::BamAlignment*>& getAlignments() const;

	/** Returns a reference to alignments of first read in current read pair.
	 *  The returned reference is valid until the next call to advance(). */
	virtual const std::vector<BamTools::BamAlignment*>& getAlignmentsFirst() const;

	/** Returns a reference to alignments of second read in current read pair.
	 *  The returned reference is valid until the next call to advance(). */
	virtual const std::vector<BamTools::BamAlignment*>& getAlignmentsSecond() const;

	/** Same as getAlignments, but transfers pointer ownerships. */
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignments();
	/** Same as getAlignmentsFirst, but transfers pointer ownerships. */
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignmentsFirst();
	/** Same as getAlignmentsSecond, but transfers pointer ownerships. */
	virtual std::auto_ptr<std::vector<BamTools::BamAlignment*> >  releaseAlignmentsSecond();

	virtual bool isSingleEnd() const;
	
	virtual bool isUnmapped() const;
	virtual bool isFirstUnmapped() const;
	virtual bool isSecondUnmapped() const;

	virtual bool hasMultipleMappingsButNoXA() const;

	virtual void enableProgressMessages(std::ostream& os, int frequency);

	virtual long long getSkippedDuplicates() const;

	virtual const BamTools::RefVector& getReferenceData() const { return bam_reader.GetReferenceData(); }

	virtual BamTools::SamHeader getHeader() const { return bam_reader.GetHeader(); }
};

#endif /* SORTEDBAMREADER_H_ */
