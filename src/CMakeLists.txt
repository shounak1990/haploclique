# Copyright 2012 Tobias Marschall
# 
# This file is part of HaploClique.
# 
# HaploClique is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HaploClique is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HaploClique.  If not, see <http://www.gnu.org/licenses/>.

add_executable ( alignmentparser
	alignmentparser.cpp
)

add_executable ( seeksingletons
	seeksingletons.cpp
)

add_executable (haploclique
	haploclique.cpp
	CliqueFinder.cpp
	Clique.cpp
	QuasispeciesEdgeCalculator.cpp
	GaussianEdgeCalculator.cpp
	CliqueWriter.cpp
	EdgeWriter.cpp
	AnyDistributionEdgeCalculator.cpp
	ReadSetZTester.cpp
	ReadSetGenericTester.cpp
	ReadGroupAwareEdgeCalculator.cpp
	ReadSetGroupWiseZTester.cpp
	ReadGroups.cpp
	AlignmentRecord.cpp
	ShortDnaSequence.cpp
	CoverageMonitor.cpp
	VariationCaller.cpp
	SingleTrackCoverageMonitor.cpp
	HistogramBasedDistribution.cpp
	Distributions.cpp
	BamHelper.cpp
	VariationIndex.cpp
	IndelLengthDistribution.cpp
	PositionSet.cpp
	Variation.cpp
	NamedDnaSequence.cpp
	Histogram.cpp
	SortedBamReader.cpp
	BamParser.cpp
)

target_link_libraries(haploclique ${Boost_LIBRARIES} BamTools)
target_link_libraries(alignmentparser ${Boost_LIBRARIES})
target_link_libraries(seeksingletons ${Boost_LIBRARIES})

set(ALL_EXE_TARGETS seeksingletons alignmentparser haploclique)

install(TARGETS ${ALL_EXE_TARGETS} RUNTIME DESTINATION bin)
