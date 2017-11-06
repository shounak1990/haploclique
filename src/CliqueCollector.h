/* Copyright 2012 Tobias Marschall
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

#ifndef CLIQUECOLLECTOR_H_
#define CLIQUECOLLECTOR_H_

#include <deque>
#include <algorithm>

#include "Clique.h"
#include "LogWriter.h"

class CliqueCollector {
private:
    std::deque<AlignmentRecord*>* superReads;
    LogWriter* lw;
    unsigned int id;
public:
    CliqueCollector(LogWriter* lw) : lw(lw), id(0) {
        superReads = new std::deque<AlignmentRecord*>;
    };

    virtual ~CliqueCollector() {
        assert(superReads->empty());
        delete superReads;
    };

    void add(std::unique_ptr<Clique> clique,std::map<std::string,int>* appearanceMap) {
        assert(clique.get() != nullptr);

//        std::cerr << "Clique " << id << std::endl;
        std::unique_ptr<std::vector<const AlignmentRecord*>> alignments = clique->getAllAlignments();

        if (lw != nullptr) {
            std::list<unsigned int> cll;            
            for (const auto& a : *alignments) {
                cll.push_back(a->getID());
            }
            lw->reportClique(this->id, cll);
        }

        AlignmentRecord* al;

        if (alignments->size() > 1) {
            al = new AlignmentRecord(alignments, this->id++,appearanceMap);
        } else {
            al = new AlignmentRecord(*(alignments->front()));
            this->id++;
        }

        superReads->push_back(al);
    };

    std::deque<AlignmentRecord*>* finish()
    {
        auto retVal = superReads;

        auto comp = [](AlignmentRecord* r1, AlignmentRecord* r2) { return r1->getIntervalStart() < r2->getIntervalStart(); };

        sort(retVal->begin(), retVal->end(), comp);

        superReads = new std::deque<AlignmentRecord*>;
        return retVal;
    };
};

#endif /* CLIQUECOLLECTOR_H_ */
