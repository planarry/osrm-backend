/*

Copyright (c) 2013, Project OSRM, Dennis Luxen, others
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef IMPORT_EDGE_H
#define IMPORT_EDGE_H

#include "../typedefs.h"

struct NodeBasedEdge
{
    bool operator<(const NodeBasedEdge &e) const;

    explicit NodeBasedEdge(NodeID source,
                           NodeID target,
                           NodeID name_id,
                           EdgeWeight weight,
                           bool forward,
                           bool backward,
                           short type,
                           bool roundabout,
                           bool in_tiny_cc,
                           bool access_restricted,
                           bool contra_flow,
                           bool is_split,
                           int length,
                           short maxload,
                           short maxheight);

    NodeID source;
    NodeID target;
    NodeID name_id;
    EdgeWeight weight;
    short type;
    bool forward : 1;
    bool backward : 1;
    bool roundabout : 1;
    bool in_tiny_cc : 1;
    bool access_restricted : 1;
    bool contra_flow : 1;
    bool is_split : 1;
    int length;
    short maxload;
    short maxheight;

    NodeBasedEdge() = delete;
};

struct EdgeBasedEdge
{

  public:
    bool operator<(const EdgeBasedEdge &e) const;

    template <class EdgeT> explicit EdgeBasedEdge(const EdgeT &myEdge);

    EdgeBasedEdge();

    explicit EdgeBasedEdge(const NodeID source,
                           const NodeID target,
                           const NodeID edge_id,
                           const EdgeWeight weight,
                           const bool forward,
                           const bool backward,
                           int length,
                           short maxload,
                           short maxheight);
    NodeID source;
    NodeID target;
    NodeID edge_id;
    EdgeWeight weight : 30;
    bool forward : 1;
    bool backward : 1;
    int length;
    short maxload;
    short maxheight;
};

typedef NodeBasedEdge ImportEdge;

#endif /* IMPORT_EDGE_H */
