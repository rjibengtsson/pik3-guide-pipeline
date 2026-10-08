from __future__ import annotations

import re
from collections import defaultdict
from typing import Iterable, Iterator, Mapping, NamedTuple

import pysam



class Hit(NamedTuple):
    """One alignment of a tile to a reference transcript."""
    query_id: str
    reference: str
    nm: int
    is_reverse: bool = False

def hits_from_sam(path: str, query_name: str) -> Iterator[Hit]:

    with pysam.AlignmentFile(path, "r") as samfile:
        for read in samfile:
            if read.query_name == query_name:
                yield Hit(
                    query_id=read.query_name,
                    reference=read.reference_name,
                    nm=read.get_tag("NM"),
                    is_reverse=read.is_reverse
                )





def main():
    query_name = "NM_005026.5:2754-2783"
    hits = list(hits_from_sam("classify-guides-wjhtb7_3/tiles.sam", query_name))
    print(hits)



if __name__ == "__main__":
    main()