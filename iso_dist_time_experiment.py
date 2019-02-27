
import time, re, sys
from czno import getMinBarrierPathDijkstra, basePairList

def experiment(filename):
    with open(filename, "r") as f:
        for s in f:
            data = s.strip().split(" ")
            assert len(data) == 3
            assert len(data[0]) == len(data[1])
            assert len(data[0]) == len(data[2])
            assert re.fullmatch(r"^[ACGU]+$", data[0]) is not None
            assert re.fullmatch(r"^[(.)]+$", data[1]) is not None
            assert re.fullmatch(r"^[(.)]+$", data[2]) is not None

            start = time.time()
            pathway = getMinBarrierPathDijkstra(data[0], data[1], data[2])
            elapsed_time = time.time() - start

            assert pathway[0] == data[1]
            assert pathway[-1] == data[2]
            bp_set1 = set(basePairList(data[1]))
            bp_set2 = set(basePairList(data[2]))
            dist = len(bp_set1 ^ bp_set2)
            assert dist == len(pathway) - 1

            print(str(len(data[0]))+" "+str(elapsed_time))

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("error")
        sys.exit()
    experiment(sys.argv[1])