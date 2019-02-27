
import glob
import tempfile
import subprocess
import re

def GetS151Dataset():
    dataset = []
    for address in glob.glob("/home/czno-gen/partial_bpseq/*.bpseq"):
        with open(address, "r") as f:
            data = []
            for line in f:
                data.append(str(line).strip().split(" "))
            sequence = ""
            structure = ""
            for d in data:
                sequence += str(d[1])
                if int(d[2]) == 0:
                    structure += "."
                elif int(d[2]) >int(d[0]):
                    structure += "("
                else:
                    structure += ")"
            dataset.append([address, sequence, structure])
    return dataset



def GetLocalOptima(sequence):
    result = []

    with open(".tmp1",'w') as file:
        file.write(sequence)
    subprocess.call('RNAsubopt -p 10000 < "%s" > "%s"' % (".tmp1", ".tmp2"), shell=True)
    subprocess.call('sed -e "1d" %s > "%s"' % (".tmp2", ".tmp3"), shell=True)
    subprocess.call('RNAlocmin -s "%s" < "%s" > "%s"' % (".tmp1", ".tmp3", ".tmp4"), shell=True)
    with open(".tmp4", "r") as f:
        for line in f:
            words = line.strip().split(" ")
            if len(words) < 2: continue
            if re.fullmatch(r"^[(.)]+$", words[1]) is None: continue
            if len(words[1]) != len(sequence): continue
            result.append(words[1])

    return result

if __name__ == '__main__':
    aaa = GetS151Dataset()
    bbb = GetLocalOptima(aaa[0][1])
    print(bbb)
    