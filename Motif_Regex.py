import re
import pandas as pd
import numpy as np

#test sequences
ZP = "MLSCNICGETVTSEPDMKAHLIVHMESEIICPFCKLSGVNYDEMCFHIETAHFEQNTLERNFERINTVQYGTSDNKKDNTLQCGMEVNSSILSGCASNHPKNSAQNLTKD" \
     "STLKHEGFYSENLTESRKFLKSREKQSSLTEIKGSVYETTYSPPECPFCGKIEEHSEDMETHVKTKHANLLDIPLEDCDQPLYDCPMCGLICTNYHILQEHVDLHLEENS" \
     "FQQGMDRVQCSGDLQLAHQLQQEEDRKRRSEESRQEIEEFQKLQRQYGLDNSGGYKQQQLRNMEIEVNRGRMPPSEFHRRKADMMESLALGFDDGKTKTSGIIEALHRYY" \
     "QNAATDVRRVWLSSVVDHFHSSLGDKGWGCGYRNFQMLLSSLLQNDAYNDCLKGMLIPCIPKIQSMIEDAWKEGFDPQGASQLNNRLQGTKAWIGACEVYILLTSLRVKC" \
     "HIVDFHKSTGPLGTHPRLFEWILNYYSSEGEGSPKVVCTSKPPIYLQHQGHSRTVIGIEEKKNRTLCLLILDPGCPSREMQKLLKQDIEASSLKQLRKSMGNLKHKQYQI" \
     "LAVEGALSLEEKLARRQASQVFTAEKIP"
Rap = "MPRRKKKVKEVSESRNLEKKDVETTSSVSVKRKRRLEDAFIVISDSDGEEPKEENGLQKTKTKQSNRAKCLAKRKIAQMTEEEQFALALKMSEQEAREVNSQEEEEEEL" \
      "LRKAIAESLNSCRPSDASATRSRPLATGPSSQSHQEKTTDSGLTEGIWQLVPPSLFKGSHISQGNEAEEREEPWDHTEKTEEEPVSGSSGSWDQSSQPVFENVNVKSF" \
      "DRCTGHSAEHTQCGKPQESTGRG"

Input = pd.read_csv("uniprot-9606.txt", sep="\t") # input dataframe to search.

'''Index(['Entry', 'Entry name', 'Status', 'Protein names', 'Gene names',
       'Organism', 'Length', 'Sequence', 'Gene ontology (GO)', 'Domain [CC]',
       'Domain [FT]'],
      dtype='object')'''

#Phi = V, I, L, F, W, Y and M.

#MIU CONSENSUS:XXDφXLAXXLXX###X
#This first regex is for the more canonical LAXXL ,motif
#MIU = re.compile('...D(V|I|L|F|W|Y|M).LA..L..(E|D|Q)(E|D|Q)(E|D|Q).')
#This regex requies hydrophobic residues in the LAXXL motif.

MIU = re.compile('...D[VILFWYM].[VILFWYM]A..[VILFWYM]..[EQD]{3}.')
#ZHA = ExxxFxxLxxxY
#ZHA = re.compile('([ED])...([VILFWYM])..([VILFWYM])...([VILFWYM])')
ZHA = re.compile('E...F..L...Y')
UIM = re.compile('[VILFWYM]..A[VILFWYM]..S')

#for i in re.findall(MIU, ZP): #regex testing
    #print(i)

def domainsearch(x, y, z):   #x is the regex, y is the col to enter (as a string)
    for i, r in z.iterrows():
        #print(r["Sequence"])

        list_search = []
        for w in re.findall(x, str(r["Sequence"])):
            if w is None:
                pass
            else:
                # print(w)
                list_search.append(w)

                df_list_search = []
                for l in list_search:
                    for m in re.finditer("{}".format(l), str(r["Sequence"])):
                        # print(m.start())
                        # print(m.group())
                        # print(m.end())
                        df_list_search.append([m.start(), m.group(), m.end()])

                print(df_list_search)
                z.loc[i, y] = str(df_list_search)


domainsearch(MIU, "MIU", Input)
domainsearch(ZHA, "ZHA", Input)
domainsearch(UIM, "UIM", Input)

Input.to_csv("motif_regex_output.csv")






