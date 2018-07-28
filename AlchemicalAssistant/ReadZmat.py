# FUNCTION  TO RETURN 0 or 1 FOR A MATCHED PATTERN
def new_func(linex):
    out = 0
    for word in linex.split():
        if(word == "Non-Bonded"):
            out = out + 1
    return(out)


def read_nbond(data, nline, oline):
    QLJDat = {}
    num = 0
    for i in range(oline, nline):
        if(len(data[i].split()) == 6):
            QLJDat[num] = data[i]
            num = num + 1
    return(QLJDat)
# SAVED THE LINE NUMBER AND ALL THE DATA HERE


def read_files(infile):
    nline = 0
    oline = 0
    data = []
    for line in infile:
        data.append(line)
        if(new_func(line) == 1):
            oline = nline
        nline += 1
    return data, nline, oline
# SAVED THE LINE NUMBER AND ALL THE DATA HERE


def GetZmatQLJDat(zmat):
    qfile = open(zmat)
    qdat, nl1, ol1 = read_files(qfile)
    QLJD = read_nbond(qdat, nl1, ol1)
    return(QLJD)


def AddImpLinestoZmat(zmat, imp_lines):
    zlines = open(zmat, 'r').readlines()
    dom_line = 0
    for no in range(len(zlines)):
        if 'Domain Definitions follow' in zlines[no]:
            dom_line = no
            break
    final_zlines = zlines[:dom_line] + imp_lines + zlines[dom_line:]
    of = open(zmat, 'w')
    for line in final_zlines:
        of.write('%s' % line)
    of.close()
#    print(final_zlines)
    return None
