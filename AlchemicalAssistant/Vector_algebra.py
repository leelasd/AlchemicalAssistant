import math
import numpy as np


def IsImproper(ty):
    dihTYPS = [160, 161, 162, 221, 277]
    if ty in dihTYPS:
        return True
    else:
        return False


def pairing_func(a, b):
    ans = (a + b) * (a + b + 1) * 0.5
    if a > b:
        ans = ans + a
    else:
        ans = ans + b
    return (int(ans))


def Vector(x, y, z):
    return (x, y, z)


def length(v):
    "Return length of a vector."
    sum = 0.0
    for c in v:
        sum += c * c
    return math.sqrt(sum)


def subtract(u, v):
    "Return difference between two vectors."
    x = u[0] - v[0]
    y = u[1] - v[1]
    z = u[2] - v[2]
    return Vector(x, y, z)


def dot(u, v):
    "Return dot product of two vectors."
    sum = 0.0
    for cu, cv in zip(u, v):
        sum += cu * cv
    return sum


def Distance(u, v):
    "Return length of a vector."
#    print(u,v)
    uv = subtract(u, v)
    lsum = 0.0
    for c in uv:
        lsum += c * c
    return math.sqrt(lsum)


def cross(u, v):
    "Return the cross product of two vectors."
    x = u[1] * v[2] - u[2] * v[1]
    y = u[2] * v[0] - u[0] * v[2]
    z = u[0] * v[1] - u[1] * v[0]
    return Vector(x, y, z)


def Mol_angle(v0, v1):
    "Return angle [0..pi] between two vectors."
    cosa = round(dot(v0, v1) / length(v0) / length(v1), 3)
    return np.arccos(cosa)


def angle(p0, p1, p2):
    "Return angle [0..pi] between two vectors."
    v0 = subtract(p0, p1)
    v1 = subtract(p2, p1)
    cosa = dot(v0, v1) / length(v0) / length(v1)
#    print(cosa)
    return 180.0 * np.arccos(round(cosa, 3)) * 7.0 / 22.0


def dihedral(p0, p1, p2, p3):
    "Return angle [0..2*pi] formed by vertices p0-p1-p2-p3."
    v01 = subtract(p0, p1)
    v32 = subtract(p3, p2)
    v12 = subtract(p1, p2)
    v0 = cross(v12, v01)
    v3 = cross(v12, v32)
    # The cross product vectors are both normal to the axis
    # vector v12, so the angle between them is the dihedral
    # angle that we are looking for.  However, since "angle"
    # only returns values between 0 and pi, we need to make
    # sure we get the right sign relative to the rotation axis
    a = Mol_angle(v0, v3)
    if dot(cross(v0, v3), v12) > 0:
        a = -a
    return a * 180.0 * 7.0 / 22.0


def ang_id(a):
    bond_a = pairing_func(a[0], a[1])
    bond_b = pairing_func(a[1], a[2])
    return(pairing_func(bond_a, bond_b))


def AtomNum2Mass(num):
    dda = {
        1:	1.00800,
        5:	10.8110,
        6:	12.0110,
        7:	14.0070,
        8:	15.9990,
        9:	18.9980,
        11:	22.9900,
        12:	24.3050,
        14:	28.0860,
        15:	30.9740,
        16:	32.0650,
        17:	35.4530,
        19:	39.0980,
        20:	40.0780,
        35:	79.9040,
        53:  126.9050,
    }
    return(dda[int(num)])


def AtomNum2Symb(num):
    dda = {
        1	: 'H',
        5	: 'B',
        6	: 'C',
        7	: 'N',
        8	: 'O',
        9	: 'F',
        11	: 'Na',
        12	: 'Mg',
        14	: 'Si',
        15	: 'P',
        16	: 'S',
        17	: 'Cl',
        19	: 'K',
        20	: 'Ca',
        35	: 'Br',
        53	: 'I', }
    return(dda[int(num)])
