import math


def write_mult_and_mesh_lines(
    I, J, K, x1, x2, y1, y2, z1, z2, MBLKS_i, MBLKS_j, MBLKS_k, name
) -> int:
    IJK = [I, J, K]
    XB = [x1, x2, y1, y2, z1, z2]
    MBLKS = [MBLKS_i, MBLKS_j, MBLKS_k]
    M0_i = MBLKS_i - 2
    M0_j = MBLKS_j - 2
    M0_k = MBLKS_k - 2
    M0 = [M0_i, M0_j, M0_k]
    I_LOWER = M0[0] - MBLKS[0] + 2
    I_UPPER = I_LOWER + MBLKS[0] - 1
    J_LOWER = M0[1] - MBLKS[1] + 2
    J_UPPER = J_LOWER + MBLKS[1] - 1
    K_LOWER = M0[2] - MBLKS[2] + 2
    K_UPPER = K_LOWER + MBLKS[2] - 1

    DX = (XB[1] - XB[0]) / IJK[0]
    DY = (XB[3] - XB[2]) / IJK[1]
    DZ = (XB[5] - XB[4]) / IJK[2]

    NX = math.ceil(IJK[0] / MBLKS[0])
    NY = math.ceil(IJK[1] / MBLKS[1])
    NZ = math.ceil(IJK[2] / MBLKS[2])

    NX = int(NX)
    NY = int(NY)
    NZ = int(NZ)
    LX = NX * DX
    LY = NY * DY
    LZ = NZ * DZ

    XB0 = (XB[0], XB[0] + LX, XB[2], XB[2] + LY, XB[4], XB[4] + LZ)
    IJK_LOC = str((NX, NY, NZ))
    XB_LOC = str(XB0)

    line = MBLKS_i * MBLKS_j * MBLKS_k
    line = str(line)

    mult_line = (
        f'&MULT ID="{name}",'
        + "DX=%s, DY=%s, DZ=%s," % (LX, LY, LZ)
        + " I_UPPER=%s, J_UPPER=%s, K_UPPER=%s " % (I_UPPER, J_UPPER, K_UPPER)
        + " /"
    )
    mesh_line = (
        "&MESH IJK="
        + IJK_LOC[1:-1]
        + ", XB="
        + XB_LOC[1:-1]
        + f', MULT_ID="{name}" '
        + "/ "
        + line
        + " Mesh"
    )

    print(mult_line)
    print(mesh_line)

    return MBLKS_i * MBLKS_j * MBLKS_k
