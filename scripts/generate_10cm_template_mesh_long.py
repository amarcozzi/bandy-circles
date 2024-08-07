from mult_writer import write_mult_and_mesh_lines

print()

"""
m1
"""
name = "m1"
x1 = -10
x2 = 10
dx = 0.1
i = int((x2 - x1) / dx)
num_x_cells_per_block = 50
MBLKS_i = i // num_x_cells_per_block

y1 = -20
y2 = 20
dy = 0.1
j = int((y2 - y1) / dy)
num_y_cells_per_block = 50
MBLKS_j = j // num_y_cells_per_block


z1 = 0
z2 = 5
dz = 0.1
k = int((z2 - z1) / dz)
num_z_cells_per_block = 50
MBLKS_k = k // num_z_cells_per_block

write_mult_and_mesh_lines(
    i, j, k, x1, x2, y1, y2, z1, z2, MBLKS_i, MBLKS_j, MBLKS_k, name
)
print()

"""
m2
"""
name = "m2"
x1 = -10
x2 = 10
dx = 0.2
i = int((x2 - x1) / dx)
num_x_cells_per_block = 50
MBLKS_i = i // num_x_cells_per_block

y1 = 20
y2 = 60
dy = 0.2
j = int((y2 - y1) / dy)
num_y_cells_per_block = 50
MBLKS_j = j // num_y_cells_per_block


z1 = 0
z2 = 5
dz = 0.2
k = int((z2 - z1) / dz)
num_z_cells_per_block = 25
MBLKS_k = k // num_z_cells_per_block

write_mult_and_mesh_lines(
    i, j, k, x1, x2, y1, y2, z1, z2, MBLKS_i, MBLKS_j, MBLKS_k, name
)
print()


"""
m3
"""
name = "m3"
x1 = -10
x2 = 10
dx = 0.2
i = int((x2 - x1) / dx)
num_x_cells_per_block = 50
MBLKS_i = i // num_x_cells_per_block

y1 = -20
y2 = 60
dy = 0.2
j = int((y2 - y1) / dy)
num_y_cells_per_block = 50
MBLKS_j = j // num_y_cells_per_block


z1 = 5
z2 = 25
dz = 0.2
k = int((z2 - z1) / dz)
num_z_cells_per_block = 50
MBLKS_k = k // num_z_cells_per_block

write_mult_and_mesh_lines(
    i, j, k, x1, x2, y1, y2, z1, z2, MBLKS_i, MBLKS_j, MBLKS_k, name
)
print()

"""
m4
"""
name = "m4"
x1 = -10
x2 = 10
dx = 0.4
i = int((x2 - x1) / dx)
num_x_cells_per_block = 50
MBLKS_i = i // num_x_cells_per_block

y1 = -20
y2 = 60
dy = 0.4
j = int((y2 - y1) / dy)
num_y_cells_per_block = 50
MBLKS_j = j // num_y_cells_per_block


z1 = 25
z2 = 45
dz = 0.4
k = int((z2 - z1) / dz)
num_z_cells_per_block = 50
MBLKS_k = k // num_z_cells_per_block

write_mult_and_mesh_lines(
    i, j, k, x1, x2, y1, y2, z1, z2, MBLKS_i, MBLKS_j, MBLKS_k, name
)
print()

# """
# m5
# """
# name = "m5"
# x1 = -10
# x2 = 10
# dx = 0.8
# i = int((x2 - x1) / dx)
# num_x_cells_per_block = 25
# MBLKS_i = i // num_x_cells_per_block

# y1 = -20
# y2 = 40
# dy = 0.8
# j = int((y2 - y1) / dy)
# num_y_cells_per_block = 50
# MBLKS_j = j // num_y_cells_per_block


# z1 = 45
# z2 = 85
# dz = 0.8
# k = int((z2 - z1) / dz)
# num_z_cells_per_block = 50
# MBLKS_k = k // num_z_cells_per_block

# write_mult_and_mesh_lines(
#     i, j, k, x1, x2, y1, y2, z1, z2, MBLKS_i, MBLKS_j, MBLKS_k, name
# )
# print()