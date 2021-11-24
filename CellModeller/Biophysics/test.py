import SPP


bac = SPP.SPP(dt_c=0.01, z_axis_c=False)
for i in range(10):
    bac.add_cell(i, i, 0.3)
for i in range(10):
    bac.step()
print(len(bac.cell_centers))
print(type(bac.cell_centers))
