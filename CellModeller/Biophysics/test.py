import SPP


bac = SPP.SPP(dt_c=0.01, z_axis_c=False)
for i in range(10):
    bac.addCell(i, i, 0, 1, 0, 0)
for i in range(10):
    bac.step()
print(len(bac.cell_centers))
print(type(bac.cell_centers))
